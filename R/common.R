if (!exists("mreadRDS")) {mreadRDS = memoise::memoise(readRDS)}
  
ewas_func2 = function(d, e, USE_PARAPPLY, model_formula, nb_fact_of_interest=1, model_func_name="modelcalllm") {
  model_func = get(model_func_name)
  y_name = rownames(attr(stats::terms(as.formula(model_formula)), "factor"))[1]
  x_values = e[colnames(d),c(rownames(attr(stats::terms(as.formula(model_formula)), "factor"))[-1], colnames(e)[1])]
  head(x_values)

  set.seed(1)
  tmp_d = x_values
  tmp_d[[y_name]] = rnorm(nrow(tmp_d))
  fake_m = lm(model_formula, tmp_d)
  expected_nb_coef = length(fake_m$coefficients)
  expected_nb_na = sum(is.na(fake_m$coefficients))

  if (USE_PARAPPLY) {
    print("ewas using parallel::parApply...")
    if (!exists("cl")) {
      nb_cores = parallel::detectCores()
      cl <<- parallel::makeCluster(nb_cores,  type="FORK")
      # parallel::stopCluster(cl)
    }
    ewas = parallel::parApply(cl, d, 1, model_func, x_values=x_values, model_formula=model_formula, nb_fact_of_interest=nb_fact_of_interest, expected_nb_coef=expected_nb_coef, expected_nb_na=expected_nb_na)
  } else {
    print("ewas using epimedtools::monitored_apply")
    ewas = epimedtools::monitored_apply(d, 1, model_func, x_values=x_values, model_formula=model_formula, nb_fact_of_interest=nb_fact_of_interest, expected_nb_coef=expected_nb_coef, expected_nb_na=expected_nb_na)
  }
  print("done")
  ewas = t(ewas)
  # head(ewas)
  # dim(ewas)
  return(ewas)
}
if (!exists("mewas_func2")) mewas_func2 = memoise::memoise(ewas_func2)
  
  
modelcalllm = function(meth, x_values, model_formula, nb_fact_of_interest, expected_nb_coef, expected_nb_na=0) {
  # meth = d["cg07164639",]
  # meth = d[140340,]

  # model
  options(contrasts=c("contr.sum", "contr.poly"))
  # options(contrasts=c("contr.treatment", "contr.poly"))

  y_name = rownames(attr(stats::terms(as.formula(model_formula)), "factor"))[1]
  tmp_d = x_values
  tmp_d[[y_name]] = meth
  head(tmp_d)
  m = try(lm(model_formula, tmp_d))
  if (attributes(m)$class == "try-error") {
    # warning(paste0("Region ", mdata$region_id, " contains probe ", probe, "that through error model. Have to fix it! "))
    set.seed(1)
    tmp_d[[y_name]] = rnorm(nrow(tmp_d))
    m = lm(model_formula, tmp_d)
    a = anova(m)
    m$coefficients[] = NA
    m$residuals[] = NA
    fstat = summary(m)$fstatistic
  } else if (length(m$coefficients) != expected_nb_coef | sum(is.na(m$coefficients)) > expected_nb_na) {
    # warning(paste0("Region ", mdata$region_id, " contains probe ", probe, "that through error model. Have to fix it! "))
    set.seed(1)
    tmp_d[[y_name]] = rnorm(nrow(tmp_d))
    m = lm(model_formula, tmp_d)
    a = anova(m)
    summary(m)$fstatistic
    m$coefficients[] = NA
    m$residuals[] = NA
    fstat = summary(m)$fstatistic
  } else {
    a = anova(m)
    fs1 = sum(a[1:nb_fact_of_interest,2]) / sum(a[1:nb_fact_of_interest,1]) / a[nrow(a),3]
    fs2 = sum(a[1:nb_fact_of_interest,1])
    # fs1 = sum(a[1,2]) / sum(a[1,1]) / a[nrow(a),3]
    # fs2 = sum(a[1,1])
    fs3 = a[nrow(a),1]
    fstat = c(fs1, fs2, fs3)
  }

  # pval for ewas
  if (is.null(fstat)) {
    # warning(paste0("Region ", mdata$region_id, " contains probe ", probe, "that through NULL fstat. Have to fix it! "))
    m$coefficients[] = NA
    lpval_fisher = NA
  } else {
    lpval_fisher = pf(fstat[[1]], fstat[[2]], fstat[[3]], lower.tail=FALSE, log.p=TRUE)/-log(10)
  }

  # ret = c(lpval_fisher=lpval_fisher, r2=summary(m)$r.squared,  nb_notna = sum(!is.na(meth)), fstat_1=fstat[[1]], fstat_2=fstat[[2]], fstat_3=fstat[[3]]) #m$coefficients[-1])
  # return(ret)
  lpv = lpval_fisher
  if (is.factor(x_values[,1])) {
    beta = -2*m$coefficients[[2]]
  } else if (is.numeric(x_values[,1])) {
    beta = m$coefficients[[2]]
  } else {
    stop("Phenotype in not a factor nor numeric (in model_call_lm).")
  }
  # coef = c(unlist(m$coefficients[2:length(levels(pheno))]))
  # coef = c(coef, pheno=-sum(coef))
  ret = c(beta=beta, lpv=lpv, coef=unlist(m$coefficients), pv=a[1:(nrow(a)-1),5])
  ret
  return(ret)
}

plot_res = function(
  roi,
  idx_probes,
  legendplace="topright",
  combp_res_probes,
  ewas,
  study_filename="",
  pheno_key,
  ADD_LEGEND=FALSE,
  rename_results=I,
  pheno_key2,
  LAYOUT=TRUE
) {
  if (!missing(study_filename)) {
    s = mreadRDS(study_filename)
  }
  d = s$data
  e = s$exp_grp

  if (missing(pheno_key2)) {
    pheno_key2 = "hb"
  }
  if (!pheno_key2 %in% colnames(e)) {
    pheno_key2 = pheno_key
  }

  if (missing(idx_probes)) {
    idx_probes = ewas[as.character(ewas[,1])==as.character(roi[[1]]) & ewas[,2]>=roi[[2]] & ewas[,2]<=roi[[3]],4]
  }

  # if (missing(pheno_key2)) {
  #   pheno_key2 = pheno_key
  # }


  idx_probes = intersect(idx_probes, rownames(d))
  d = d[idx_probes,]
  if ("tissue" %in% colnames(e)) {
    idx_sample = rownames(e)[order(e[["tissue"]], e[[pheno_key]], e[[pheno_key2]])]
  } else if ("hb" %in% colnames(e)) {
    idx_sample = rownames(e)[order(e[[pheno_key]], e[["hb"]])]
  } else {
    idx_sample = rownames(e)[order(e[[pheno_key]])]
  }

  if (length(idx_probes) > 1) {
    # layout(matrix(c(1, 2, 2, 2, 2), 1))
    if (LAYOUT) {
      layout(matrix(c(
        c(1,2,2,2,2,4,4,4),
        c(1,2,2,2,2,3,3,3),
        c(1,2,2,2,2,3,3,3),
        c(1,2,2,2,2,3,3,3)
        ), 4, byrow=TRUE), respect=TRUE)
    }

    # PHENO
    par(mar=c(5.7, 4.1, 4.1, 0))
    if (is.factor(e[,pheno_key])) {
      col=(1:4)[as.numeric(e[idx_sample,pheno_key])]
    } else {
      col=1
    }
    plot(
      e[idx_sample,pheno_key2], 1:length(idx_sample),
      xlab=rename_results(pheno_key2), ylab=rename_results(pheno_key),
      yaxt="n", yaxs = "i", main="",
      col=col,
      pch=16
    )
    # axis(1, at=c(min(e[idx_sample,pheno_key2]), max(e[idx_sample,pheno_key2])), signif(c(min(e[idx_sample,pheno_key2], na=rm=TRUE), max(e[idx_sample,pheno_key2])), 3))

    # METH
    colors = c("cyan", "black", "red")
    cols = colorRampPalette(colors)(1000)
    breaks = seq(0, 1, length.out = length(cols) + 1)
    main = paste0(rename_results(study_filename), " ", roi[[1]], ":", roi[[2]], "-", roi[[3]])
    par(mar=c(5.7, 0, 4.1,0))
    image(d[idx_probes,idx_sample], col=cols, breaks=breaks, xaxt="n", yaxt="n", main=main)
    axis(1, (1:nrow(d[idx_probes,idx_sample]) - 1)/(nrow(d[idx_probes,idx_sample]) - 1), rownames(d[idx_probes,idx_sample]), las = 2)

    # DATAFRAME
    df = lapply(idx_sample, function(i){
      df = lapply(idx_probes, function(j){
        # if (is.null(confounder)) {
          list(meth=d[j,i], probes=as.character(j), pheno=e[i,pheno_key])
        # } else {
        #   list(meth=d[j,i], probes=as.character(j), pheno=e[i,pheno_key], confounder=e[i,confounder])
        # }
      })
      df = do.call(rbind, df)
      df
    })
    df = do.call(rbind, df)
    df = data.frame(lapply(data.frame(df, stringsAsFactors=FALSE), unlist), stringsAsFactors=FALSE)
    # df
    df$probes = factor(df$probes, levels=idx_probes)
    if (is.factor(df$pheno)) {
      df$pheno = factor(df$pheno, levels=levels(df$pheno)[levels(df$pheno)%in%unique(as.character(df$pheno))])
    }

    # # effect
    # options(contrasts=c("contr.treatment", "contr.poly"))
    # if (is.null(confounder)) {
    #   m = lm(meth~pheno+probes, df)
    # } else {
    #   m = lm(meth~pheno+probes+confounder, df)
    # }
    # m$coefficients
    # main = paste0(levels(df$pheno)[2], " effect = " ,   signif(m$coefficients[[2]],3))
    main = ""

    # BOXPLOT
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    if (is.factor(e[[pheno_key]])) {
      par(mar=c(5.7, 4.1, 0, 0))
      # if (is.null(confounder)) {
        boxplot(meth~pheno+probes, df, las=2, col=1:length(unique(na.omit(df$pheno))), ylim=c(0,1), cex.axis=0.5,
          ylab="methylation",
          #, yaxt="n",
          cex.axis=0.5
        )
        # axis(4)
        # mtext("methylation", side=4, line=3)

      # } else {
        # boxplot(meth~pheno+confounder+probes, df, las=2, col=1:length(unique(na.omit(df$pheno))), ylim=c(0,1), main=main)
      # }
    } else {
      plot(0,type='n',axes=FALSE,ann=FALSE)
    }

    # EWAS
    if (!missing(ewas)) {
      # combp
      if (!missing(combp_res_probes)) {
        idx_probes = ewas[as.character(ewas[,1])==as.character(roi[[1]]) & ewas[,2]>=roi[[2]] & ewas[,2]<=roi[[3]],4]
        sub_ewas = ewas[ewas[,4]%in%idx_probes, ]
        sub_ewas = sub_ewas[!duplicated(paste0(sub_ewas[,1], ":", sub_ewas[,2])), ]
        rownames(sub_ewas) = paste0(sub_ewas[,1], ":", sub_ewas[,2])
        head(sub_ewas)
        dim(sub_ewas)
        pval_ewas = combp_res_probes[paste0(combp_res_probes[,1], ":", combp_res_probes[,2]) %in% rownames(sub_ewas),4]
        pval_slk =  combp_res_probes[paste0(combp_res_probes[,1], ":", combp_res_probes[,2]) %in% rownames(sub_ewas),5]
        qval_slk =  combp_res_probes[paste0(combp_res_probes[,1], ":", combp_res_probes[,2]) %in% rownames(sub_ewas),6]
        pval_ewas[pval_ewas==0] = 10^-45
        pval_slk [pval_slk ==0] = 10^-45
        qval_slk [qval_slk ==0] = 10^-45
      } else {
        pval_ewas = 10^-ewas[rownames(ewas) %in% idx_probes, 2]
        # pval_ewas[pval_ewas==0] = 10^-45
        qval_slk = pval_slk = pval_ewas
      }

      # layout(matrix(c(2,1,1,1,1), 1))
      par(mar=c(5.1, 2.1, 0, 0))
      par(mar=c(0, 4.1, 4.1, 0))
      x = 1:length(-log10(pval_ewas))
      plot(x, -log10(pval_ewas), col="red", xaxt="n",
        xlab="", ylab="-log10(pv)",
        # yaxt="n",
        # main=paste0("meth~", gene),
        ylim=c(0, min(45, max(-log10(pval_slk), -log10(pval_ewas)))),
        type="l", lty=3
      )
      # axis(4)
      # mtext("-log10(pv)", side=4, line=3)
      axis(3, at=x, labels=names(pval_ewas),las=2, cex.axis = 0.5)
      lines(-log10(pval_slk), col="blue"  , type="l", lty=3)
      lines(-log10(qval_slk), col="purple", type="l", lty=3)

      # # add Student pvals
      # if (length(gene_symbols)>1) {
      #   for (g in gene_symbols) {
      #     lines(sub_ewas[,paste0("lpval_student_", g)], col=pals::glasbey()[which(gene_symbols%in%g)], type="l")
      #   }
      # }
      # # add DMR
      # abline(h=-log10(as.numeric(pval_tresh)), col="black", lwd=1, lty=2)
      # for (i in 1:nrow(combp_res_region)) {
      #   x1 = c(which(sub_ewas[,2] == combp_res_region[i,2]), which(sub_ewas[,3] == combp_res_region[i,3]))
      #   y1 = c(-log10(as.numeric(pval_tresh)), -log10(as.numeric(pval_tresh)))
      #   lines(x1,y1, type="o", col="green", pch=18, lwd=4)
      # }
      # add legend

      if (ADD_LEGEND) {
        col = c("red","blue", "purple", "black", "green")
        lwd = c(1,1,1,1,4)
        lty = c(3,3,3,2,1)
        legend=c("pval Fisher", "pval SLK", "qval SLK",  "threshold", "DMR")
        legend(legendplace, legend=legend, col=col, lwd=lwd, lty=lty)
        par(mar=c(0,0,0,0), mgp=c(3, 1, 0), las=0)
        plot.new()
        par(mar=c(0,0,0,0), mgp=c(3, 1, 0), las=0)
      }
    } else {
      plot(0,type='n',axes=FALSE,ann=FALSE)
    }


    # return(mdata)
  } else {
    par(mar=c(0, 0, 0, 0), mgp=c(3, 1, 0), las=0)
    plot(0,type='n',axes=FALSE,ann=FALSE)
    plot(0,type='n',axes=FALSE,ann=FALSE)
    plot(0,type='n',axes=FALSE,ann=FALSE)
    plot(0,type='n',axes=FALSE,ann=FALSE)
  }
  par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
}





# if (!exists("mreadtablegz")) {mreadtablegz = memoise::memoise(function(file, ...){read.table(gzfile(file), ...)})}
# if (!exists("mread.table")) {mread.table = memoise::memoise(function(file, ...){read.table(file, ...)})}
# if (!exists("mread.xlsx")) {mread.xlsx = memoise::memoise(openxlsx::read.xlsx)}
#
#
#
#
#
#
#
#
#
#
#
#
#
# # library(RefFreeEWAS)
# call_rfe = function(Dsub, k) {
#   # load the packages you want
#   # exectime_tmp <- system.time({
#   #   foo = RefFreeCellMix(Dsub,k=5,iters=9)
#   #   Test = foo$Mu
#   #   Aest = t(foo$Omega)
#   #   Dest <- Test %*% Aest
#   #   better with this code
#     foo = RefFreeEWAS::RefFreeCellMixArray(Dsub,K=k,iters=9,dist.method="manhattan")
#     Aest = t(foo[[as.character(k)]]$Omega)
#     Test = foo[[as.character(k)]]$Mu
#   # })
#
#   # exectime = exectime_tmp[3]
#   # Dest = D
#   # Dest[] = 0
#   # dir.create(res_dir, showWarnings=FALSE, recursive=TRUE)
#   # save(Aest, Dest, exectime, file = paste0(res_dir, "/res.RData"))
#   # layout(1, respect=TRUE)
#   # proportion.heatmap(Aest, clusterCols=TRUE, clusterRows=TRUE)
#   return(Aest)
# }
#
#
#
#
#
#
#
#
#
# rename_results = function(vec) {
#   vec = gsub("dmr_"                , ""        ,vec)
#   vec = gsub("modelcalllm_"        , ""        ,vec)
#   vec = gsub(".regions-t.bed"      , ""        ,vec)
#   vec = gsub("meth~cmsnocmsrinc"   , "CMS"        ,vec)
#   vec = gsub("meth~citynocms"           , "alti"        ,vec)
#   vec = gsub("\\+tissue"             , ""        ,vec)
#   vec = gsub("_0.00001"            , ""        ,vec)
#   vec = gsub("_0.00005"            , ""        ,vec)
#   vec = gsub("_0.0001"            , ""        ,vec)
#   vec = gsub("_0.0005"            , ""        ,vec)
#   vec = gsub("_0.001"            , ""        ,vec)
#   vec = gsub("_0.005"            , ""        ,vec)
#   vec = gsub("_0.01"            , ""        ,vec)
#   vec = gsub("_0.05"            , ""        ,vec)
#   vec = gsub("study_exp5300_mergemeth_" , "BS"        ,vec)
#   vec = gsub("study_exp5300_meth_" , ""        ,vec)
#   vec = gsub(".rds"                , ""        ,vec)
#   vec = gsub("n_probes"            , "n"       ,vec)
#   vec = gsub("z_sidak_p"           , "p"       ,vec)
#   vec = gsub("cmsnolima"          , "cnl"     ,vec)
#   vec = gsub("cmslimactrl"        , "CMS"     ,vec)
#   vec = gsub("citynocms"          , "city"    ,vec)
#   vec = gsub("convoluted"          , "cnv"     ,vec)
#     vec = gsub("cnv"               , ""        ,vec)
#   vec = gsub("rnbeads"             , "rnb"     ,vec)
#   vec = gsub("rnb_rfe"             , "rnb_rfe" ,vec)
#   vec = gsub("corrected"           , "dec"     ,vec)
#   vec = gsub("corrected"           , "dec"     ,vec)
#   vec = gsub("saliva_"             , "saliva"  ,vec)
#   vec = gsub("blood_"              , "blood"   ,vec)
#   vec = gsub("_hb_"                , "_hmb_"   ,vec)
#   vec
#
# }
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# proportion.heatmap<-function(
#     Ahat,
#     sample.characteristic=NULL,
#     clusterCols=FALSE,
#     clusterRows=FALSE
#     ){
#
#   if(!is.null(sample.characteristic)){
#     data.ch<-sample.characteristic
#
#     if(is.numeric(data.ch)){
#       palette<-grey.colors(start=0.9, end=0.1, min(20, length(data.ch)))
#       data.cols<-palette[cut(data.ch,min(20, length(data.ch)),labels=FALSE)]
#       data.cols[is.na(data.cols)]<-"#ffffff"
#     }else{
#       if(is.character(data.ch)){
#         data.ch<-factor(data.ch, levels=unique(data.ch))
#       }
#       palette<-rainbow(length(levels(data.ch)))
#       data.cols<-palette[as.integer(data.ch)]
#     }
#   }else{
#     data.cols<-NA
#   }
#
#   N_COL_BINS<-50
#   hm_args<-list(x=Ahat,
#       scale="none", trace="none", density.info="none",
#       keysize=1, key.xlab="proportion",  key.title="proportion",
#       Colv=clusterCols, Rowv=clusterRows,
#       dendrogram=c("none", "column", "row", "both")[1+clusterCols+2*clusterRows],
#       cexRow=2,
#       #margins=if(!all(is.na(data.cols))) c(7,12) else c(5,5),
#       margins=c(7,12),
#       #col=c("white",gray.colors(N_COL_BINS-1, start=0.9, end=0.3), "black"),
#       col=gray.colors(N_COL_BINS-1, start=1, end=0, gamma=1.5),
#       #breaks=c(0.0001,(1:N_COL_BINS)*(1/N_COL_BINS),0.9999),
#       breaks=(1:N_COL_BINS)*(1/N_COL_BINS))
#   if(!all(is.na(data.cols))) {
#     hm_args$ColSideColors=data.cols
#   }
#   hm = do.call(gplots::heatmap.2, hm_args)
#
#   if(!all(is.na(data.cols))){
#
#     if(is.factor(data.ch)){
#       var_levels<-levels(data.ch)
#     }else if(is.numeric(data.ch)){
#       var_levels<-levels(cut(data.ch,min(20, length(data.ch))))
#     }
#     nc<-1+length(var_levels)%/%10
#     xcoord<-0.85-0.05*nc
#     ycoord<-1.15-0.05*nc #/(length(data.ch)%%10+1)
#     legend(x=xcoord, y=ycoord, var_levels, col=palette, pch=15,
#         ncol=nc,
#         xpd=TRUE, horiz=FALSE, bty="n", title=attr(data.ch, "name"))
#   }
#   return(hm)
# }
#
#
#
#
#
# build_feats_epic_grch38 = function(pf_orig, pf_chr_colname, pf_pos_colname) {
#   # params
#   extend_region_dist = 1000
#   # meth pf
#   if (missing(pf_orig)) {
#     pf_orig = mreadRDS("~/projects/datashare/platforms/EPIC.hg38.manifest.full.fch.rds")
#     pf_chr_colname = "seqnames"
#     pf_pos_colname = "start"
#   }
#   pf_orig = pf_orig[pf_orig[,pf_pos_colname]>0,]
#   pf_orig = pf_orig[order(pf_orig[[pf_chr_colname]],pf_orig[[pf_pos_colname]]), ]
#   ## index meth probes by chr
#   chrs = unique(pf_orig[[pf_chr_colname]])
#   chrs_indexed_methpf = lapply(chrs, function(chr) {
#     print(chr)
#     idx = rownames(pf_orig)[!is.na(pf_orig[[pf_chr_colname]]) & pf_orig[[pf_chr_colname]]==chr]
#     ret = pf_orig[idx,]
#     return(ret)
#   })
#   names(chrs_indexed_methpf) = chrs
#
#   fat_feat = lapply(unique(pf_orig[,1]), function(chr) {
#     d = pf_orig[pf_orig[,1]==chr,]
#     i = intervals::Intervals(c(d[,2], d[,3]), type="Z")
#     # enlarge your fat feat
#     l = extend_region_dist
#     c = intervals::close_intervals( intervals::contract( intervals::reduce(intervals::expand(i, l)), l) )
#     dim(c)
#     df = data.frame(chr, c[,1], c[,2])
#     return(df)
#   })
#   fat_feat = do.call(rbind, fat_feat)
#   dim(fat_feat)
#   fat_feat[,4] = paste0(fat_feat[,1], ":", fat_feat[,2], "-", fat_feat[,3])
#   fat_feat[,5] = fat_feat[,3] - fat_feat[,2]
#   fat_feat[,6] = "+"
#   fat_feat = fat_feat[fat_feat[,5]>1,]
#   rownames(fat_feat) = fat_feat[,4]
#   colnames(fat_feat)  = c("chr", "start", "end", "id", "score", "strand")
#   dim(fat_feat)
#   head(fat_feat)
#
#   ## index probes by feat name
#   print("# indexing probes by feat name")
#   feat_indexed_probes = epimedtools::monitored_apply(fat_feat, 1, function(feat) {
#     # feat = fat_feat[3,]
#     # print(feat)
#     chr = feat[[1]]
#     len = as.numeric(feat[[5]])
#     meth_platform = chrs_indexed_methpf[[chr]]
#     ret = dmprocr::get_probe_names(feat, meth_platform, pf_chr_colname, pf_pos_colname, 0, len)
#     # meth_platform[ret,1:3]
#     # feat
#     return(ret)
#   })
#
#   nb_probes = sapply(feat_indexed_probes, length)
#   fat_feat$score = nb_probes[rownames(fat_feat)]
#   # saveRDS(fat_feat, "feats_epic_grch38.rds")
#   # return("feats_epic_grch38.rds")
#   # layout(matrix(1:2, 1), respect=TRUE)
#   # barplot(table(sapply(feat_indexed_probes, length)), las=2, main="#probes")
#   # plot(fat_feat$end - fat_feat$start, fat_feat$score, main="#probes")
#
#   # options(scipen=999)
#   # probes = unique(unlist(feat_indexed_probes))
#   # pf = pf_orig[probes, c(pf_chr_colname, pf_pos_colname)]
#   # pf[,3] = pf[,2] + 1
#   # pf[,4] = probes
#   # pf[,5] = 1
#   # pf[,6] = "+"
#   # write.table(pf  , file="pf.bed"  , sep="\t", quote=FALSE,row.names=FALSE, col.names=FALSE)
#   return(feat_indexed_probes)
# }
#
# if (!exists("mbuild_feats_epic_grch38")) {mbuild_feats_epic_grch38 = memoise::memoise(build_feats_epic_grch38)}
#
#
#
#
# model_call_rlm = function(meth, pheno, confounder) {
#   # print(meth)
#   # print(pheno)
#   # print(confounder)
#   #
#   # meth       <<- meth
#   # pheno      <<- pheno
#   # confounder <<- confounder
#
#   if (missing(confounder)) {
#     m = MASS::rlm(meth~pheno, maxit=400)
#   } else {
#     m = MASS::rlm(meth~pheno+confounder, data.frame(meth,pheno,confounder), maxit=400)
#   }
#   cf = lmtest::coeftest(m, vcov=sandwich::vcovHC(m, type="HC0"))
#   z = cf[2,1]/cf[2,2]
#   lpv = -(log10(2) + pnorm(abs(z), lower.tail=FALSE, log.p=TRUE)/log(10))
#   # wtest = survey::regTermTest(m, "pheno", null=NULL, df=Inf, method=c("Wald"))
#   beta = m$coefficients[[2]]
#   return(c(beta=beta, lpv=lpv))
# }
#
#
#
#
#
# # ewas = mewas_func(d=d, e, USE_PARAPPLY, pheno_key=pheno_key, confounder=confounder, model_func_name=model_func_name)
# # ewas = mewas_func(d=d[1:100,], e, USE_PARAPPLY, pheno_key=pheno_key, confounder=confounder, model_func_name=model_func_name)
#
#
#
#
# # model_call(d[1,], pheno=e[colnames(d),"tissue"], confounder=e[colnames(d),"Individu"])
# # ewas = epimedtools::monitored_apply(d, 1, model_call, pheno=e[colnames(d),"tissue"], confounder=e[colnames(d),"Individu"])
# # model_call(meth, pheno, confounder)
# # meth=d[1,]
# # model_call(meth, pheno, confounder)
#
#
#
# # ewas = mewas_func2(d=d, e=e, USE_PARAPPLY=USE_PARAPPLY, model_formula=model_formula, model_func_name=model_func_name, nb_fact_of_interest=nb_fact_of_interest)
# # ewas = mewas_func2(d=d[140341:140350,], e=e, USE_PARAPPLY=FALSE, model_formula=model_formula, model_func_name=model_func_name, nb_fact_of_interest=nb_fact_of_interest)
# # dim(ewas)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# model_call_smoke3p = function(meth, pheno, conf=NULL) {
#   if (is.null(conf)) {
#     m = try(lm(meth~pheno))
#     if (attributes(m)$class == "try-error") {
#       # warning(paste0("Region ", mdata$region_id, " contains probe ", probe, "that through error model. Have to fix it! "))
#       m = lm(rnorm(length(meth))~pheno)
#       m$coefficients[] = NA
#       m$residuals[] = NA
#     }
#     if (length(m$coefficients) != length(na.omit(unique(pheno)))) {
#       # warning(paste0("Region ", mdata$region_id, " contains probe ", probe, "that through error model. Have to fix it! "))
#       m = lm(rnorm(length(meth))~pheno)
#       m$coefficients[] = NA
#       m$residuals[] = NA
#     }
#   } else {
#     m = try(lm(meth~pheno+conf))
#     if (attributes(m)$class == "try-error") {
#       # warning(paste0("Region ", mdata$region_id, " contains probe ", probe, "that through error model. Have to fix it! "))
#       m = lm(rnorm(length(meth))~pheno+conf)
#       m$coefficients[] = NA
#       m$residuals[] = NA
#       fstat = summary(m)$fstatistic
#     } else if (length(m$coefficients) != length(na.omit(unique(pheno)))+length(na.omit(unique(conf)))-1) {
#       # warning(paste0("Region ", mdata$region_id, " contains probe ", probe, "that through error model. Have to fix it! "))
#       m = lm(rnorm(length(meth))~pheno+conf)
#       summary(m)$fstatistic
#       m$coefficients[] = NA
#       m$residuals[] = NA
#       fstat = summary(m)$fstatistic
#     } else {
#       fstat = list(anova(m)[1,4], anova(m)[1,1], anova(m)[3,1])
#     }
#   }
#   if (is.null(fstat)) {
#     # warning(paste0("Region ", mdata$region_id, " contains probe ", probe, "that through NULL fstat. Have to fix it! "))
#     m$coefficients[] = NA
#     lpval_fisher = NA
#   } else {
#     lpval_fisher = pf(fstat[[1]], fstat[[2]], fstat[[3]], lower.tail=FALSE, log.p=TRUE)/-log(10)
#   }
#
#   # ret = c(lpval_fisher=lpval_fisher, r2=summary(m)$r.squared,  nb_notna = sum(!is.na(meth)), fstat_1=fstat[[1]], fstat_2=fstat[[2]], fstat_3=fstat[[3]]) #m$coefficients[-1])
#   # return(ret)
#   lpv = lpval_fisher
#   beta = m$coefficients[[2]]
#   coef = c(unlist(m$coefficients), pheno=-sum(m$coefficients[-1]))
#   ret = c(beta=beta, lpv=lpv, coef=coef)
#   ret
#   return(ret)
# }
#
#
#
#
#
# model_call_var = function(meth, pheno, conf=NULL) {
#   m = try(var.test(meth[pheno%in%unique(na.omit(pheno))[1]], meth[pheno%in%unique(na.omit(pheno))[2]]) )
#   # set.seed(1)
#   # m =     var.test(rnorm(10,0,1)                           , rnorm(10,0,2)                           )
#   # m$p.value
#   # pv = 1-pf(m$statistic[[1]], m$parameter[1], m$parameter[2], lower.tail=FALSE)
#   # pv
#   # -log10(pv)
#   # -log10(2)+pf(m$statistic[[1]], m$parameter[1], m$parameter[2], lower.tail=FALSE, log.p=TRUE)/-log(10)
#
#   if (attributes(m)$class == "try-error") {
#     lpv = 0
#     beta = NA
#     coef = NA
#   } else {
#     # lpv = max(0, -log10(2)+pf(m$statistic[[1]], m$parameter[1], m$parameter[2], lower.tail=FALSE, log.p=TRUE)/-log(10))
#     lpv = -log10(m$p.value)
#     beta = m$estimate[[1]]
#     coef = NA
#   }
#   ret = c(beta=beta, lpv=lpv, coef=coef)
#   ret
#
#   return(ret)
# }
#
# # idx_probes = rownames(d)[395315:395317]
# # ewas_func(d=d[idx_probes,], e, USE_PARAPPLY, pheno_key=pheno_key, confounder=confounder, model_func_name="model_call_var")
# #
# # meth = d[idx_probes,][1,]
#
#
#
# # if (!exists("mget_meth_by_gene")) mget_meth_by_gene = memoise::memoise(get_meth_by_gene)
