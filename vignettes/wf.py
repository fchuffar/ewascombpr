import os 

localrules: target

rule target:
    threads: 1
    message: "-- Rule target completed. --"
    input: 
      "/home/chuffarf/projects/ewascombpr/vignettes/02_ewas_combp_study_GSE40279.rds_modelcalllm_meth~age__y_+gender_1e-30.html" ,
    shell:"""
pwd
          """

rule ewas_combp:
    input: 
      # rmd_script="{rmd_script_prefix}.Rmd",
      study="{prefix}/{study_filename}.rds",
      rmd_script="{prefix}/02_ewas_combp.Rmd",      
    output: 
      html = "{prefix}/02_ewas_combp_{study_filename}.rds_{modelcall}_meth~{model_formula}_{pval_thresh}.html"      ,           
    threads: 31
    shell:"""
export PATH="/summer/epistorage/opt/bin:$PATH"
export PATH="/summer/epistorage/miniconda3/bin:$PATH"
cd {wildcards.prefix}
echo "study_filename='{wildcards.study_filename}.rds'; model_func_name = '{wildcards.modelcall}' ; model_formula='meth~{wildcards.model_formula}'; pval_thresh='{wildcards.pval_thresh}'; rmarkdown::render('{input.rmd_script}', output_file=paste0('02_ewas_combp_', study_filename, '_', model_func_name, '_', model_formula, '_', pval_thresh, '.html'))" | Rscript -
"""


rule build_study_geo:
    input: 
      rmd_script="{prefix}/01_build_study_geo.Rmd",      
    output: 
      html = "{prefix}/01_build_study_{gse}.html"      ,           
      study = "{prefix}/study_{gse}.rds"      ,           
    threads: 1
    shell:"""
export PATH="/summer/epistorage/opt/bin:$PATH"
export PATH="/summer/epistorage/miniconda3/bin:$PATH"
cd {wildcards.prefix}
echo "gse = '{wildcards.gse}'; rmarkdown::render('{input.rmd_script}', output_file=paste0('01_build_study_', gse ,'.html'))" | Rscript -
"""


