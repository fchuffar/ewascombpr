# rule ewas_neighb:
#     input:
#       bed_ewas = "{prefix}/ewas4combp_{study_filename}.rds_{modelcall}_meth~{model_formula}.bed",
#       rmd_script="{prefix}/01_ewas_neighb.Rmd",
#     output:
#       study_rds = "{prefix}/{study_filename}_{model_func_name}_{model_formula}_ewas{newas}_nn{neighb}.rds",
#       html =      "{prefix}/01_ewas_neighb_{study_filename}_{model_func_name}_{model_formula}_ewas{newas}_nn{neighb}.html"      ,
#     threads: 31
#     shell:"""
# export PATH="/summer/epistorage/opt/bin:$PATH"
# export PATH="/summer/epistorage/miniconda3/bin:$PATH"
# cd {wildcards.prefix}
# echo "study_filename='{wildcards.study_filename}.rds' ; model_func_name='{wildcards.model_func_name}' ; model_formula='{wildcards.model_formula}' ; newas={wildcards.newas}; neighb={wildcards.neighb}; rmarkdown::render('{input.rmd_script}', output_file='{output.html}')" | Rscript -
# """



rule ewas:
    input: 
      # rmd_script="{rmd_script_prefix}.Rmd",
      study="{prefix}/{study_filename}.rds",
      rmd_script="{prefix}/01_ewas.Rmd",      
    output: 
      html = "{prefix}/01_ewas_{study_filename}.rds_{modelcall}_meth~{model_formula}.html"      ,           
      bed_ewas = "{prefix}/ewas4combp_{study_filename}.rds_{modelcall}_meth~{model_formula}.bed",           
      rds_ewas = "{prefix}/ewas_{study_filename}.rds_{modelcall}_meth~{model_formula}.rds"      ,           
    threads: 32
    shell:"""
export PATH="/summer/epistorage/opt/bin:$PATH"
export PATH="/summer/epistorage/miniconda3/bin:$PATH"
cd {wildcards.prefix}
echo "study_filename='{wildcards.study_filename}.rds'; model_func_name = '{wildcards.modelcall}' ; model_formula='meth~{wildcards.model_formula}' ; rmarkdown::render('{input.rmd_script}', output_file=paste0('01_ewas_', study_filename, '_', model_func_name, '_', model_formula, '.html'))" | Rscript -
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
      html = "{prefix}/01_build_study_GSE{gse}.html"      ,           
      study = "{prefix}/study_GSE{gse}.rds"      ,           
    threads: 1
    shell:"""
export PATH="/summer/epistorage/opt/bin:$PATH"
export PATH="/summer/epistorage/miniconda3/bin:$PATH"
cd {wildcards.prefix}
echo "gse = 'GSE{wildcards.gse}'; rmarkdown::render('{input.rmd_script}', output_file='{output.html}')" | Rscript -
"""


rule build_rfe_study:
    input:
      rmd_script="{prefix}/01_RefFreeEWAS.Rmd",
      study_orig="{prefix}/study_{studysuffix}.rds"      ,
    output:
      study_rfe="{prefix}/study_rfe{studysuffix}.rds"      ,
    threads: 1
    shell:"""
export PATH="/summer/epistorage/opt/bin:$PATH"
export PATH="/summer/epistorage/miniconda3/bin:$PATH"
cd {wildcards.prefix}
echo "study_filename='{input.study_orig}' ; confounder=NULL ; sd_thresh=0.05 ; output_file=paste0('05_RefFreeEWAS_', sd_thresh, '_', study_filename, '.html') ; rmarkdown::render('{input.rmd_script}', output_file=output_file)" | Rscript -
"""
