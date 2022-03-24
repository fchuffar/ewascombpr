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


