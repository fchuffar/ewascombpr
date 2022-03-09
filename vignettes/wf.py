import os 

localrules: target

rule target:
    threads: 1
    message: "-- Rule target completed. --"
    input: 
      os.getcwd()+"/02_ewas_combp_study_GSE40279.rds_modelcalllm_meth~age__y_+gender_1e-30.html" ,
    shell:"""
pwd
          """

include: "basic_rules.py"

