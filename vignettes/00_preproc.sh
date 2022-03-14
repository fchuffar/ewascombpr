cd ~/projects/expedition_5300/results/GSE162606/vignettes
source config
rsync -auvP ~/projects/${project}/results/${study}/ cargo:~/projects/${project}/results/${study}/
  
# launch default pipeline
snakemake -s wf.py -pn

# launch default pipeline
cp wf.py wf_${project}.py
cp basic_rules.py ${project}_rules.py

