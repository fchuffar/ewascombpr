cd ~/projects/expedition_5300/results/cmslimapuno/vignettes
source config
echo ${study}
echo ${project}
rsync -auvP ~/projects/${project}/results/${study}/ cargo:~/projects/${project}/results/${study}/
  
# launch default pipeline
snakemake -s wf.py -pn

# launch default pipeline
cp wf.py wf_${study}.py
cp basic_rules.py ${study}_rules.py

snakemake -s wf_${study}.py -pn