cd ~/projects/ewascombpr/results/ewascombpr/vignettes
source config
echo ${study}
echo ${project}
rsync -auvP ~/projects/${project}/results/${study}/ cargo:~/projects/${project}/results/${study}/
  
# launch default pipeline
snakemake -s wf.py -pn

# launch default pipeline
cp wf.py 00_local_wf.py
cp rules.py 00_local_rules_.py

snakemake -s 00_local_wf.py -pn