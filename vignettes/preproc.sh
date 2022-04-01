cd ~/projects/ewascombpr/results/ewascombpr/vignettes
source config
echo ${study}
echo ${project}
rsync -auvP ~/projects/${project}/results/${study}/ cargo:~/projects/${project}/results/${study}/

# launch default pipeline
snakemake -s wf.py -pn

# launch default pipeline
cp preproc.sh 00_preproc.sh
cp wf.py 00_wf_local.py
cp rules.py 00_rules_local.py

snakemake -s 00_wf_local.py -pn