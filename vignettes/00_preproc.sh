cd ~/projects/smoke3p/results/ewascombpr/vignettes
source config
rsync -auvP ~/projects/${project}/results/ewascombpr/ cargo:~/projects/${project}/results/ewascombpr/
  
# launch default pipeline
snakemake -s wf.py -pn

# launch default pipeline
cp wf.py wf_${project}.py
cp basic_rules.py ${project}_rules.py

