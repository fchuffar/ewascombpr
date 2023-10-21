cd ~/projects/ewascombpr/vignettes
rsync -auvP ~/projects/ewascombpr/ cargo:~/projects/ewascombpr/
source config
echo ${study}
echo ${project}
rsync -auvP ~/projects/${project}/results/${study}/ cargo:~/projects/${project}/results/${study}/ --dry-run

# launch default pipeline
snakemake --cores 1 -s wf.py -pn

# launch custom pipeline
cp 00_preproc.sh 00_preproc_local.sh
cp wf.py 00_wf_local.py
cp rules.py 00_rules_local.py
cp config config_local

snakemake --cores 1 -s 00_wf_local.py -pn



source ~/conda_config.sh 
cd ~/projects/breast/results/04.2_ewascombpr/vignettes
echo $PYTHONPATH
# PYTHONPATH="${PYTHONPATH}:/summer/epistorage/opt/combined-pvalues/
Sys.setenv(PYTHONPATH = "/summer/epistorage/opt/combined-pvalues/")
/summer/epistorage/opt/combined-pvalues/cpv/comb-p pipeline -c 5 --seed 1e-30 --dist 1000 -p dmr_study_TCGA-BRCA.rds_modelcalllm_meth~gec_1e-30 --region-filter-p 0.05 --region-filter-n 2 ewas4combp_study_TCGA-BRCA.rds_modelcalllm_meth~gec.bed


export PYTHONPATH="/summer/epistorage/opt/combined-pvalues/"
pip install toolshed interlap




s = readRDS("study_r0_ewas2000_nn1000_GSE42861.rds")
idx = readRDS("df_r0_ewas2000_nn1000_GSE42861.rds")
colnames(idx)
s$data = s$data[rownames(s$data)%in%colnames(idx),]
s$platform = s$platform[rownames(s$data),]
saveRDS(s, "study_r0_ewas2000_nn1000_GSE42861.rds")

