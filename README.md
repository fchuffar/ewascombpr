# ewascombpr
A line that pipes epimedtools DNA methylation study to an ewas then to combp.


# Usage
  
```
cd ~/projects/ewascombpr/

# push scripts on the server
rsync -auvP ~/projects/ewascombpr/ cargo:~/projects/ewascombpr/

# build GSE40279 methylation study
cd ~/projects/ewascombpr/
echo 'rmarkdown::render("vignettes/01_build_study_geo.Rmd")' | Rscript -

# analyse GSE40279 methylation study (meth~age+gender)
cd ~/projects/ewascombpr/
echo 'rmarkdown::render("vignettes/02_ewas_combp.Rmd")' | Rscript -

# build and analyse GSE40279 methylation study with sankemake
cd ~/projects/ewascombpr/vignettes/
snakemake -s wf.py -pn
```