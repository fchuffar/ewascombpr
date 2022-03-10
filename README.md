# ewascombpr
A line that pipes epimedtools DNA methylation study to an ewas then to combp.

# Quick launch

## Under your terminal

```
cd ~/projects
git clone https://github.com/fchuffar/ewascombpr.git
cd ewascombpr/vignettes
R 
```

## Under R

```
# build Hannum 2013 study
rmarkdown::render("vignettes/01_build_study_geo.Rmd")
# run pipeline
rmarkdown::render("vignettes/02_ewas_combp.Rmd")
```



# Usage
  
```
cd ~/projects/ewascombpr/

# push scripts on the server
rsync -auvP ~/projects/ewascombpr/ cargo:~/projects/ewascombpr/

# build GSE40279 methylation study
echo 'rmarkdown::render("vignettes/01_build_study_geo.Rmd")' | Rscript -

# analyse GSE40279 methylation study (meth~age+gender)
echo 'rmarkdown::render("vignettes/02_ewas_combp.Rmd")' | Rscript -

# build and analyse GSE40279 methylation study with sankemake
cd vignettes/
snakemake -s wf.py -pn
```
