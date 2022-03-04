# ewascombpr
A line that pipes epimedtools DNA methylation study to an ewas then to combp


# Usage
  
```
cd ~/projects/ewascombpr/
rsync -auvP ~/projects/ewascombpr/ cargo:~/projects/ewascombpr/

echo 'rmarkdown::render("vignettes/05_sensitivity_analysis.Rmd")' | Rscript -


```