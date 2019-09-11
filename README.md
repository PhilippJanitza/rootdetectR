# rootdetectR
How To Install rootdetectR

Make sure you have the latest version of R and RStudio.
__Caution:__ Windows user will need to install [RTools](cran.r-project.org/bin/windows/Rtools/) to be able to biuld the package.

To build the vignette on you computer you will need to install the packages knitr and rmarkdown.
```r
install.packages(c('knitr', 'rmarkdown'))
```

To install rootdetectR from github the package devtools is required. 

```r
# devtools needs to be available to download packages from github
if('devtools' %in% rownames(installed.packages()) == FALSE) {install.packages('devtools')}
library('devtools')

#install from github
install_github("PhilippJanitza/rootdetectR", build_vignettes = TRUE)
```


