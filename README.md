# rootdetectR
Install rootdetectR

```r
# devtools needs to be available to download packages from github
if('devtools' %in% rownames(installed.packages()) == FALSE) {install.packages('devtools')}
library('devtools')

#install from github
install_github("PhilippJanitza/rootdetectR", build_vignettes = TRUE)
```
