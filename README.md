# scDC: Single cell differential composition analysis

`scDC` is a package to perform differential composition analysis on scRNA-seq data.


## Installation

```
## Some CRAN packages required by scDC
install.packages(c("parallel", "DescTools", "lme4", "reshape2", "ggridges", 
"lme4", "mice", "broom.mixed"))

## Some BioConductor packages required by scDC
BiocManager::install(c("scran"))

## Some Github packages required by scDC
devtools::install_github("taiyunkim/scClustBench")


## Installing scDC 
devtools::install_github("SydneyBioX/scDC")
```


## Vignette

You can find the vignette at our website: https://sydneybiox.github.io/scDC/index.html.


## Authors

Yingxin Lin
Yue Cao
