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

## Citation

Cao, Y., Lin, Y., Ormerod, J. T., Yang, P., Yang, J. Y. H., & Lo, K. K. (2019). scDC: single cell differential composition analysis. In BMC Bioinformatics (Vol. 20, Issue S19). Springer Science and Business Media LLC. https://doi.org/10.1186/s12859-019-3211-9 
