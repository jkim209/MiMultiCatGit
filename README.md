# MiMultiCatGit

Title: MiMultiCat: A unified web cloud analytic platform for user-friendly and interpretable microbiome data mining using tree-based methods


Version: 1.0.1

Maintainer: Jihun Kim <toujours209@gmail.com>

Description: MiMultiCat is a unified web cloud analytic platform for user-friendly and interpretable microbiome data mining using tree-based methods. Machine learning is a promising approach to help such an effort especially due to the high complexity of microbiome data. However, many of the current machine learning algorithms are in a “black box”. They are hard to understand and interpret. Clinicians, public health practitioners or biologists are not also usually skilled at computer programming, and they do not always have a high-end computing device. MiMultiCat employs tree-based learning methods, including 1) decision tree, 2) random forest and 3) gradient boosting, that are both well understood and suited to human microbiome studies, for both classification and regression problems through covariate-adjusted or unadjusted analysis. 

NeedsCompilation: No

Depends: R(≥ 4.1.0)
**(1.4.0.1)**
Imports: Bioconductor ('BiocManager', 'biomformat', 'phyloseq'); CRAN ('shiny', 'rmarkdown', 'seqinr', 'shinydashboard', 'tidyverse', 'plotly', 'shinyWidgets', 'shinyjs', 'googleVis', 'xtable',
            'DT', 'htmltools', 'phangorn', 'bios2mds', 'zip', 'ape', 'zCompositions', 'compositions', 'stringr', 'caret', 'ggplot2', 'randomForest', 
            'data.table', 'xgboost', 'SHAPforxgboost', 'fontawesome', 'grid', 'ggplotify', 'remotes', 'doParallel', 'reshape2', 'fossil', 
            'proxy', 'ecodist', 'GUniFrac', 'picante', 'FSA', 'tibble', 'forestplot', 'VGAM', 'rgl', 'MiRKAT'); GitHub: ('dashboardthemes, 'edarf', 'MiRKATMC', 'vegan')

License: GPL 1, GPL 2 

## URLs

* Web application (online implementation): http://mimulticat.micloud.kr
* GitHub repository (local implementation): https://github.com/jkim209/MiMultiCatGit

## References

* Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)


## Prerequites

#### Notice: For the local implementation, you do not need to install all the pre-requite R packages individually. You only need to install the 'shiny' package, and then run a simple command in 'Launch App' below. Then, all the pre-requisite R packages will be installed and imported automatically. 
#### For MAC users, please make sure that you already dowload and install xQuartz (https://www.xquartz.org/) to your device before the local implementation. 


shiny
```
install.packages('shiny')
```

## Launch App

```
library(shiny)

runGitHub('MiMultiCatGit', 'jkim209', ref = 'main')
```

## Troubleshooting Tips

If you have any problems for using MiMultiCat, please report in issues (https://github.com/jkim209/MiMultiCatGit/issues) or email Jihun Kim (toujours209@gmail.com).
