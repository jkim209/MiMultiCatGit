# MiMultiCatGit

Title: MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses


Version: 1.0.0

Maintainer: Jihun Kim <toujours209@gmail.com>

Description: MiMultiCat is a unified cloud platform for the analysis of microbiome data with multi-categorical responses. The two key features of MiMultiCat are as follows. First, MiMultiCat streamlines a long sequence of microbiome data preprocessing and analytic procedures on user-friendly web interfaces; as such, it is easy to use for many people in various disciplines (e.g., biology, medicine, public health). Second, MiMultiCat performs both association testing and prediction modeling extensively. For association testing, MiMultiCat handles both ecological (e.g., alpha- and beta-diversity) and taxonomical (e.g., phylum, class, order, family, genus, species) contexts through covariate-adjusted or unadjusted analysis. For prediction modeling, MiMultiCat employs random forest and gradient boosting algorithms that are well-suited to microbiome data with nice visual interpretations. 

NeedsCompilation: No

Depends: R(≥ 4.1.0)

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
