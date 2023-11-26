ls.pkg <- c('shiny', 'BiocManager', 'rmarkdown', 'seqinr', 'shinydashboard', 'tidyverse', 'plotly', 'shinyWidgets', 'shinyjs', 'googleVis', 'xtable',
            'DT', 'htmltools', 'phangorn', 'bios2mds', 'zip', 'ape', 'zCompositions', 'compositions', 'stringr', 'caret', 'ggplot2', 'randomForest', 
            'data.table', 'xgboost', 'SHAPforxgboost', 'fontawesome', 'grid', 'ggplotify', 'remotes', 'doParallel', 'reshape2', 'fossil', 
            'proxy', 'ecodist', 'GUniFrac', 'picante', 'FSA', 'tibble', 'forestplot', 'VGAM', 'rgl', 'MiRKAT')

new.pkg <- ls.pkg[!(ls.pkg %in% installed.packages()[,"Package"])]
if(length(new.pkg)) install.packages(new.pkg, repos = 'https://cloud.r-project.org/')
if(packageVersion("zCompositions") != "1.4.0.1") remotes::install_version("zCompositions", version = "1.4.0.1", repos = "http://cran.us.r-project.org")

if(!require('phyloseq')) BiocManager::install('phyloseq')
if(!require('biomformat')) remotes::install_github('joey711/biomformat')
if(!require('dashboardthemes')) remotes::install_github('nik01010/dashboardthemes', force = TRUE)
if(!require('edarf')) remotes::install_github('zmjones/edarf', subdir = 'pkg')
if(!require('MiRKATMC')) remotes::install_github("Zhiwen-Owen-Jiang/MiRKATMC")
if(!require('vegan')) remotes::install_github("vegandevs/vegan")

library(shiny)
library(BiocManager)
library(seqinr)
library(shinydashboard)
library(tidyverse)
library(plotly)
library(shinyWidgets)
library(shinyjs)
library(googleVis)
library(xtable)
library(DT)
library(htmltools)
library(phangorn)
library(bios2mds)
library(zip)
library(ape)
library(zCompositions)
library(compositions)
library(stringr)
library(caret)
library(ggplot2)
library(randomForest)
library(data.table)
library(xgboost)
library(SHAPforxgboost)
library(fontawesome)
library(grid)
library(ggplotify)
library(remotes)
library(FSA)
library(VGAM)

source("Source/MiDataProc.Data.Upload.R")
source("Source/MiDataProc.Data.Input.R")
source("Source/MiDataProc.Alpha.Diversity.R")
source("Source/MiDataProc.Beta.Diversity.R")
source("Source/MiDataProc.Taxa.R")
source("Source/MiDataProc.ML.Models.R")
source("Source/MiDataProc.ML.RF.R")
source("Source/MiDataProc.ML.XGB.R")
source("Source/setSliderColor.R")
options(scipen=999)

# COMMENTS ------
{
  ## HOME COMMENTS -----
  
  TITLE = p("MiMultiCat: A Unified Cloud Platform for the Analysis of Microbiome Data with Multi-Categorical Responses", style = "font-size:17pt")
  HOME_COMMENT = p(strong("MiMultiCat", style = "font-size:12pt"), "is a unified cloud platform for the analysis of microbiome data with multi-categorical responses. 
                   The two key features of MiMultiCat are as follows. First, MiMultiCat streamlines a long sequence of microbiome data preprocessing and analytic procedures on user-friendly web interfaces; 
                   as such, it is easy to use for many people in various disciplines (e.g., biology, medicine, public health). 
                   Second, MiMultiCat performs both association testing and prediction modeling extensively. 
                   For association testing, MiMultiCat handles both ecological (e.g., alpha- and beta-diversity) and taxonomical (e.g., phylum, class, order, family, genus, species) contexts through covariate-adjusted or unadjusted analysis. 
                   For prediction modeling, MiMultiCat employs random forest and gradient boosting algorithms that are well-suited to microbiome data with nice visual interpretations.", style = "font-size:12pt")
  HOME_COMMENT2 = p(strong("URLs:"), "Web server (online implementation):", tags$a(href = "http://mimulticat.micloud.kr", "http://mimulticat.micloud.kr"), 
                    "; GitHub repository (local implementation):", tags$a(href = "https://github.com/jkim209/MiMultiCatGit", "https://github.com/jkim209/MiMultiCatGit"), style = "font-size:12pt")
  HOME_COMMENT3 = p(strong("Maintainers:"), "Jihun Kim (", tags$a(href = "toujours209@gmail.com", "toujours209@gmail.com"), ")", style = "font-size:12pt")
  HOME_COMMENT4 = p(strong("Reference:"), "Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", style = "font-size:12pt")
  
  # INPUT COMMENTS -----
  
  INPUT_PHYLOSEQ_COMMENT1 = p("Description:", br(), br(), 
                              "This should be an '.Rdata' or '.rds' file, and the data should be in the 'phyloseq' format (see ", 
                              htmltools::a(tags$u("https://bioconductor.org/packages/release/bioc/html/phyloseq.html"), style = "color:red3"),
                              "). The phyloseq object should contain all the four necessary data, feature (OTU or ASV) table, taxonomic table, meta/sample information, and phylogenetic tree",br(), br(),
                              
                              "Details:",br(), br(), 
                              "1) The feature table should contain counts, where rows are features (OTUs or ASVs) 
                              and columns are subjects (row names are feature IDs and column names are subject IDs).", br(), br(),
                              
                              "2) The taxonomic table should contain taxonomic names, where rows are features and columns are seven taxonomic ranks 
                              (row names are feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' or 
                              'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species').", br(), br(),
                              
                              "3) The metadata/sample information should contain variables for the subjects about host phenotypes, medical interventions, disease status or environmental/behavioral factors, 
                              where rows are subjects and columns are variables (row names are subject IDs, and column names are variable names).", br(), br(),
                              
                              "4) The phylogenetic tree should be a rooted tree. Otherwise, MiMultiCat automatically roots the tree through midpoint rooting (phangorn::midpoint). 
                              The tip labels of the phylogenetic tree are feature IDs.", br(), br(),
                              
                              "* The features should be matched and identical across feature table, taxonomic table and phylogenetic tree. 
                              The subjects should be matched and identical between feature table and metadata/sample information. 
                              MiMultiCat will analyze only the matched features and subjects.", style = "font-size:11pt")
  
  INPUT_PHYLOSEQ_COMMENT2 = p("You can download example microbiome data 'biom.Rdata' in the  'phyloseq' format. For more details about 'phyloseq', see ", 
                              htmltools::a(tags$u("https://bioconductor.org/packages/release/bioc/html/phyloseq.html"), style = "color:red3"), br(), br(), 
                              
                              "> setwd('/yourdatadirectory/')", br(), br(), 
                              "> load(file = 'biom.Rdata')", br(), br(), 
                              "> library(phyloseq)", br(), br(), 
                              " > otu.tab <- otu_table(biom)",  br(), 
                              " > tax.tab <- tax_table(biom)", br(), 
                              " > sam.dat <- sample_data(biom)", br(), 
                              " > tree <- phy_tree(biom)",br(), br(),
                              
                              "You can check if the features are matched and identical across feature table, taxonomic table and phylogenetic tree, 
                              and the subjects are matched and identical between feature table and metadata/sample information using following code.", br(), br(), 
                              
                              " > identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                              " > identical(colnames(otu.tab), tree$tip.label)", br(),
                              " > identical(colnames(otu.tab), rownames(sam.dat))", style = "font-size:11pt", br(), br(),
                              
                              strong("Reference:"), "Goodrich JK, Waters JL, Poole AC, Sutter JL, Koren O, Blekhman R, et al. Human genetics shape the gut microbiome. Cell. 2014:159(4):789-799.")
  
  INPUT_INDIVIDUAL_DATA_COMMENT = p("Description:", br(), br(), 
                                    "1) The feature table (.txt or .csv) should contain counts, where rows are features (OTUs or ASVs) and columns are subjects (row names are feature IDs and column names are subject IDs). 
                                    Alternatively, you can upload .biom file processed by QIIME.", br(), br(),
                                    
                                    "2) The taxonomic table (.txt) should contain taxonomic names, where rows are features and columns are seven taxonomic ranks 
                                    (row names are feature IDs and column names are 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species' or 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'). 
                                    Alternatively, you can upload .tsv file processed by QIIME.", br(), br(),
                                    
                                    "3) The metadata/sample information (.txt or .csv) should contain variables for the subjects about host phenotypes, medical interventions, disease status or environmental/behavioral factors, 
                                    where rows are subjects and columns are variables (row names are subject IDs, and column names are variable names).", br(), br(),
                                    
                                    "4) The phylogenetic tree (.tre or .nwk) should be a rooted tree. Otherwise, MiMultiCat automatically roots the tree through midpoint rooting (phangorn::midpoint). 
                                    The tip labels of the phylogenetic tree are feature IDs." ,br(), br(),
                                    
                                    "* The features should be matched and identical across feature table, taxonomic table and phylogenetic tree. 
                                    The subjects should be matched and identical between feature table and metadata/sample information. 
                                    MiMultiCat will analyze only the matched features and subjects.", style = "font-size:11pt")
  
  INPUT_INDIVIDUAL_DATA_COMMENT2 = p("You can download example microbiome data 'biom.zip'. This zip file contains four necessary data, 
                                     feature table (otu.tab.txt), taxonomic table (tax.tab.txt), and metadata/sample information (sam.dat.txt), and phylogenetic tree (tree.tre).", br(), br(),
                                     
                                     "> setwd('/yourdatadirectory/')", br(), br(), 
                                     "> otu.tab <- read.table(file = 'otu.tab.txt', check.names = FALSE)", br(), 
                                     "> tax.tab <- read.table(file = 'tax.tab.txt', check.names = FALSE)", br(), 
                                     "> sam.dat <- read.table(file = 'sam.dat.txt', check.names = FALSE)", br(),
                                     "> tree <- read.table(file = 'tree.tre')", br(), br(),
                                     
                                     "You can check if the features are matched and identical across feature table, taxonomic table and phylogenetic tree, 
                                     and the subjects are matched and identical between feature table and metadata/sample information using following code.", br(), br(), 
                                     
                                     " > identical(rownames(otu.tab), rownames(tax.tab))", br(), 
                                     " > identical(rownames(otu.tab), tree$tip.label)", br(), 
                                     " > identical(colnames(otu.tab), rownames(sam.dat))", style = "font-size:11pt", br(), br(),
                                     
                                     strong("Reference:"), "Goodrich JK, Waters JL, Poole AC, Sutter JL, Koren O, Blekhman R, et al. Human genetics shape the gut microbiome. Cell. 2014:159(4):789-799.")
  
  # QC COMMENTS -----
  
  QC_KINGDOM_COMMENT = p("A microbial kingdom to be analyzed. Default is 'Bacteria' for 16S data. 
                         Alternatively, you can type 'Fungi' for ITS data or any other kingdom of interest for shotgun metagenomic data.", style = "font-size:11pt")
  QC_LIBRARY_SIZE_COMMENT1 = p("Remove units that have low library sizes (total read counts). Default is 3,000.", style = "font-size:11pt")
  QC_LIBRARY_SIZE_COMMENT2 = p("Library size: The total read count per unit.", style = "font-size:11pt")
  QC_MEAN_PROP_COMMENT1 = p("Remove features (OTUs or ASVs) that have low mean relative abundances (Unit: %). Default is 0.02% (0.0002).",style = "font-size:11pt")
  QC_MEAN_PROP_COMMENT2 = p("Mean proportion: The average of relative abundances (i.e., proportions) per feature.", style = "font-size:11pt")
  QC_TAXA_NAME_COMMENT1 = p("Remove taxonomic names in the taxonomic table that are completely matched with the specified character strings. 
                            Multiple character strings should be separated by a comma. 
                            Default is \"\", \"metagenome\", \"gut metagenome\", \"mouse gut metagenome\".",
                            style = "font-size:11pt")
  QC_TAXA_NAME_COMMENT2 = p("Remove taxonomic names in the taxonomic table that are partially matched with the specified character strings 
                            (i.e., taxonomic names that contain the specified character strings). Multiple character strings should be separated by a comma. 
                            Default is \"uncultured\", \"incertae\", \"Incertae\", \"unidentified\", \"unclassified\", \"unknown\".",
                            style = "font-size:11pt")
  
  # MiMultiCat Reference -----
  
  MiMultiCat_REFERENCE <- p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)")
  
  # ALPHA COMMENTS -----
  
  ALPHA_COMMENT = p("Calculate alpha-diversity indices: Richness (Observed), Shannon (Shannon, 1948), Simpson (Simpson, 1949), Inverse Simpson (Simpson, 1949), 
                    Fisher (Fisher et al., 1943), Chao1 (Chao, 1984), ACE (Chao and Lee, 1992), ICE (Lee and Chao, 1994), PD (Faith, 1992).")
  ALPHA_REFERENCES = p("1. Chao A, Lee S. Estimating the number of classes via sample coverage. J Am Stat Assoc. 1992:87:210-217.", br(),
                       "2. Chao A. Non-parametric estimation of the number of classes in a population. Scand J Stat. 1984:11:265-270.", br(),
                       "3. Faith DP. Conservation evaluation and phylogenetic diversity. Biol Conserv. 1992:61:1-10.", br(),
                       "4. Fisher RA, Corbet AS, Williams CB. The relation between the number of species and the number of individuals 
                       in a random sample of an animal population. J Anim Ecol. 1943:12:42-58.", br(),
                       "5. Lee S, Chao A. Estimating population size via sample coverage for closed capture-recapture models. Biometrics. 1994:50:1:88-97.", br(),
                       "6. Shannon CE. A mathematical theory of communication. Bell Syst Tech J. 1948:27:379-423 & 623-656.", br(),
                       "7. Simpson EH. Measurement of diversity. Nature 1949:163:688.", br())
  ALPHA_ANOVAF_REFERENCE = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                             "2. Tukey JW. Commparing Individual Means in the Analysis of Variance. Biometrics. 1949;5(2):99-114", br())
  ALPHA_KRUSKAL_REFERENCE = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                              "2. Kruskal WH, Wallis WA. Use of Ranks in One-Criterion Variance Analysis. Journal of the American Statistical Association. 1952;47(260):583-621", br(), 
                              "3. Dunn OH. Multiple Comparisons Using Rank Sums. Technometrics. 1964;6(3):241-252")
  ALPHA_PROPODDS_REFERENCE = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                               "2. McCullaph P. Regression models for ordinal data. J R Stat Soc Series B. 1980;42(2):109-142.")
  ALPHA_MULTINOM_REFERENCE <- p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)")
  
  # BETA COMMENTS -----
  
  BETA_COMMENT = p("Calculate beta-diversity indices: Jaccard dissimilarity (Jaccard, 1912), Bray-Curtis dissimilarity (Bray and Curtis, 1957), Unweighted UniFrac distance 
                   (Lozupone and Knight, 2005), Generalized UniFrac distance (Chen et al., 2012), Weighted UniFrac distance (Lozupone et al., 2007).")
  BETA_REFERENCES = p("1. Bray JR, Curtis JT. An ordination of the upland forest communities of Southern Wisconsin. Ecol Monogr. 1957;27(32549).", br(),
                      "2. Chen J, Bittinger K, Charlson ES, Hoffmann C, Lewis J, Wu GD., et al. Associating microbiome composition with environmental 
                      covariates using generalized UniFrac distances. Bioinformatics. 2012;28(16):2106-13.", br(),
                      "3. Jaccard P. The distribution of the flora in the alpine zone. New Phytol. 1912;11(2):37-50.", br(),
                      "4. Lozupone CA, Hamady M, Kelley ST, Knight R. Quantitative and qualitative β-diversity measures lead to 
                      different insights into factors that structure microbial communities. Appl Environ Microbiol. 2007;73(5):1576-85.", br(),
                      "5. Lozupone CA, Knight R. UniFrac: A new phylogenetic method for comparing microbial communities. Appl Environ Microbiol. 2005;71(12):8228-35.")
  BETA_DA_REFERENCE = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                        "2.Jiang Z, He M, Chen J, Zhao N, Zhan X. MiRKAT-MC: A Distance-Based Microbiome Kernel Association Test With Multi-Categorical Outcomes. Front Genet. 2022;13:841764.")
  BETA_PERMANOVA_REFERENCE = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                               "2. Anderson MA. new method for non-parametric multivariate analysis of variance. Austral Ecology. 2001;26(1):32-46.", br(),
                               "3. McArdle BH, Anderson MJ. Fitting multivariate models to community data: A comment on distance-based redundancy analysis. Ecology. 2001;82(1):290-297.")
  
  # DT COMMENTS -----
  
  DATA_TRANSFORM_COMMENT = p("Transform the data into four different formats (1) CLR (centered log ratio) (Aitchison, 1982), (2) Count (Rarefied) (Sanders, 1968), (3) Proportion, (4) Arcsine-root 
                             for each taxonomic rank (phylum, class, order, familiy, genus, species).")
  DATA_TRANSFORM_REFERENCE = p("1. Aitchison J. The statistical analysis of compositional data. J R Stat Soc B. 1982;44(2):139-77", br(),
                               "2. Sanders HL. Marine benthic diversity: A comparative study. Am Nat. 1968;102:243-282.")
  
  # RF COMMENTS -----
  
  RF_REFERENCE = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                   "2. Breiman L. Random forests. Mach Learn. 2001;45:5-32", br())
  RF_REFERENCE_CLR = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                       "2. Breiman L. Random forests. Mach Learn. 2001;45:5-32", br(),
                       "3. Aitchison J. The statistical analysis of compositional data. J R Stat Soc B. 1982;44(2):139-77")
  RF_REFERENCE_RC = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                      "2. Breiman L. Random forests. Mach Learn. 2001;45:5-32", br(),
                      "3. Sanders HL. Marine benthic diversity: A comparative study. Am Nat. 1968;102:243-282.")
  
  # XGB COMMENTS -----
  
  XGB_REFERENCE = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                    "2. Friedman JH. Greedy function approximation: A gradient boosting machine. Ann Stat. 2001;29(5):1189-1232",br(),
                    "3. Chen T, Guestrin C. XGBoost: A scalable tree boosting system. in Proc the 22nd ACM SIGKDD Int Conf KDD. ACM. 2016;785-794", br(),
                    "4. Lundberg SM, Lee SI. A unified approach to interpreting model predictions. in Proc Adv Neural Inf Process Syst. 2017;4765-4774.")
  XGB_REFERENCE_CLR = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                        "2. Friedman JH. Greedy function approximation: A gradient boosting machine. Ann Stat. 2001;29(5):1189-1232",br(),
                        "3. Chen T, Guestrin C. XGBoost: A scalable tree boosting system. in Proc the 22nd ACM SIGKDD Int Conf KDD. ACM. 2016;785-794", br(),
                        "4. Lundberg SM, Lee SI. A unified approach to interpreting model predictions. in Proc Adv Neural Inf Process Syst. 2017;4765-4774.",br(),
                        "5. Aitchison J. The statistical analysis of compositional data. J R Stat Soc B. 1982;44(2):139-77")
  XGB_REFERENCE_RC = p("1. Kim J, Jang H, Koh H. MiMultiCat: A unified cloud platform for the analysis of microbiome data with multi-categorical responses. (Under review)", br(),
                       "2. Friedman JH. Greedy function approximation: A gradient boosting machine. Ann Stat. 2001;29(5):1189-1232",br(),
                       "3. Chen T, Guestrin C. XGBoost: A scalable tree boosting system. in Proc the 22nd ACM SIGKDD Int Conf KDD. ACM. 2016;785-794", br(),
                       "4. Lundberg SM, Lee SI. A unified approach to interpreting model predictions. in Proc Adv Neural Inf Process Syst. 2017;4765-4774.",br(),
                       "5. Sanders HL. Marine benthic diversity: A comparative study. Am Nat. 1968;102:243-282.")
}

# UI ---------------------------------------------------------------------------
{
  ui = dashboardPage(
    title = "MiMultiCat",
    dashboardHeader(title = span(TITLE, style = "float:left;font-size: 20px"), titleWidth = "100%"),
    dashboardSidebar(
      tags$script(JS("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';")),
      sidebarMenu(id = "side_menu",
                  menuItem("Home", tabName = "home"),
                  menuItem("Data Processing",
                           menuSubItem("Data Input", tabName = "step1", 
                                       icon = fontawesome::fa("upload", margin_left = "0.3em", margin_right = "0.1em")),
                           menuSubItem("Quality Control", tabName = "step2", 
                                       icon = fontawesome::fa("chart-bar", margin_left = "0.3em")),
                           menuSubItem(span(span("Diversity Calculation /", br()), 
                                            span("Data Transformation", style = "margin-left: 23px")), tabName = "divdtCalculation", 
                                       icon = fontawesome::fa("calculator", margin_left = "0.3em", margin_right = "0.25em"))),
                  menuItem("Association",
                           menuSubItem("Alpha Diversity", tabName = "alphaDiv", 
                                       icon = fontawesome::fa("font", margin_left = "0.2em", margin_right = "0.1em")),
                           menuSubItem("Beta Diversity", tabName = "betaDiv", 
                                       icon = fontawesome::fa("bold", margin_left = "0.3em", margin_right = "0.15em")),
                           menuSubItem("Taxonomic Analysis", tabName = "taxa", 
                                       icon = fontawesome::fa("diagram-project"))),
                  menuItem("Prediction",
                           menuSubItem("Random Forest", tabName = "rf", 
                                       icon = fontawesome::fa("network-wired")),
                           menuSubItem("Gradient Boosting", tabName = "xgb", 
                                       icon = fontawesome::fa("diagram-project")))
                  )
      ),
    
    dashboardBody(
      
      # THEME -----
      shinyDashboardThemeDIY(
        
        # general
        appFontFamily = "Arial"
        ,appFontColor = "rgb(0,0,0)"
        ,primaryFontColor = "rgb(0,0,0)"
        ,infoFontColor = "rgb(0,0,0)"
        ,successFontColor = "rgb(0,0,0)"
        ,warningFontColor = "rgb(0,0,0)"
        ,dangerFontColor = "rgb(0,0,0)"
        ,bodyBackColor = "rgb(255,255,255)"
        
        # header
        ,logoBackColor = "rgb(65,186,81)"
        
        ,headerButtonBackColor = "rgb(65,186,81)"
        ,headerButtonIconColor = "rgb(65,186,81)"
        ,headerButtonBackColorHover = "rgb(65,186,81)"
        ,headerButtonIconColorHover = "rgb(0,0,0)"
        
        ,headerBackColor = "rgb(65,186,81)"
        ,headerBoxShadowColor = "#aaaaaa"
        ,headerBoxShadowSize = "0px 0px 0px"
        
        # sidebar
        ,sidebarBackColor = "rgb(24,31,41)"
        ,sidebarPadding = 0
        
        ,sidebarMenuBackColor = "transparent"
        ,sidebarMenuPadding = 0
        ,sidebarMenuBorderRadius = 0
        
        ,sidebarShadowRadius = ""
        ,sidebarShadowColor = "0px 0px 0px"
        
        ,sidebarUserTextColor = "rgb(24,31,41)"
        
        ,sidebarSearchBackColor = "rgb(255, 255, 255)"
        ,sidebarSearchIconColor = "rgb(24,31,41)"
        ,sidebarSearchBorderColor = "rgb(24,31,41)"
        
        ,sidebarTabTextColor = "rgb(210,210,210)"
        ,sidebarTabTextSize = 14
        ,sidebarTabBorderStyle = "none"
        ,sidebarTabBorderColor = "none"
        ,sidebarTabBorderWidth = 0
        
        ,sidebarTabBackColorSelected = "rgb(45,52,63)"
        ,sidebarTabTextColorSelected = "rgb(252,255,255)"
        ,sidebarTabRadiusSelected = "0px"
        
        ,sidebarTabBackColorHover = "rgb(67,75,86)"
        ,sidebarTabTextColorHover = "rgb(252,255,255)"
        ,sidebarTabBorderStyleHover = "none"
        ,sidebarTabBorderColorHover = "none"
        ,sidebarTabBorderWidthHover = 0
        ,sidebarTabRadiusHover = "0px"
        
        # boxes
        ,boxBackColor = "rgb(245,245,245)"
        ,boxBorderRadius = 3
        ,boxShadowSize = "0px 0px 0px"
        ,boxShadowColor = "rgba(0,0,0,0)"
        ,boxTitleSize = 16
        ,boxDefaultColor = "rgb(210,214,220)"
        ,boxPrimaryColor = "rgb(35, 49, 64)"
        ,boxInfoColor = "rgb(65,186,81)"
        ,boxSuccessColor = "rgb(102,199,115)"
        ,boxWarningColor = "rgb(244,156,104)"
        ,boxDangerColor = "rgb(255,88,55)"
        
        ,tabBoxTabColor = "rgb(255,255,255)"
        ,tabBoxTabTextSize = 14
        ,tabBoxTabTextColor = "rgb(0,0,0)"
        ,tabBoxTabTextColorSelected = "rgb(35, 49, 64)"
        ,tabBoxBackColor = "rgb(255,255,255)"
        ,tabBoxHighlightColor = "rgb(65,186,81)"
        ,tabBoxBorderRadius = 0
        
        # inputs
        ,buttonBackColor = "rgb(245,245,245)"
        ,buttonTextColor = "rgb(0,0,0)"
        ,buttonBorderColor = "rgb(24,31,41)"
        ,buttonBorderRadius = 3
        
        ,buttonBackColorHover = "rgb(227,227,227)"
        ,buttonTextColorHover = "rgb(100,100,100)"
        ,buttonBorderColorHover = "rgb(200,200,200)"
        
        ,textboxBackColor = "rgb(255,255,255)"
        ,textboxBorderColor = "rgb(200,200,200)"
        ,textboxBorderRadius = 0
        ,textboxBackColorSelect = "rgb(245,245,245)"
        ,textboxBorderColorSelect = "rgb(200,200,200)"
        
        # tables
        ,tableBackColor = "rgb(255, 255, 255)"
        ,tableBorderColor = "rgb(245, 245, 245)"
        ,tableBorderTopSize = 1
        ,tableBorderRowSize = 1
        
      ),
      
      # STYLE -----
      
      ## CONTENTS -----
      tags$head(tags$style(HTML(".content { padding-top: 2px;}"))),
      
      ## PROGRESS BAR -----
      tags$head(tags$style(HTML('.progress-bar {background-color: (102,199,115);}'))),
      
      ## PRETTY RADIO BUTTON -----
      tags$head(tags$style(HTML('
      .pretty input:checked~.state.p-primary label:after, .pretty.p-toggle .state.p-primary label:after {
    background-color: #66C773!important;}'))),
      
      ## SLIDER -----
      setSliderColor(rep("#66C773", 100), seq(1, 100)),
      chooseSliderSkin("Flat"),
      
      tags$script(src = "fileInput_text.js"),
      useShinyjs(),
      tabItems(
        
        # Home -----
        tabItem(tabName = "home",
                div(id = "homepage", br(), HOME_COMMENT, 
                    div(tags$img(src="MiMultiCat_Home_Img.png", height = 720, width = 600), style = "text-align: center"),
                    HOME_COMMENT2, HOME_COMMENT3, HOME_COMMENT4)),
        
        ## 0. DATA INPUT -----
        
        tabItem(tabName = "step1", br(),
                fluidRow(column(width = 6,
                                box(
                                  width = NULL, status = "info", solidHeader = TRUE,
                                  title = strong("Data Input", style = "color:white"),
                                  selectInput("inputOption", h4(strong("Data type")), 
                                              c("Choose one" = "", "Phyloseq", "Individual Data"), width = '30%'),
                                  div(id = "optionsInfo", 
                                      tags$p("You can choose phyloseq or individual data.", style = "font-size:11pt"), 
                                      tags$p("", style = "margin-bottom:-8px"), style = "margin-top: -15px"),
                                  uiOutput("moreOptions")
                                  )
                                ),
                         column(width = 6, style='padding-left:0px', 
                                uiOutput("addDownloadinfo")
                                )
                         )
                ),
        
        ## 1-1. QC ----
        
        tabItem(tabName = "step2", br(),
            fluidRow(column(width = 3,  style = "padding-left:+15px",
                            # Quality Control
                            box(
                              width = NULL, status = "info", solidHeader = TRUE,
                              title = strong("Quality Control", style = "color:white"),
                              textInput("kingdom", h4(strong("Kingdom")), value = "Bacteria"),
                              QC_KINGDOM_COMMENT,
                              
                              tags$style(type = 'text/css', '#slider1 .irs-grid-text {font-size: 1px}'),
                              tags$style(type = 'text/css', '#slider2 .irs-grid-text {font-size: 1px}'),
                              
                              sliderInput("slider1", h4(strong("Library size")), 
                                          min=0, max=10000, value = 3000, step = 1000),
                              QC_LIBRARY_SIZE_COMMENT1,
                              QC_LIBRARY_SIZE_COMMENT2,
                              
                              sliderInput("slider2", h4(strong("Mean proportion")), 
                                          min = 0, max = 0.1, value = 0.02, step = 0.001,  post  = " %"),
                              QC_MEAN_PROP_COMMENT1,
                              QC_MEAN_PROP_COMMENT2,
                              
                              br(),
                              p(" ", style = "margin-bottom: -20px;"),
                              
                              h4(strong("Errors in taxonomic names")),
                              textInput("rem.str", label = "Complete match", value = ""),
                              QC_TAXA_NAME_COMMENT1,
                              
                              textInput("part.rem.str", label = "Partial match", value = ""),
                              QC_TAXA_NAME_COMMENT2,
                              
                              actionButton("run", (strong("Run!")), class = "btn-info"), 
                              p(" ", style = "margin-bottom: +10px;"), 
                              p(strong("Attention:"),"You have to click this Run button to perform diversity calculation, 
                                data transformation and further analyses.", style = "margin-bottom:-10px"), br()),
                            
                            uiOutput("moreControls"),
                            uiOutput("qc_reference")),
                     
                     column(width = 9, style = "padding-left:+10px",
                            box(
                              width = NULL, status = "info", solidHeader = TRUE,
                              fluidRow(width = 12,
                                       status = "info", solidHeader = TRUE, 
                                       valueBoxOutput("sample_Size", width = 3),
                                       valueBoxOutput("OTUs_Size", width = 3),
                                       valueBoxOutput("phyla", width = 3),
                                       valueBoxOutput("classes", width = 3)),
                              
                              fluidRow(width = 12, 
                                       status = "info", solidHeader = TRUE,
                                       valueBoxOutput("orders", width = 3),
                                       valueBoxOutput("families", width = 3),
                                       valueBoxOutput("genera", width = 3),
                                       valueBoxOutput("species", width = 3)),
                              
                              fluidRow(style = "position:relative",
                                       tabBox(width = 6, title = strong("Library Size", style = "color:black"), 
                                              tabPanel("Histogram",
                                                       plotlyOutput("hist"),
                                                       sliderInput("binwidth", "# of Bins:",
                                                                   min = 0, max = 100, value = 50, width = "100%")),
                                              tabPanel("Box Plot", 
                                                       plotlyOutput("boxplot"))),
                                       tabBox(width = 6, title = strong("Mean Proportion", style = "color:black"), 
                                              tabPanel("Histogram",
                                                       plotlyOutput("hist2"),
                                                       sliderInput("binwidth2", "# of Bins:",
                                                                   min = 0, max = 100, value = 50, width = "100%")),
                                              tabPanel("Box Plot",
                                                       plotlyOutput("boxplot2")))))))
            ),
        
        ## 1-2. DC & DT -----
        
        tabItem(tabName = "divdtCalculation", br(),
                fluidRow(column(width = 6, style = "padding-left:+15px",
                                box(title = strong("Diversity Calculation & Data Transformation", style = "color:white"), 
                                    width = NULL, status = "info", solidHeader = TRUE, 
                                    
                                    strong(p("Diversity Calculation", style = "font-size:12pt")),
                                    ALPHA_COMMENT, 
                                    BETA_COMMENT, 
                                    
                                    strong(p("Data Transformation", style = "font-size:12pt")),
                                    DATA_TRANSFORM_COMMENT, 
                                    actionButton("divdtCalcRun", (strong("Run!")), class = "btn-info")),
                                uiOutput("divdtDownload")),
                         
                         column(width = 6, style='padding-left:0px',
                                box(title = strong("References", style = "color:white"), 
                                    width = NULL, status = "info", solidHeader = TRUE,
                                    strong(p("Alpha Diversity", style = "font-size:12pt")), ALPHA_REFERENCES,
                                    strong(p("Beta Diversity", style = "font-size:12pt")), BETA_REFERENCES,
                                    strong(p("Data Transformation", style = "font-size:12pt")), DATA_TRANSFORM_REFERENCE)))
              ),
        
        ## 2. ALPHA ------
        
        tabItem(tabName = "alphaDiv", br(),
                tabPanel(
                  title = NULL,
                  sidebarLayout( 
                    position = "left",
                    sidebarPanel(width = 3,
                                 shinyjs::hidden(
                                   uiOutput("alpha_primvars"),
                                   uiOutput("alpha_var_type"),
                                   uiOutput("alpha_primvars_ref"),
                                   uiOutput("alpha_primvars_rename"),
                                   uiOutput("alpha_covariate"),
                                   uiOutput("alpha_method"),
                                   uiOutput("alpha_run"),
                                   uiOutput("alpha_downloadTable"),
                                   uiOutput("alpha_references"))),
                    
                    mainPanel(width = 9,
                              fluidPage(width = NULL, uiOutput("alpha_display_results"), uiOutput("alpha_display_pair")),
                              uiOutput("barPanel"), br(), br())))),
        
        ## 3. BETA ------
        
        tabItem(tabName = "betaDiv", br(),
                tabPanel(
                  title = NULL,
                  sidebarLayout( 
                    position = "left",
                    sidebarPanel(width = 3,
                                 shinyjs::hidden(
                                   uiOutput("beta_primvars"),
                                   uiOutput("beta_var_type"),
                                   uiOutput("beta_primvars_ref"),
                                   uiOutput("beta_primvars_rename"),
                                   uiOutput("beta_covariate"),
                                   uiOutput("beta_method"),
                                   uiOutput("beta_run"),
                                   uiOutput("beta_downloadTabUI"),
                                   uiOutput("beta_reference"))),
                    
                    mainPanel(width = 9,
                              fluidPage(width = 12, uiOutput("beta_nom_results")),
                              uiOutput("beta_barPanel"), br(), br())))),
        
        ## 4. TAXA ------
        
        tabItem(tabName = "taxa", br(),
                tabPanel(
                  title = NULL,
                  sidebarLayout( 
                    position = "left",
                    sidebarPanel(width = 3,
                                 shinyjs::hidden(
                                   uiOutput("taxa_primvars"),
                                   uiOutput("taxa_var_type"),
                                   uiOutput("taxa_primvars_ref"),
                                   uiOutput("taxa_primvars_rename"),
                                   uiOutput("taxa_covariate"),
                                   uiOutput("taxa_method"),
                                   uiOutput("taxa_run"),
                                   uiOutput("taxa_downloadTable"),
                                   uiOutput("taxa_references"))),
                    
                    mainPanel(width = 9,
                              fluidPage(width = NULL, 
                                        uiOutput("taxa_display_global"),
                                        uiOutput("taxa_display_forest"),
                                        div(id = "taxa_display_area", style='height:550px;overflow-y: scroll;', uiOutput("taxa_display")), 
                                        uiOutput("taxa_display_dend")),
                              br(), 
                              br(), 
                              br(), 
                              br(), 
                              uiOutput("taxa_barPanel"))))),
        
        ## 5-1. RF ------
        
        tabItem(tabName = "rf", br(),
                sidebarLayout(
                  position = "left",
                  sidebarPanel(width = 3,
                               shinyjs::hidden(
                                 uiOutput("rf_nom_data_input"),
                                 uiOutput("rf_nom_data_input_opt"),
                                 uiOutput("rf_nom_train_setting"),
                                 uiOutput("rf_nom_downloadTabUI"),
                                 uiOutput("rf_nom_reference")
                                 )),
                  mainPanel(width = 9,
                            fluidRow(width = 12, uiOutput("rf_nom_results"))))
                ),
        
        ## 5-2. XGB ------
        
        tabItem(tabName = "xgb", br(),
                sidebarLayout(
                  position = "left",
                  sidebarPanel(width = 3,
                               shinyjs::hidden(
                                 uiOutput("xgb_nom_data_input"),
                                 uiOutput("xgb_nom_data_input_opt"),
                                 uiOutput("xgb_nom_train_setting"),
                                 uiOutput("xgb_nom_downloadTabUI"),
                                 uiOutput("xgb_nom_reference")
                                 )),
                  mainPanel(width = 9,
                            fluidRow(width = 12, uiOutput("xgb_nom_results"))))
                )
      )
    )
  )
}

# SERVER -----------------------------------------------------------------------
server = function(input, output, session){
  options(shiny.maxRequestSize=30*1024^2)
  
  ## REACTIVE VALUES -------
  
  # User input & Quality control
  infile = reactiveValues(biom = NULL, qc_biom = NULL, rare_biom = NULL, na_omit_biom = NULL)
  chooseData = reactiveValues(sam.dat = NULL, mon.sin.rev.bin.con = NULL, prim_vars = NULL, alpha.div = NULL,
                              alpha.div.rare = NULL, alpha.div.qc = NULL, taxa.out = NULL, tax.tab = NULL)
  taxa.results = reactiveValues(bin.var = NULL, cov.var = NULL, id.var = NULL, taxa = NULL, taxa.bin.sum.out = NULL,
                                con.var = NULL, taxa.con.sum.out = NULL, lib.size = NULL)
  
  # Diversity calculation
  ds.Ks <- reactiveValues(res = NULL)
  
  # Reactive responses & analysis
  alpha_response <- reactiveValues(cat = NULL, cat.name = NULL)
  beta_response <- reactiveValues(cat = NULL, cat.name = NULL)
  taxa_response <- reactiveValues(cat = NULL, cat.name = NULL)
  rf_response <- reactiveValues(cat = NULL)
  xgb_response <- reactiveValues(cat = NULL)
  xgb.model.input.cla <- reactiveValues(eval = NULL)
  
  # Histogram colors on data processing tab
  rcol = reactiveValues(selected = "lightblue")
  
  ## EX DATA ------

  env.v1 <- new.env()
  nm.v1 <- load(file = "Data/biom.Rdata", env.v1)[1]
  biom.v1 <- env.v1[[nm.v1]]
  
  otu.tab.v1 <- otu_table(biom.v1)
  tax.tab.v1 <- tax_table(biom.v1)
  sam.dat.v1 <- sample_data(biom.v1)
  tree.v1 <- phy_tree(biom.v1)
  
  output$biom.v1 <- downloadHandler(
    filename = function() {
      paste("biom.v1.Rdata", sep = "")
    },
    content = function(file1) {
      save(biom.v1, file = file1)
    })
  output$download.biom.v1 <- downloadHandler(
    filename = function() {
      paste("biom",".zip", sep = "")
    },
    content <- function(fname) {
      temp <- setwd(tempdir())
      on.exit(setwd(temp))
      dataFiles = c("otu.tab.txt", "tax.tab.txt", "sam.dat.txt", "tree.tre")
      write.table(otu.tab.v1, "otu.tab.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.table(tax.tab.v1, "tax.tab.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.table(data.frame(sam.dat.v1), "sam.dat.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
      write.tree(tree.v1, "tree.tre")
      zip(zipfile=fname, files=dataFiles)
    })
  
  ## 0. DATA INPUT -----------
  observeEvent(input$inputOption,{
    observe({
      if (input$inputOption == "Phyloseq") {
        shinyjs::hide(id = "optionsInfo")
        output$moreOptions <- renderUI({
          tagList(
            tags$style("
                       .btn-file {
                       border-top-left-radius: 5px !important; border-bottom-left-radius: 5px !important; border-left-style: solid !important; border-left-width: 1px !important;
                       border-top-right-radius: 0px !important; border-bottom-right-radius: 0px !important; border-right-width: 0px !important;
                       }"
            ),
            fileInput("phyloseqData", strong("Please upload your 'phyloseq' data (.Rdata, .rds)", style = "color:black"), 
                      accept = c(".Rdata", ".rds"), width = '80%'), div(style = "margin-top: -15px"),
            actionButton('Load_Phyloseq_Data', 'Upload', class = "btn-info"), 
            p(" ", style = "margin-bottom: +10px;"), 
            p(strong("Attention:"), "You have to click this Upload button to perform following data processing and downstream data analyses."),br(),
            shinyjs::hidden(
              shiny::div(id = "phyloseqUpload_error",
                         shiny::tags$p("Please upload a Rdata file!!",
                                       style = "color: red; font-weight: bold; padding-top: 5px;", class = "text-center"))
            ),
            INPUT_PHYLOSEQ_COMMENT1,
            p("", style = "margin-bottom:-8px")
          )
        })
        
        output$addDownloadinfo <- renderUI({
          tagList(
            box(title = strong("Example Data", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                downloadButton("biom.v1", "Example Data", width = '30%', style = "color:black; background-color: red2"),
                br(),br(),
                INPUT_PHYLOSEQ_COMMENT2,
                p("", style = "margin-bottom:-8px")
            )
          )
        })
      }
      else if (input$inputOption == "Individual Data") {
        shinyjs::hide(id = "optionsInfo")
        output$moreOptions <- renderUI({
          tagList(
            tags$style("
                       .btn-file {
                       border-top-left-radius: 5px !important; border-bottom-left-radius: 5px !important; border-left-style: solid !important; border-left-width: 1px !important;
                       border-top-right-radius: 0px !important; border-bottom-right-radius: 0px !important; border-right-width: 0px !important;
                       }"
            ),
            fileInput("otuTable", strong("Please upload your feature (OTU or ASV) table (.txt, .csv, .biom)", style = "color:black"), 
                      accept = c(".txt", ".csv", ".biom"), width = '80%'), div(style = "margin-top: -15px"),
            fileInput("taxTable", strong("Please upload your taxonomic table (.txt, .tsv)", style = "color:black"), 
                      accept = c(".txt", ".tsv"), width = '80%'), div(style = "margin-top: -15px"),
            fileInput("samData", strong("Please upload your metadata/sample information (.txt, .csv)", style = "color:black"), 
                      accept = c(".txt", ".csv"), width = '80%'), div(style = "margin-top: -15px"),
            fileInput("tree", strong("Please upload your phylogenetic tree (.tre, .nwk)", style = "color:black"), 
                      accept = c(".tre", ".nwk"), width = '80%'), div(style = "margin-top: -15px"),
            actionButton('Load_Individual_Data', 'Upload', class = "btn-info"), br(),br(),
            shinyjs::hidden(
              shiny::div(id = "textfilesUpload_error",
                         shiny::tags$p("Please upload txt and tre files!!",
                                       style = "color: red; font-weight: bold; padding-top: 5px;", class = "text-center"))
            ),
            INPUT_INDIVIDUAL_DATA_COMMENT, 
            p("", style = "margin-bottom:-8px")
          )
        })
        
        output$addDownloadinfo <- renderUI({
          tagList(
            box(title = strong("Example Data", style = "color:white"), 
                width = NULL, status = "info", solidHeader = TRUE,
                downloadButton("download.biom.v1", "Example Data", width = '30%', style = "color:black; background-color: red2"),
                br(),br(),
                INPUT_INDIVIDUAL_DATA_COMMENT2,
                p("", style = "margin-bottom:-8px")
            )
          )
        })
      }
    })
    
  }, ignoreInit = TRUE, once = TRUE, ignoreNULL = TRUE)
  
  observe({
    toggleState("Load_Phyloseq_Data", !is.null(input$phyloseqData))
    toggleState("Load_Individual_Data", 
                !(is.null(input$otuTable) | is.null(input$taxTable) | is.null(input$samData)))
    toggleState("run", !is.null(infile$biom))
    toggleState("skip", !is.null(infile$biom))
    toggleState("slider1", !is.null(infile$biom))
    toggleState("slider2", !is.null(infile$biom))
    toggleState("kingdom", !is.null(infile$biom))
    toggleState("binwidth", !is.null(infile$biom))
    toggleState("binwidth2", !is.null(infile$biom))
    
    toggleState("divdtCalcRun", !is.null(infile$rare_biom))
    toggleState("datTransRun", !is.null(infile$rare_biom))
  })
  
  observeEvent(input$Load_Phyloseq_Data, {
    
    if (!is.null(input$phyloseqData)) {
      dataInfile  = reactive({
        phyloseq.data = input$phyloseqData
        ext <- tools::file_ext(phyloseq.data$datapath)
        
        req(phyloseq.data)
        if (ext == "Rdata") {
          phyloseq.dataPath = phyloseq.data$datapath
          e = new.env()
          name <- load(phyloseq.dataPath, envir = e)
          data <- e[[name]]
          
          if (sum(sapply(sample_data(data),is.factor))!=0) {
            sample_data(data)[,which(sapply(sample_data(data), is.factor))] = lapply(sample_data(data)[,which(sapply(sample_data(data), is.factor))], as.character)
          }
          
          colnames(tax_table(data)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
          
          if (sum(colnames(otu_table(data)) %in% rownames(sample_data(data))) < sum(rownames(otu_table(data)) %in% rownames(sample_data(data)))) {
            otu_table(data) = t(otu_table(data))
          }
          
          return(data)
        } 
        else if (ext == "rds") {
          phyloseq.dataPath = phyloseq.data$datapath
          data <- readRDS(phyloseq.dataPath)
          
          if (sum(sapply(sample_data(data),is.factor))!=0) {
            sample_data(data)[,which(sapply(sample_data(data), is.factor))] = lapply(sample_data(data)[,which(sapply(sample_data(data), is.factor))], as.character)
          }
          colnames(tax_table(data)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
          
          if (sum(colnames(otu_table(data)) %in% rownames(sample_data(data))) < sum(rownames(otu_table(data)) %in% rownames(sample_data(data)))) {
            otu_table(data) = t(otu_table(data))
          }
          
          return(data)
        } 
        else {
          shinyjs::toggle(id = "phyloseqUpload_error", anim = TRUE, time = 1, animType = "fade")
          shinyjs::delay(5000, shinyjs::toggle(id = "phyloseqUpload_error", anim = TRUE, time = 1, animType = "fade"))
          return(NULL)
        }
      })
    } else {
      return(NULL)
    }
    
    if (is.null(dataInfile)) {
      infile$biom <- NULL
      infile$qc_biom <- NULL
      infile$rare_biom = NULL
    } 
    else {
      infile$biom <- dataInfile()
      infile$qc_biom <- dataInfile()
      infile$rare_biom = NULL
    }
    
    updateTabsetPanel(session, "side_menu",
                      selected = "step2")
    rcol$selected = "lightblue"
    
    if (!is.null(infile$biom)) QC$resume()
  })
  observeEvent(input$Load_Individual_Data, {
    shinyjs::disable("Load_Individual_Data")
    
    withProgress(
      message = 'Calculation in progress',
      detail = 'This may take a while...', value = 0, {
        incProgress(3/10, message = "File Check")
        if (!is.null(input$otuTable) & !is.null(input$taxTable) & !is.null(input$samData) & !is.null(input$tree)) {
          dataInfile  = reactive({
            otu.table = input$otuTable
            ext1 <- tools::file_ext(otu.table$datapath)
            
            tax.table = input$taxTable
            ext2 <- tools::file_ext(tax.table$datapath)
            
            sam.data = input$samData
            ext3 <- tools::file_ext(sam.data$datapath)
            
            tree.data = input$tree
            ext4 <- tools::file_ext(tree.data$datapath)
            
            req(otu.table, tax.table, sam.data, tree.data)
            if ((ext1 == "txt"| ext1 == "csv" | ext1 == "biom") & (ext2 == "txt" | ext2 == "tsv") &
                (ext3 == "txt" | ext3 == "csv") & (ext4 == "tre" | ext4 == "nwk")) {
              otu.table.path = otu.table$datapath
              tax.table.path = tax.table$datapath
              sam.data.path = sam.data$datapath
              tree.data.path = tree.data$datapath
              
              if (ext1 == "txt") {
                otu.tab <- read.table(otu.table.path, header=TRUE, check.names = FALSE, sep = "\t")
              } 
              else if (ext1 == "csv") {
                otu.tab <- read.csv(otu.table.path, check.names = FALSE)
                rownames(otu.tab) = otu.tab[,1];otu.tab = otu.tab[,-1]
              }
              else if (ext1 == "biom") {
                biom <- read_biom(otu.table.path)
                otu.tab <- as.matrix(biom_data(biom))
              }
              
              if (ext2 == "txt") {
                tax.tab <- read.table(tax.table.path, header=TRUE, check.names = FALSE, sep = "\t")
              } 
              else if (ext2 == "tsv") {
                tax.tab <- read.table(tax.table.path, header=TRUE, sep="\t")
                tax.tab = preprocess.tax.tab(tax.tab)
              }
              
              if (ext3 == "txt") {
                sam.dat <- read.table(sam.data.path, header=TRUE, check.names = FALSE, sep = "\t")
              } 
              else if (ext3 == "csv") {
                sam.dat <- read.csv(sam.data.path, check.names = FALSE)
                rownames(sam.dat) = sam.dat[,1]
                sam.dat = sam.dat[,-1]
              }
              
              if (ext4 == "tre") {
                tree <- read.tree(file = tree.data.path)
              } else if (ext4 == "nwk") {
                tree <- read.tree(file = tree.data.path) 
              }
              
              otu.tab <- otu_table(otu.tab, taxa_are_rows = TRUE)
              tax.tab <- tax_table(as.matrix(tax.tab))
              sam.dat <- sample_data(sam.dat)
              tree <- phy_tree(tree)
              
              if (sum(colnames(otu.tab) %in% rownames(sam.dat)) < sum(rownames(otu.tab) %in% rownames(sam.dat))) {
                otu.tab = t(otu.tab)
              }
              
              incProgress(3/10, message = "Validating")
              shiny::validate(
                if (biom.check.samples(otu.tab, sam.dat)) {
                  if (biom.check.otu(otu.tab, tax.tab, tree)) {
                    showNotification(h4("Error: There is no common samples among OTU/feature table and Sample Data. And
                                        there is no common OTUs among OTU/feature table and taxonomic table."),
                                     type = "error")
                  } 
                  else {
                    showNotification(h4("Error: There is no common samples among OTU/feature table and Sample Data"),
                                     type = "error")
                  }
                } 
                else if (biom.check.otu(otu.tab, tax.tab, tree)) {
                  showNotification(h4("Error: There is no common OTUs among OTU/feature table and taxonomic table."),
                                   type = "error")
                } 
                else {
                  NULL
                }
              )
              
              incProgress(1/10, message = "Merging")
              biomData <- merge_phyloseq(otu.tab, tax.tab, sam.dat, tree)
              return(biomData)
            } 
            else {
              shinyjs::toggle(id = "textfilesUpload_error", anim = TRUE, time = 1, animType = "fade")
              shinyjs::delay(5000, shinyjs::toggle(id = "textfilesUpload_error", anim = TRUE, time = 1, animType = "fade"))
              return(NULL)
            }
          })
        } 
        else {
          return(NULL)
        }
        
        if (is.null(dataInfile)) {
          infile$biom <- NULL
          infile$qc_biom <- NULL
          infile$rare_biom = NULL
        } 
        else {
          infile$biom <- dataInfile()
          infile$qc_biom <- dataInfile()
          infile$rare_biom = NULL
        }
        
        updateTabsetPanel(session, "side_menu",
                          selected = "step2")
        rcol$selected = "lightblue"
        
        if (!is.null(infile$biom)) QC$resume()
      })
    shinyjs::enable("Load_Individual_Data")
  })
  
  ## 1-1. QC -----------
  # This reactive expression stores the input data from either the individual data or phyloseq data
  QC = observe(suspended = T,{
    taxa.results$lib.size <- lib.size.func(infile$biom)$lib.size
    
    # Plots graphs using example infile data
    output$hist <- renderPlotly({
      lib_size = lib.size.func(infile$qc_biom)$lib.size
      plot_ly(x = ~lib_size, nbinsx = input$binwidth,
              type = "histogram",
              marker = list(color = rcol$selected, line = list(color = "black", width = 2))) %>%
        layout(
          yaxis = list(title = "Frequency", zeroline = FALSE),
          xaxis = list(title = "Library Size", zeroline = FALSE))
    })
    
    output$hist2 <- renderPlotly({
      mean_prop = mean.prop.func(infile$qc_biom)$mean.prop
      plot_ly(x = ~mean_prop, nbinsx = input$binwidth2,
              type = "histogram",
              marker = list(color = rcol$selected, line = list(color = "black", width = 2))) %>%
        layout(
          yaxis = list(title = "Frequency", zeroline = FALSE),
          xaxis = list(title = "Mean Proportion", zeroline = FALSE))
    })
    
    output$boxplot<- renderPlotly({
      lib_size = lib.size.func(infile$qc_biom)$lib.size
      
      plot_ly(x = ~lib_size, type = "box", notched=TRUE, name = "Library Size",
              color = ~"lib_size", colors = rcol$selected, line = list(color = 'black'))%>%
        layout(
          yaxis = list(title = "", zeroline = FALSE),
          xaxis = list(title = "", zeroline = FALSE), showlegend = FALSE)
    })
    
    output$boxplot2<- renderPlotly({
      mean_prop = mean.prop.func(infile$qc_biom)$mean.prop
      
      plot_ly(x = ~mean_prop, type = "box", notched=TRUE, name = "Mean Proportion",
              color = ~"mean_prop", colors = rcol$selected, line = list(color = 'black'))%>%
        layout(
          yaxis = list(title = "", zeroline = FALSE),
          xaxis = list(title = "", zeroline = FALSE), showlegend = FALSE)
    })
    
    ## Number of Taxonomic Rank for biom before QC
    num_tax.rank = reactive({
      tax.tab = tax_table(infile$qc_biom)
      num.tax.rank(tax.tab)
    })
    
    ## Fills value boxes using example biom data
    output$sample_Size <- renderValueBox({
      valueBox(
        value = tags$p(paste0(lib.size.func(infile$qc_biom)$num.sams), style = "font-size: 75%;"),
        "Sample Size", icon = icon("user-circle"), color = "fuchsia")
    })
    
    output$OTUs_Size <- renderValueBox({
      valueBox(
        value = tags$p(paste0(lib.size.func(infile$qc_biom)$num.otus), style = "font-size: 75%;"),
        "Number of Features", icon = icon("dna"), color = "aqua")
    })
    
    output$phyla <- renderValueBox({
      num.phyla = num_tax.rank()[1]
      valueBox(
        value = tags$p(paste0(num.phyla), style = "font-size: 75%;"),
        "Number of Phyla", icon = icon("sitemap"), color = "orange")
    })
    
    output$classes <- renderValueBox({
      num.classes = num_tax.rank()[2]
      valueBox(
        value = tags$p(paste0(num.classes), style = "font-size: 75%;"),
        "Number of Classes", icon = icon("sitemap"), color = "purple")
    })
    
    output$orders <- renderValueBox({
      num.orders = num_tax.rank()[3]
      valueBox(
        value = tags$p(paste0(num.orders), style = "font-size: 75%;"),
        "Number of Orders", icon = icon("sitemap"), color = "blue")
    })
    
    output$families <- renderValueBox({
      num.families = num_tax.rank()[4]
      valueBox(
        value = tags$p(paste0(num.families), style = "font-size: 75%;"),
        "Number of Families", icon = icon("sitemap"), color = "red")
    })
    
    output$genera <- renderValueBox({
      num.genera = num_tax.rank()[5]
      valueBox(
        value = tags$p(paste0(num.genera), style = "font-size: 75%;"),
        "Number of Genera", icon = icon("sitemap"), color = "lime")
    })
    
    output$species <- renderValueBox({
      num.species = num_tax.rank()[6]
      valueBox(
        value = tags$p(paste0(num.species), style = "font-size: 75%;"),
        "Number of Species", icon = icon("sitemap"), color = "teal" )
    })
    
    ## This event handler checks whether there is an input file and updates the slider options accordingly
    maxi.slider1 = as.numeric(lib.size.func(infile$qc_biom)$lib.size.sum["3rd quartile"])
    # max.mean.prop = as.numeric(mean.prop.func(infile$qc_biom)$mean.prop.sum["3rd quartile"])
    # maxi.slider2 = round(max.mean.prop, digits = 6)
    # if (maxi.slider2 < 2e-05) {
    #   maxi.slider2 = 2e-05
    # }
    # 
    updateSliderInput(session, "slider1", min = 0, max = round(maxi.slider1,-3))
    # updateSliderInput(session, "slider2", min = 0, max = maxi.slider2*100)
  })
  
  ## 2. ALPHA ---------------------------
  
  output$alpha_primvars <- renderUI({
    tagList(
      p(" ", style = "margin-top: 25px;"),
      h4(strong("Response variable", style = "color:black")),
      p("E.g., a multi-categorical variable on the host's health or disease status", style = "font-size:10pt"),
      p(" ", style = "margin-bottom: +15px;"),
      selectInput("alpha_prim", label = NULL, 
                  choices = c("Choose one" = "",  bmc.col.check(chooseData$sam.dat, type = "Multinomial")), 
                  width = '70%')
    )
  })
  
  output$alpha_var_type <- renderUI({
    tagList(
      prettyRadioButtons("alpha_var_type", label = h4(strong("Variable type", style = "color:black")), animation = "jelly",
                         c("Nominal", "Ordinal"), icon = icon("check"), selected = "Nominal",width = '70%')
    )
  })
  
  observeEvent(input$alpha_prim, {
    
    if(input$alpha_prim != ""){
      
      shinyjs::show("alpha_primvars_ref")
      shinyjs::show("alpha_primvars_rename")
      
      observeEvent(input$alpha_var_type, {
        
        if(input$alpha_var_type == "Nominal"){
          
          alpha_response$cat = category.names(chooseData$sam.dat, input$alpha_prim)
          alpha_nom_num_cat = length(alpha_response$cat)
          
          output$alpha_primvars_rename <- renderUI({
            tagList(
              h4(strong("Reorder/rename categories", style = "color:black")),
              p("You can reorder/rename groups/levels of the chosen response variable. MiMultiCat keeps up to 8 characters on graphs.", style = "font-size:10pt"),
              lapply(1:alpha_nom_num_cat, function(i){
                fluidRow(
                  column(7, style = list("padding-right: 0px;"),
                         textInput(paste0("alphaCat_", i), label = (paste0("Group/Level ", i, ": ", alpha_response$cat[i])), value = alpha_response$cat[i], width = '80%')),
                  column(5, style = list("padding-left: 0px;"),
                         selectInput(paste0("alphaCat_ord", i), label = "Level:", choices = c("Reference" = 1, 2:alpha_nom_num_cat), selected = i, width = '70%')))
              }))
          }) 

        }
        else{
          
          alpha_response$cat = category.names(chooseData$sam.dat, input$alpha_prim)
          alpha_nom_num_cat = length(alpha_response$cat)
          
          output$alpha_primvars_rename <- renderUI({
            tagList(
              h4(strong("Reorder/rename categories", style = "color:black")),
              p("You can reorder/rename groups/levels of the chosen response variable. MiMultiCat keeps up to 8 characters on graphs.", style = "font-size:10pt"),
              lapply(1:alpha_nom_num_cat, function(i){
                fluidRow(
                  column(7, style = list("padding-right: 0px;"),
                         textInput(paste0("alphaCat_", i), label = (paste0("Group/Level ", i, ": ", alpha_response$cat[i])), value = alpha_response$cat[i], width = '80%')),
                  column(5, style = list("padding-left: 0px;"),
                         selectInput(paste0("alphaCat_ord", i), label = "Level:", choices = 1:alpha_nom_num_cat, selected = i, width = '70%')))
              }))
          }) 
        }
      })
      
    }
    else{ # input$alpha_prim == ""
      shinyjs::hide("alpha_primvars_ref")
      shinyjs::hide("alpha_primvars_rename")
    }
  })
  
  output$alpha_covariate <- renderUI({
    tagList(
      p(" ", style = "margin-top: 25px;"),
      h4(strong("Covariate(s)", style = "color:black")),
      p("Potential confounders (e.g., age, gender) to be adjusted for.", style = "font-size:10pt"),
      p(" ", style = "margin-bottom: +15px;"),
      prettyRadioButtons("alpha_covariate_yn",label = NULL, icon = icon("check"),
                         animation = "jelly", c("None", "Covariate(s)"), selected = "None", width = '70%'),
      
      shinyjs::hidden(
        shiny::div(id = "alpha_covariate_list", style = "margin-left: 2%",
                   prettyCheckboxGroup("alpha_covariate_options"," Please select covariate(s)",
                                       choices = get.cov.col(chooseData$sam.dat)[!get.cov.col(chooseData$sam.dat) %in% c(input$alpha_prim)], width = '70%')))
    )
  })
  
  observeEvent(input$alpha_covariate_yn, {
    observeEvent(input$alpha_var_type, {
      if(input$alpha_var_type == "Nominal"){
        # Nominal / No covariate(s)
        if(input$alpha_covariate_yn == "None"){
          output$alpha_method <- renderUI({
            tagList(
              selectInput("alpha_method", label = h4(strong("Method", style = "color:black")), 
                          choices = c("ANOVA", "Kruskal-Wallis", "Multinomial Logistic Regression"), width = '98%')
            )
          })
        }
        # Nominal / With covariate(s)
        else{
          output$alpha_method <- renderUI({
            tagList(
              selectInput("alpha_method", label = h4(strong("Method", style = "color:black")), 
                          choices = c("Multinomial Logistic Regression"), width = '98%')
            )
          })
        }
      }
      
      else if(input$alpha_var_type == "Ordinal"){
        # Ordinal
        output$alpha_method <- renderUI({
          tagList(
            selectInput("alpha_method", label = h4(strong("Method", style = "color:black")), 
                        choices = c("Proportional Odds Model"), width = '98%')
          )
        })
      }
    })
  })
  
    observeEvent(input$alpha_covariate_yn, {

      if(input$alpha_covariate_yn == "Covariate(s)"){
        shinyjs::show("alpha_covariate_list")
      }
      else{
        shinyjs::hide("alpha_covariate_list")
      }
    })
  
  output$alpha_run <- renderUI({
    tagList(
      actionButton("alpha_runButton", (strong("Run!")), class = "btn-info")
    )
  })
  
  observe({
    alpha.check.order <- c()
    for(i in 1:length(alpha_response$cat)){
      alpha.check.order <- c(alpha.check.order, input[[paste0("alphaCat_ord",i)]])
    }
    alpha.check.order <- length(unique(alpha.check.order))
    print(alpha.check.order)
    toggleState("alpha_runButton", !(input$alpha_prim == "") 
                & (input$alpha_covariate_yn == "None" | (input$alpha_covariate_yn == "Covariate(s)" & length(input$alpha_covariate_options) != 0)) 
                & (((input$alpha_var_type == "Ordinal") & (alpha.check.order == length(alpha_response$cat))) | ((input$alpha_var_type == "Nominal") & (alpha.check.order == length(alpha_response$cat)))))
  })
  
  ## 3. BETA ---------------------------
  
  output$beta_primvars <- renderUI({
    tagList(
      p(" ", style = "margin-top: 25px;"),
      h4(strong("Response variable", style = "color:black")),
      p("E.g., a multi-categorical variable on the host's health or disease status", style = "font-size:10pt"),
      p(" ", style = "margin-bottom: +15px;"),
      selectInput("beta_prim", label = NULL, 
                  choices = c("Choose one" = "",  bmc.col.check(chooseData$sam.dat, type = "Multinomial")), 
                  width = '70%')
    )
  })
  
  output$beta_var_type <- renderUI({
    tagList(
      prettyRadioButtons("beta_var_type", label = h4(strong("Variable type", style = "color:black")), animation = "jelly",
                         c("Nominal", "Ordinal"), icon = icon("check"), selected = "Nominal",width = '70%')
    )
  })
  
  observeEvent(input$beta_prim, {
    
    if(input$beta_prim != ""){
      
      shinyjs::show("beta_primvars_ref")
      shinyjs::show("beta_primvars_rename")
      
      observeEvent(input$beta_var_type, {
        
        if(input$beta_var_type == "Nominal"){
          
          beta_response$cat = category.names(chooseData$sam.dat, input$beta_prim)
          beta_nom_num_cat = length(beta_response$cat)
          
          output$beta_primvars_rename <- renderUI({
            tagList(
              h4(strong("Reorder/rename categories", style = "color:black")),
              p("You can reorder/rename groups/levels of the chosen response variable. MiMultiCat keeps up to 8 characters on graphs.", style = "font-size:10pt"),
              lapply(1:beta_nom_num_cat, function(i){
                fluidRow(
                  column(7, style = list("padding-right: 0px;"),
                         textInput(paste0("betaCat_", i), label = (paste0("Group/Level ", i, ": ", beta_response$cat[i])), value = beta_response$cat[i], width = '80%')),
                  column(5, style = list("padding-left: 0px;"),
                         selectInput(paste0("betaCat_ord", i), label = "Level:", choices = c("Reference" = 1, 2:beta_nom_num_cat), selected = i, width = '70%')))
              }))
          }) 
          
        }
        else{
          
          beta_response$cat = category.names(chooseData$sam.dat, input$beta_prim)
          beta_nom_num_cat = length(beta_response$cat)
          
          output$beta_primvars_rename <- renderUI({
            tagList(
              h4(strong("Reorder/rename categories", style = "color:black")),
              p("You can reorder/rename groups/levels of the chosen response variable. MiMultiCat keeps up to 8 characters on graphs.", style = "font-size:10pt"),
              lapply(1:beta_nom_num_cat, function(i){
                fluidRow(
                  column(7, style = list("padding-right: 0px;"),
                         textInput(paste0("betaCat_", i), label = (paste0("Group/Level ", i, ": ", beta_response$cat[i])), value = beta_response$cat[i], width = '80%')),
                  column(5, style = list("padding-left: 0px;"),
                         selectInput(paste0("betaCat_ord", i), label = "Level:", choices = 1:beta_nom_num_cat, selected = i, width = '70%')))
              }))
          }) 
        }
      })
      
    }
    else{ # input$beta_prim == ""
      shinyjs::hide("beta_primvars_ref")
      shinyjs::hide("beta_primvars_rename")
    }
  })
  
  output$beta_covariate <- renderUI({
    tagList(
      p(" ", style = "margin-top: 25px;"),
      h4(strong("Covariate(s)", style = "color:black")),
      p("Potential confounders (e.g., age, gender) to be adjusted for.", style = "font-size:10pt"),
      p(" ", style = "margin-bottom: +15px;"),
      prettyRadioButtons("beta_covariate_yn",label = NULL, icon = icon("check"),
                         animation = "jelly", c("None", "Covariate(s)"), selected = "None", width = '70%'),
      
      shinyjs::hidden(
        shiny::div(id = "beta_covariate_list", style = "margin-left: 2%",
                   prettyCheckboxGroup("beta_covariate_options"," Please select covariate(s)",
                                       choices = get.cov.col(chooseData$sam.dat)[!get.cov.col(chooseData$sam.dat) %in% c(input$beta_prim)], width = '70%')))
    )
  })
  
  observeEvent(input$beta_covariate_yn, {
    observeEvent(input$beta_var_type, {
      if(input$beta_var_type == "Nominal"){
        # Nominal / No covariate(s)
        if(input$beta_covariate_yn == "None"){
          output$beta_method <- renderUI({
            tagList(
              selectInput("beta_method", label = h4(strong("Method", style = "color:black")), 
                          choices = c("PERMANOVA", "MiRKAT-MC"), width = '98%')
            )
          })
        }
        # Nominal / With covariate(s)
        else{
          output$beta_method <- renderUI({
            tagList(
              selectInput("beta_method", label = h4(strong("Method", style = "color:black")), 
                          choices = c("MiRKAT-MC"), width = '98%')
            )
          })
        }
      }
      
      else if(input$beta_var_type == "Ordinal"){
        # Ordinal
        output$beta_method <- renderUI({
          tagList(
            selectInput("beta_method", label = h4(strong("Method", style = "color:black")), 
                        choices = c("MiRKAT-MC"), width = '98%')
          )
        })
      }
    })
  })
  
  observeEvent(input$beta_covariate_yn, {
    
    if(input$beta_covariate_yn == "Covariate(s)"){
      shinyjs::show("beta_covariate_list")
    }
    else{
      shinyjs::hide("beta_covariate_list")
    }
  })
  
  output$beta_run <- renderUI({
    tagList(
      actionButton("beta_runButton", (strong("Run!")), class = "btn-info")
    )
  })
  
  observe({
    beta.check.order <- c()
    for(i in 1:length(beta_response$cat)){
      beta.check.order <- c(beta.check.order, input[[paste0("betaCat_ord",i)]])
    }
    beta.check.order <- length(unique(beta.check.order))
    
    toggleState("beta_runButton", !(input$beta_prim == "") 
                & (input$beta_covariate_yn == "None" | (input$beta_covariate_yn == "Covariate(s)" & length(input$beta_covariate_options) != 0)) 
                & (((input$beta_var_type == "Ordinal") & (beta.check.order == length(beta_response$cat))) | ((input$beta_var_type == "Nominal") & (beta.check.order == length(beta_response$cat)))))
  })
  
  ## 4. TAXA ---------------------------
  
  output$taxa_primvars <- renderUI({
    tagList(
      prettyRadioButtons("dataType_taxa", label = h4(strong("Data format", style = "color:black")), animation = "jelly",
                         c("CLR (Default)", "Proportion", "Arcsine-root",  "Count (Rarefied)"), icon = icon("check"), selected = "CLR (Default)",width = '70%'),
      
      p(" ", style = "margin-top: 25px;"),
      p(" ", style = "margin-bottom: -7px;"),
      h4(strong("Response variable", style = "color:black")),
      p("E.g., a multi-categorical variable on the host's health or disease status", style = "font-size:10pt"),
      p(" ", style = "margin-bottom: +15px;"),
      selectInput("taxa_prim", label = NULL, 
                  choices = c("Choose one" = "",  bmc.col.check(chooseData$sam.dat, type = "Multinomial")), 
                  width = '70%')
    )
  })
  
  output$taxa_var_type <- renderUI({
    tagList(
      prettyRadioButtons("taxa_var_type", label = h4(strong("Variable type", style = "color:black")), animation = "jelly",
                         c("Nominal", "Ordinal"), icon = icon("check"), selected = "Nominal",width = '70%')
    )
  })
  
  observeEvent(input$taxa_prim, {
    
    if(input$taxa_prim != ""){
      
      shinyjs::show("taxa_primvars_ref")
      shinyjs::show("taxa_primvars_rename")
      
      observeEvent(input$taxa_var_type, {
        
        if(input$taxa_var_type == "Nominal"){
          
          taxa_response$cat = category.names(chooseData$sam.dat, input$taxa_prim)
          taxa_nom_num_cat = length(taxa_response$cat)
          
          output$taxa_primvars_rename <- renderUI({
            tagList(
              h4(strong("Reorder/rename categories", style = "color:black")),
              p("You can reorder/rename groups/levels of the chosen response variable. MiMultiCat keeps up to 8 characters on graphs.", style = "font-size:10pt"),
              lapply(1:taxa_nom_num_cat, function(i){
                fluidRow(
                  column(7, style = list("padding-right: 0px;"),
                         textInput(paste0("taxaCat_", i), label = (paste0("Group/Level ", i, ": ", taxa_response$cat[i])), value = taxa_response$cat[i], width = '80%')),
                  column(5, style = list("padding-left: 0px;"),
                         selectInput(paste0("taxaCat_ord", i), label = "Level:", choices = c("Reference" = 1, 2:taxa_nom_num_cat), selected = i, width = '70%')))
              }))
          }) 
          
        }
        else{
          
          taxa_response$cat = category.names(chooseData$sam.dat, input$taxa_prim)
          taxa_nom_num_cat = length(taxa_response$cat)
          
          output$taxa_primvars_rename <- renderUI({
            tagList(
              h4(strong("Reorder/rename categories", style = "color:black")),
              p("You can reorder/rename groups/levels of the chosen response variable. MiMultiCat keeps up to 8 characters on graphs.", style = "font-size:10pt"),
              lapply(1:taxa_nom_num_cat, function(i){
                fluidRow(
                  column(7, style = list("padding-right: 0px;"),
                         textInput(paste0("taxaCat_", i), label = (paste0("Group/Level ", i, ": ", taxa_response$cat[i])), value = taxa_response$cat[i], width = '80%')),
                  column(5, style = list("padding-left: 0px;"),
                         selectInput(paste0("taxaCat_ord", i), label = "Level:", choices = 1:taxa_nom_num_cat, selected = i, width = '70%')))
              }))
          }) 
        }
      })
      
    }
    else{ # input$taxa_prim == ""
      shinyjs::hide("taxa_primvars_ref")
      shinyjs::hide("taxa_primvars_rename")
    }
  })
  
  output$taxa_covariate <- renderUI({
    tagList(
      p(" ", style = "margin-top: 25px;"),
      h4(strong("Covariate(s)", style = "color:black")),
      p("Potential confounders (e.g., age, gender) to be adjusted for.", style = "font-size:10pt"),
      p(" ", style = "margin-bottom: +15px;"),
      prettyRadioButtons("taxa_covariate_yn",label = NULL, icon = icon("check"),
                         animation = "jelly", c("None", "Covariate(s)"), selected = "None", width = '70%'),
      
      shinyjs::hidden(
        shiny::div(id = "taxa_covariate_list", style = "margin-left: 2%",
                   prettyCheckboxGroup("taxa_covariate_options"," Please select covariate(s)",
                                       choices = get.cov.col(chooseData$sam.dat)[!get.cov.col(chooseData$sam.dat) %in% c(input$taxa_prim)], width = '70%')))
    )
  })
  
  observeEvent(input$taxa_covariate_yn, {
    observeEvent(input$taxa_var_type, {
      if(input$taxa_var_type == "Nominal"){
        # Nominal / No covariate(s)
        if(input$taxa_covariate_yn == "None"){
          output$taxa_method <- renderUI({
            tagList(
              selectInput("taxa_method", label = h4(strong("Method", style = "color:black")), 
                          choices = c("ANOVA", "Kruskal-Wallis", "Multinomial Logistic Regression"), width = '98%')
            )
          })
        }
        # Nominal / With covariate(s)
        else{
          output$taxa_method <- renderUI({
            tagList(
              selectInput("taxa_method", label = h4(strong("Method", style = "color:black")), 
                          choices = c("Multinomial Logistic Regression"), width = '98%')
            )
          })
        }
      }
      
      else if(input$taxa_var_type == "Ordinal"){
        # Ordinal
        output$taxa_method <- renderUI({
          tagList(
            selectInput("taxa_method", label = h4(strong("Method", style = "color:black")), 
                        choices = c("Proportional Odds Model"), width = '98%')
          )
        })
      }
    })
  })
  
  observeEvent(input$taxa_covariate_yn, {
    
    if(input$taxa_covariate_yn == "Covariate(s)"){
      shinyjs::show("taxa_covariate_list")
    }
    else{
      shinyjs::hide("taxa_covariate_list")
    }
  })
  
  output$taxa_run <- renderUI({
    tagList(
      p(" ", style = "margin-top: +25px;"),
      h4(strong("Taxonomic ranks", style = "color:black")),
      p("Taxonomic ranks to be analyzed. ‘Phylum - Genus (16S)’ can be used for 16S data; ‘Phylum - Species (Metagenomics)’ can be used for shotgun metagenomic data.", style = "font-size:10pt"),
      p(" ", style = "margin-bottom: +15px;"),
      prettyRadioButtons("include_species_taxa", label = NULL, status = "primary", animation = "jelly",
                         c("Phylum - Genus (16S)", "Phylum - Species (Metagenomics)"), selected = "Phylum - Genus (16S)",
                         icon = icon("check"), width = '80%'),
      actionButton("taxa_runButton", (strong("Run!")), class = "btn-info")
    )
  })
  
  observe({
    taxa.check.order <- c()
    for(i in 1:length(taxa_response$cat)){
      taxa.check.order <- c(taxa.check.order, input[[paste0("taxaCat_ord",i)]])
    }
    taxa.check.order <- length(unique(taxa.check.order))
    
    toggleState("taxa_runButton", !(input$taxa_prim == "") 
                & (input$taxa_covariate_yn == "None" | (input$taxa_covariate_yn == "Covariate(s)" & length(input$taxa_covariate_options) != 0)) 
                & (((input$taxa_var_type == "Ordinal") & (taxa.check.order == length(taxa_response$cat))) | ((input$taxa_var_type == "Nominal") & (taxa.check.order == length(taxa_response$cat)))))
  })
  
  ## 5-1. RF ---------------------------
  
  output$rf_nom_data_input <- renderUI({
    tagList(
      prettyRadioButtons("rf_nom_dataType", label = h4(strong("Data format", style = "color:black")), animation = "jelly",
                         c("CLR (Default)", "Count (Rarefied)", "Proportion", "Arcsine-root"), icon = icon("check"), selected = "CLR (Default)",width = '70%'),
      
      p(" ", style = "margin-top: 25px;"),
      h4(strong("Response variable", style = "color:black")),
      p("E.g., a multi-categorical variable on the host's health or disease status", style = "font-size:10pt"),
      p(" ", style = "margin-bottom: +15px;"),
      selectInput("rf_nom_response", label = NULL, 
                  choices = bmc.col.check(chooseData$sam.dat, type = "Multinomial"), 
                  width = '70%'))
  })
  
  observeEvent(input$rf_nom_response, {
    
    if(input$rf_nom_response != ""){
      rf_response$cat = category.names(chooseData$sam.dat, input$rf_nom_response)
      rf_response_length = length(rf_response$cat)
      
      output$rf_nom_data_input_opt <- renderUI({
        tagList(
          h4(strong("Reorder/rename categories", style = "color:black")),
          p("You can reorder/rename group/levels of the chosen response variable. MiMultiCat keeps up to 8 characters on graphs.", style = "font-size:10pt"),
          lapply(1:rf_response_length, function(i){
            fluidRow(
              column(7, style = list("padding-right: 0px;"),
                     textInput(paste0("rf_nom_label", i), label = (paste0("Group/Level ", i, ": ", rf_response$cat[i])), value = rf_response$cat[i], width = '80%')),
              column(5, style = list("padding-left: 0px;"),
                     selectInput(paste0("rf_nom_ord", i), label = "Level:", choices = c("Reference" = 1, 2:rf_response_length), selected = i, width = '70%')))
          }))
      })
    }
    else{
      shinyjs::hide("rf_nom_data_input_opt")
    }
  })
  
  output$rf_nom_train_setting <- renderUI({
    tagList(
      p(" ", style = "margin-top: 25px;"),
      h4(strong("# folds", style = "color:black")),
      p("The number of non-overlapping folds of the data to be used in cross-validations.", style = "font-size:10pt"),
      p(" ", style = "margin-bottom: +15px;"),
      selectInput("rf_nom_nfold", label = NULL, 
                  c("Choose one" = "", c(5, 10)), selected = 5, width = '70%'),
      
      p(" ", style = "margin-top: 25px;"),
      h4(strong("# trees", style = "color:black")),
      p("The number of bagged trees to be aggregated (Default: 1,000).", style = "font-size:10pt"),
      p(" ", style = "margin-bottom: +15px;"),
      selectInput("rf_nom_ntree", label = NULL, 
                  c("Choose one" = "", c(1000, 3000, 5000, 10000)), selected = 1000, width = '70%'),
      
      p(" ", style = "margin-top: 25px;"),
      h4(strong("# taxa to be displayed", style = "color:black")),
      p("The maximum number of taxa to be displayed in importance and partial dependence plots (Default: 20).", style = "font-size:10pt"),
      p(" ", style = "margin-bottom: +15px;"),
      sliderInput("rf_nom_var_num", label = NULL, min = 5, max = 20, value = 20, step = 5),
      
      prettyRadioButtons("rf_nom_include_species", label = h4(strong("Taxonomic ranks", style = "color:black")), animation = "jelly",
                         c("Phylum - Genus (16S)", "Phylum - Species (Metagenomics)"), selected = "Phylum - Genus (16S)",
                         icon = icon("check"), width = '80%'),
      
      actionButton("rf_nom_runButton", (strong("Run!")), class = "btn-info"))
  })
  
  observe({
    rf.check.order <- c()
    for(i in 1:length(rf_response$cat)){
      rf.check.order <- c(rf.check.order, input[[paste0("rf_nom_ord",i)]])
    }
    rf.check.order <- length(unique(rf.check.order))
    print(rf.check.order)
    toggleState("rf_nom_runButton", (rf.check.order == length(rf_response$cat))) # !(input$rf_nom_response == "") & 
  })
  
  ## 5-2. XGB ---------------------------
  
  output$xgb_nom_data_input <- renderUI({
    tagList(
      prettyRadioButtons("xgb_nom_dataType", label = h4(strong("Data format", style = "color:black")), icon = icon("check"), animation = "jelly",
                         c("CLR (Default)", "Count (Rarefied)", "Proportion", "Arcsine-root"), selected = "CLR (Default)",width = '70%'),
      
      p(" ", style = "margin-top: 25px;"),
      h4(strong("Response variable", style = "color:black")),
      p("E.g., a multi-categorical variable on the host's health or disease status", style = "font-size:10pt"),
      p(" ", style = "margin-bottom: +15px;"),
      selectInput("xgb_nom_response", label = NULL, 
                  choices = bmc.col.check(chooseData$sam.dat, type = "Multinomial"),
                  width = '70%'))
  })
  
  observeEvent(input$xgb_nom_response, {
    
    if(input$xgb_nom_response != ""){
      xgb_response$cat = category.names(chooseData$sam.dat, input$xgb_nom_response)
      xgb_response_length = length(xgb_response$cat)
      
      output$xgb_nom_data_input_opt <- renderUI({
        tagList(
          h4(strong("Reorder/rename categories", style = "color:black")),
          p("You can reorder/rename group/levels of the chosen response variable. MiMultiCat keeps up to 8 characters on graphs.", style = "font-size:10pt"),
          lapply(1:xgb_response_length, function(i){
            fluidRow(
              column(7, style = list("padding-right: 0px;"),
                     textInput(paste0("xgb_nom_label", i), label = (paste0("Group/Level ", i, ": ", xgb_response$cat[i])), value = xgb_response$cat[i], width = '80%')),
              column(5, style = list("padding-left: 0px;"),
                     selectInput(paste0("xgb_nom_ord", i), label = "Level:", choices = c("Reference" = 1, 2:xgb_response_length), selected = i, width = '70%')))
          }))
      })
    }
    else{
      shinyjs::hide("xgb_nom_data_input_opt")
    }
  })
  
  output$xgb_nom_train_setting <- renderUI({
    tagList(
      p(" ", style = "margin-top: 25px;"),
      h4(strong("# folds", style = "color:black")),
      p("The number of non-overlapping folds of the data to be used in cross-validations.", style = "font-size:10pt"),
      p(" ", style = "margin-bottom: +15px;"),
      selectInput("xgb_nom_nfold", label = NULL, 
                  c("Choose one" = "", c(5, 10)), selected = 5, width = '70%'),
      
      p(" ", style = "margin-top: 25px;"),
      h4(strong("Max # iterations", style = "color:black")),
      p("The maximum number of iterations (updates) in the boosting process. A large maximum number of iterations (e.g., 10,000) is recommended 
        to boost the tree sufficiently, but it can be at the cost of heavy computation. (Default: 1,000)", style = "font-size:10pt"),
      p(" ", style = "margin-bottom: +15px;"),
      selectInput("xgb_nom_nrounds", label = NULL, 
                  c("Choose one" = "", c(1000, 3000, 5000, 10000)), selected = 1000, width = '70%'),
      
      p(" ", style = "margin-top: 25px;"),
      h4(strong("Learning rate", style = "color:black")),
      p("The rate of a newly fitted tree to be reflected into the aggregation (update). A low learning rate (e.g., 0.001) is recommended to elaborate the boosting process through slow learning, but it is at the cost of heavy computations. (Default: 0.005)", style = "font-size:10pt"),
      p(" ", style = "margin-bottom: +15px;"),
      selectInput("xgb_nom_eta", label = NULL, 
                  c("Choose one" = "", c(0.001, 0.005, 0.01, 0.05)), selected = 0.005, width = '70%'),
      
      prettyRadioButtons("xgb_nom_penalty", label = h4(strong("Regularization", style = "color:black")), animation = "jelly",
                         c("Yes (Default)", "No"), selected = "Yes (Default)", icon = icon("check"), width = "70%"),
      
      p(" ", style = "margin-top: 25px;"),
      h4(strong("# taxa to be displayed", style = "color:black")),
      p("The maximum number of taxa per taxonomic rank to be displayed in plots.", style = "font-size:10pt"),
      p(" ", style = "margin-bottom: +15px;"),
      sliderInput("xgb_nom_var_num", label = NULL, min = 5, max = 20, value = 20, step = 5),
      
      prettyRadioButtons("xgb_nom_include_species", label = h4(strong("Taxonomic ranks", style = "color:black")), animation = "jelly",
                         c("Phylum - Genus (16S)", "Phylum - Species (Metagenomics)"), selected = "Phylum - Genus (16S)",
                         icon = icon("check"), width = '80%'),
      actionButton("xgb_nom_runButton", (strong("Run!")), class = "btn-info"))
  })
  
  observe({
    xgb.check.order <- c()
    for(i in 1:length(xgb_response$cat)){
      xgb.check.order <- c(xgb.check.order, input[[paste0("xgb_nom_ord",i)]])
    }
    xgb.check.order <- length(unique(xgb.check.order))
    print(xgb.check.order)
    toggleState("xgb_nom_runButton", (xgb.check.order == length(xgb_response$cat))) # !(input$xgb_nom_response == "") & 
  })
  
  ## Run Buttons ------
  
  ## 1-1. QC -----
  
  observeEvent(input$run, {

    withProgress(
      message = 'Calculation in progress', 
      detail = 'This may take a while...', value = 0, {
        
        incProgress(1/10, message = "Data Trimming in progress")
        
        if (nchar(input$part.rem.str) == 0) {
          rem.tax.complete <- rem.tax.d
          rem.tax.partial <- rem.tax.str.d
        } else {
          rem.tax.complete <- unique(c(unlist(strsplit(input$rem.str, split = ",")), rem.tax.d))
          rem.tax.partial <- unique(c(unlist(strsplit(input$part.rem.str, split = ",")), rem.tax.str.d))
        }
        
        tax.tab <- tax_table(infile$biom)
        
        if (input$kingdom != "all") {
          ind <- is.element(tax.tab[,1], input$kingdom)
          shiny::validate(
            if (sum(ind) == 0) {
              showNotification(h4(paste("Error: Please select valid Kingdom. Available kingdoms are:", 
                                        paste(c(na.omit(unique(tax.tab[,1])) ,"and all"), collapse = ", "))),
                               type = "error")
            } else {
              NULL
            }
          )
        }
        
        shinyjs::disable("run")
        shinyjs::disable("slider1")
        shinyjs::disable("slider2")
        shinyjs::disable("kingdom")
        shinyjs::disable("skip")
        shinyjs::disable("binwidth")
        shinyjs::disable("binwidth2")
        shinyjs::disable("rem.str")
        shinyjs::disable("part.rem.str")
        
        rcol$selected = "rgba(255, 0, 0, 0.6)"
        
        tree.exists <- !is.null(access(infile$biom, "phy_tree"))
        print(tree.exists)
        infile$qc_biom = biom.clean(infile$biom, 
                                    input$kingdom, 
                                    lib.size.cut.off = input$slider1, 
                                    mean.prop.cut.off = input$slider2/100,
                                    rem.tax = rem.tax.complete, rem.tax.str = rem.tax.partial,
                                    tree.exists = tree.exists)
        
        incProgress(3/10, message = "Rarefying in progress")
        lib_size.sum = lib.size.func(infile$qc_biom)$lib.size.sum
        infile$rare_biom = rarefy.func(infile$qc_biom, 
                                       cut.off = lib_size.sum["Minimum"],
                                       multi.rarefy = 1,
                                       tree.exists = tree.exists)
        
        incProgress(2/10, message = "Saving File in progress")
        
        chooseData$sam.dat = sample_data(infile$qc_biom)
        chooseData$mon.sin.rev.bin.con = is.mon.sin.rev.bin.con(chooseData$sam.dat)
        chooseData$prim_vars = pri.func(chooseData$sam.dat, chooseData$mon.sin.rev.bin.con)
        chooseData$tax.tab = tax_table(infile$rare_biom)
        
        output$moreControls <- renderUI({
          tagList(
            box(title = strong("Download Data", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                span(textOutput("text"), style="font-size:13pt"),
                h5("Data after Quality Control"),
                downloadButton("downloadData2", "Download", width = '50%', style = "color:black; background-color: red3"),br(),
                h5("Data after Quality Control and Rarefaction"),
                downloadButton("downloadData3", "Download", width = '50%', style = "color:black; background-color: red3"),br(),
                p("For your reference, you can download the data files above for the phyloseq object (biom.after.qc) after QC and
                    (rare.biom.after.qc) after QC and rarefaction.",
                  style = "font-size:11pt")
            )
          )
        })
        
        output$text <- renderText({"You are all set! You can proceed to data analysis!"})
        
        biom.after.qc = infile$qc_biom
        output$downloadData2 <- downloadHandler(
          filename = function() {
            paste("biom.after.qc.Rdata")
          },
          content = function(file1) {
            save(biom.after.qc, file = file1)
          })
        
        rare.biom.after.qc = infile$rare_biom
        output$downloadData3 <- downloadHandler(
          filename = function() {
            paste("rare.biom.after.qc.Rdata")
          },
          content = function(file1) {
            save(rare.biom.after.qc, file = file1)
          })
        
        incProgress(1/10, message = "Done")
        shinyjs::enable("run")
        shinyjs::enable("slider1")
        shinyjs::enable("slider2")
        shinyjs::enable("kingdom")
        shinyjs::enable("skip")
        shinyjs::enable("binwidth")
        shinyjs::enable("binwidth2")
        shinyjs::enable("rem.str")
        shinyjs::enable("part.rem.str")
      })
  })
  
  ## 1-2. DC & DT -----
  
  observeEvent(input$divdtCalcRun, {
    withProgress(
      message = 'Calculation in progress', 
      detail = 'This may take a while...', value = 0, {
        shinyjs::disable("divdtCalcRun")
        
        ### DC -----
        
        incProgress(2/10, message = "Calculating Diversity")
        # rare.otu.tab <- otu_table(infile$rare_biom)
        
        chooseData$alpha.div.rare = alpha.v1.func(infile$rare_biom)
        chooseData$alpha.div.qc = alpha.v1.func(infile$qc_biom)
        chooseData$alpha.div = chooseData$alpha.div.rare
        
        incProgress(3/10, message = "Calculating Distance")
        
        ds.Ks$res = Ds.Ks.func(infile$rare_biom, infile$qc_biom, "withTree")
        
        ## DT -----
        
        incProgress(2/10, message = "Transformation in progress")
        rare.otu.tab <- otu_table(infile$rare_biom)
        rare.tax.tab <- tax_table(infile$rare_biom)
        no.rare.otu.tab <- otu_table(infile$qc_biom)
        no.rare.tax.tab <- tax_table(infile$qc_biom)
        
        chooseData$taxa.out = tax.trans(no.rare.otu.tab, no.rare.tax.tab, rare.otu.tab, rare.tax.tab)
        
        chooseData$taxa.names.out = taxa.names.rank(chooseData$taxa.out[[1]])
        
        chooseData$tax.tab = rare.tax.tab
        
        incProgress(2/10, message = "Saving")
        
        output$divdtDownload <- renderUI({
          tagList(
            box(title = strong("Download Data", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                
                span(textOutput("text"), style="font-size:15pt"),
                strong(p("Alpha- and Beta-diversity Data",style = "font-size:12pt")),
                h5("Alpha Diversity", HTML('&emsp;'), HTML('&emsp;'), "Beta Diversity"),
                downloadButton("alphaDiv", "Download", width = '50%', style = "color:black; background-color: red3"), HTML('&emsp;'),
                downloadButton("betaDiv", "Download", width = '50%', style = "color:black; background-color: red3"), br(), br(),
                
                span(textOutput("text"), style="font-size:15pt"),
                strong(p("Taxonomic Abundance Data",style = "font-size:12pt")),
                h5("Count", HTML('&emsp;'), HTML('&emsp;'), HTML('&emsp;'), HTML('&emsp;'), HTML('&emsp;'), "Count (Rarefied)"),
                downloadButton("taxadataCount", "Download", width = '50%', style = " background-color: red3"), HTML('&emsp;'),
                downloadButton("taxadataRareCount", "Download", width = '50%', style = " background-color: red3"), br(), 
                h5("Proportion", HTML('&emsp;'), HTML('&emsp;'), HTML('&emsp;'), HTML('&nbsp;'), "CLR"),
                downloadButton("taxadataProp", "Download", width = '50%', style = " background-color: red3"), HTML('&emsp;'),
                downloadButton("taxadataCLR", "Download", width = '50%', style = " background-color: red3"), br(),
                h5("Arcsine-root"),
                downloadButton("taxadataArc", "Download", width = '50%', style = " background-color: red3"), br(), p("", style = "margin-bottom:5px")
            )
          )
        })
        
        alpha.div = chooseData$alpha.div
        
        output$alphaDiv <- downloadHandler(
          filename = function() {
            paste("Alpha.Diversity.txt")
          },
          content = function(alpha.file) {
            write.table(chooseData$alpha.div, file = alpha.file, row.names = TRUE, col.names = TRUE, sep = "\t")
          })
        
        output$betaDiv <- downloadHandler(
          filename = function() {
            paste("Beta.Diversity.zip")
          },
          content <- function(fname) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Jaccard.txt", "Bray.Curtis.txt", "U.UniFrac.txt" ,"G.UniFrac.txt", "W.UniFrac.txt")
            
            write.table(as.data.frame(ds.Ks$res$Ds$Jaccard), file = "Jaccard.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(ds.Ks$res$Ds$Bray.Curtis), file = "Bray.Curtis.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(ds.Ks$res$Ds$U.UniFrac), file = "U.UniFrac.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(ds.Ks$res$Ds$G.UniFrac), file = "G.UniFrac.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(ds.Ks$res$Ds$W.UniFrac), file = "W.UniFrac.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=fname, files=dataFiles)
          }
        )
        
        count_biom = chooseData$taxa.out$count
        rare_biom = chooseData$taxa.out$rare.count
        prop_biom = chooseData$taxa.out$prop
        clr_biom = chooseData$taxa.out$clr
        arc_biom = chooseData$taxa.out$arcsin
        
        output$taxadataCount <- downloadHandler(
          
          filename = function() {
            paste("Count.Data.zip")
          },
          content = function(count.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(count_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(count_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=count.file, files=dataFiles)
          }
        )
        
        output$taxadataRareCount <- downloadHandler(
          
          filename = function() {
            paste("Rarefied.Count.Data.zip")
          },
          content = function(rare.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(rare_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(rare_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=rare.file, files=dataFiles)
          }
        )
        output$taxadataProp <- downloadHandler(
          filename = function() {
            paste("Proportion.Data.zip")
          },
          content = function(prop.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(prop_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(prop_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=prop.file, files=dataFiles)
          }
        )
        output$taxadataCLR <- downloadHandler(
          filename = function() {
            paste("CLR.Transformed.Data.zip")
          },
          content = function(clr.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(clr_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(clr_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=clr.file, files=dataFiles)
          }
        )
        output$taxadataArc <- downloadHandler(
          filename = function() {
            paste("Arcsin.Transformed.Data.zip")
          },
          content = function(arc.file) {
            temp <- setwd(tempdir())
            on.exit(setwd(temp))
            dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
            write.table(as.data.frame(arc_biom$phylum), file = "Phylum.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$class), file = "Class.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$order), file = "Order.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$family), file = "Family.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$genus), file = "Genus.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            write.table(as.data.frame(arc_biom$species), file = "Species.txt", row.names = TRUE, col.names = TRUE, sep = "\t")
            zip(zipfile=arc.file, files=dataFiles)
          }
        )
        incProgress(1/10, message = "Done")
        shinyjs::enable("divdtCalcRun")
        
        shinyjs::show("alpha_primvars")
        shinyjs::show("alpha_primvars_ref")
        shinyjs::show("alpha_primvars_rename")
        shinyjs::show("alpha_var_type")
        shinyjs::show("alpha_covariate")
        shinyjs::show("alpha_method")
        shinyjs::show("alpha_run")
        shinyjs::show("alpha_display_results")
        shinyjs::show("alpha_display_pair")
        shinyjs::show("alpha_downloadTable")
        shinyjs::show("alpha_references")
        
        shinyjs::show("beta_primvars")
        shinyjs::show("beta_primvars_ref")
        shinyjs::show("beta_primvars_rename")
        shinyjs::show("beta_var_type")
        shinyjs::show("beta_covariate")
        shinyjs::show("beta_method")
        shinyjs::show("beta_run")
        shinyjs::show("beta_downloadTabUI")
        shinyjs::show("beta_reference")
        
        shinyjs::show("taxa_primvars")
        shinyjs::show("taxa_primvars_ref")
        shinyjs::show("taxa_primvars_rename")
        shinyjs::show("taxa_var_type")
        shinyjs::show("taxa_covariate")
        shinyjs::show("taxa_method")
        shinyjs::show("taxa_run")
        shinyjs::show("taxa_downloadTable")
        shinyjs::show("taxa_references")
        
        shinyjs::show("rf_nom_data_input")
        shinyjs::show("rf_nom_data_input_opt")
        shinyjs::show("rf_nom_covariate")
        shinyjs::show("rf_nom_train_setting")
        shinyjs::show("rf_nom_downloadTabUI")
        shinyjs::show("rf_nom_reference")
        
        shinyjs::show("xgb_nom_data_input")
        shinyjs::show("xgb_nom_data_input_opt")
        shinyjs::show("xgb_nom_covariate")
        shinyjs::show("xgb_nom_train_setting")
        shinyjs::show("xgb_nom_downloadTabUI")
        shinyjs::show("xgb_nom_reference")
      })
  })
  
  ## 2. ALPHA -------------------
  
  observeEvent(input$alpha_runButton, {
    shinyjs::disable("alpha_primvars")
    shinyjs::disable("alpha_primvars_ref")
    shinyjs::disable("alpha_primvars_rename")
    shinyjs::disable("alpha_var_type")
    shinyjs::disable("alpha_covariate")
    shinyjs::disable("alpha_method")
    shinyjs::disable("alpha_nom_downloadTabUI")
    shinyjs::disable("alpha_references")
    
    prim_length = length(names(table(chooseData$sam.dat[,input$alpha_prim])))
    
    for (i in 1:prim_length){
      shinyjs::disable(paste0("alphaCat_", i))
    }
    
    withProgress(
      message = "Calculation in progress",
      detail = "This may take a while...", value = 0, {
        
        incProgress(2/10, message = "Calculating distance")
        
        alpha.div <- chooseData$alpha.div
        sam_dat_alpha <- chooseData$sam.dat
        alpha.categors <- try(alpha.bin.cat.func(sam_dat_alpha, input$alpha_prim), silent = TRUE)
        
        alpha.bin.categos = c() 
        alpha.cat.order <- c()
        
        for (i in 1:prim_length){
          name_ind <- eval(parse(text = paste0("input$alphaCat_", i)))
          alpha.bin.categos <- c(alpha.bin.categos, name_ind)
          
          name_ind <- eval(parse(text = paste0("input$alphaCat_ord", i)))
          alpha.cat.order <- c(alpha.cat.order, name_ind)
        }
        
        sam_dat_alpha <- try(bin.cat.recode.func.mult(sam_dat_alpha, input$alpha_prim, alpha.categors, alpha.bin.categos), silent = TRUE)
        
        # Covariate
        if(input$alpha_covariate_yn == "None"){
          covariate = NULL
        } 
        else{
          covariate = input$alpha_covariate_options
        }
        
        alpha.div <- try(no.missing.alpha(alpha.div, sam_dat_alpha, input$alpha_prim, covariate), silent = TRUE)
        sam_dat_alpha <- try(no.missing.var(sam_dat_alpha, input$alpha_prim, covariate), silent = TRUE)
        
        incProgress(3/10, message = "Calculating")
        
        if(input$alpha_var_type == "Nominal"){
          names(alpha.cat.order) <- alpha.bin.categos
          
          if(input$alpha_method == "Multinomial Logistic Regression"){
            sam_dat_alpha[[input$alpha_prim]] <- factor(sam_dat_alpha[[input$alpha_prim]], levels = names(c(sort(alpha.cat.order)[-1], sort(alpha.cat.order)[1]))) # Multicategory variable
          }
          else{
            sam_dat_alpha[[input$alpha_prim]] <- factor(sam_dat_alpha[[input$alpha_prim]], levels = names(sort(alpha.cat.order))) # Multicategory variable
          }
        }
        else{
          sam_dat_alpha[[input$alpha_prim]] <- factor(sam_dat_alpha[[input$alpha_prim]]) # Multicategory variable
          names(alpha.cat.order) <- alpha.bin.categos
          sam_dat_alpha[[input$alpha_prim]] <- ordered(sam_dat_alpha[[input$alpha_prim]], levels = names(sort(alpha.cat.order)))
        }
        
        num_cat <- length(category.names(sam_dat_alpha, input$alpha_prim))
        
        # covariate as.numeric or as.character
        if(input$alpha_covariate_yn == "Covariate(s)"){
          
          for(cov in covariate){
            if(col.str.check(sam_dat_alpha, cov) == "numeric"){
              sam_dat_alpha[[cov]] <- as.numeric(sam_dat_alpha[[cov]])
            } 
            else if(col.str.check(sam_dat_alpha, cov) == "factor"){
              sam_dat_alpha[[cov]] <- as.factor(sam_dat_alpha[[cov]])
            }
          }
        }
        
        ### ANOVA ---------
        
        if (input$alpha_method == "ANOVA"){
          
          data.out <- try(data_mani(alpha.div, sam_dat_alpha, input$alpha_prim, covariate), silent = TRUE) 
          
          incProgress(3/10, message = "F test")
          
          F.test.out <- tryCatch(alpha.f.pair.mult.overall(data.out), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          tukey.out <- tryCatch(alpha.f.pair.mult.tukey(data.out), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          table.out <- list(ANOVA_F = F.test.out, Tukey_HSD = tukey.out)
          
          incProgress(3/10, message = "Graph")
          
          p.val.f.test <- try(F.test.out[,"P.value"], silent = TRUE)
          
          output$alpha_display_results = renderUI({
            tagList(
              box(title = strong("ANOVA F-test (Global Test)", style = "color:white"),  
                  align = "center", width = NULL, solidHeader = TRUE,  status = "info", 
                  plotOutput("box_plots", height = 800, width = 650)
              )
            )
          })
          
          output$box_plots = try(renderPlot({alpha.bin.f.mult(alpha.div, data.out, p.val.f.test)}), silent = TRUE)
          
          shinyjs::show("barPanel")
          shinyjs::hide("alpha_display_pair")
          
          output$barPanel = renderUI({
            navbarPage(
              title = "Tukey's HSD (Pairwise Comparison)",
              tabPanel("Observed", dataTableOutput("alpha_1")),
              tabPanel("Shannon", dataTableOutput("alpha_2")),
              tabPanel("Simpson", dataTableOutput("alpha_3")),
              tabPanel("InvSimpson", dataTableOutput("alpha_4")),
              tabPanel("Fisher", dataTableOutput("alpha_5")),
              tabPanel("Chao1", dataTableOutput("alpha_6")),
              tabPanel("ACE", dataTableOutput("alpha_7")),
              tabPanel("ICE", dataTableOutput("alpha_8")),
              tabPanel("PD", dataTableOutput("alpha_9")),
              position = "static-top"
            )
          })
          
          output$alpha_1  = try(renderDataTable({data.table(tukey.out[["Observed"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_2  = try(renderDataTable({data.table(tukey.out[["Shannon"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_3  = try(renderDataTable({data.table(tukey.out[["Simpson"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_4  = try(renderDataTable({data.table(tukey.out[["InvSimpson"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_5  = try(renderDataTable({data.table(tukey.out[["Fisher"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_6  = try(renderDataTable({data.table(tukey.out[["Chao1"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_7  = try(renderDataTable({data.table(tukey.out[["ACE"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_8  = try(renderDataTable({data.table(tukey.out[["ICE"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_9  = try(renderDataTable({data.table(tukey.out[["PD"]])}, options = list(dom = 't')), silent = TRUE)
          
          
          output$alpha_references <- renderUI({
            tagList(
              p(" ", style = "margin-top: 20px;"),
              box(title = strong("References", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  p(ALPHA_ANOVAF_REFERENCE, style = "font-size:11pt")
              )
            )
          })
          
          output$alpha_downloadTable <- renderUI({
            tagList(
              p(" ", style = "margin-top: +20px;"),
              box(title = strong("Download Output Table", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  p("You can download the data analysis outputs.",
                    style = "font-size:11pt"), 
                  h5("Data Analysis Outputs"),
                  downloadButton("downloadTabl_alpha", "Download", width = '50%', style = "color:black; background-color: red3")
              )
            )
          })
          
          output$downloadTabl_alpha <- downloadHandler(
            filename = function() {
              paste("Alpha.DA.Output.txt")
            },
            content = function(file) {
              sink(file); print(table.out); sink()    
            }
          )
          
        }
        
        ## Kruskal - Wallis ---------
        
        else if (input$alpha_method == "Kruskal-Wallis"){
          
          data.out <- try(data_mani(alpha.div, sam_dat_alpha, input$alpha_prim, covariate), silent = TRUE) 
          
          incProgress(3/10, message = "Kruskal-Wallis test")
          
          test.out <- tryCatch(alpha.kruskal.pair.mult.overall(data.out), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          dunn.out <- tryCatch(alpha.f.pair.mult.dunn(data.out), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          table.out <- list(KW = test.out, DUNN = dunn.out)
          
          incProgress(3/10, message = "Graph")
          
          p.val.test <- try(test.out[,"P.value"], silent = TRUE)
          
          output$alpha_display_results = renderUI({
            tagList(
              box(title = strong("Kruskal-Wallis test (Global Test)", style = "color:white"),  
                  align = "center", width = NULL, solidHeader = TRUE, status = "info", 
                  plotOutput("box_plots", height = 800, width = 650)
              )
            )
          })
        
          
          output$box_plots = try(renderPlot({alpha.bin.f.mult(alpha.div, data.out, p.val.test)}), silent = TRUE)
          
          shinyjs::show("barPanel")
          shinyjs::hide("alpha_display_pair")
          
          output$barPanel = renderUI({
            navbarPage(
              title = "Dunn's test (Pairwise Comparison)",
              tabPanel("Observed", dataTableOutput("alpha_1")),
              tabPanel("Shannon", dataTableOutput("alpha_2")),
              tabPanel("Simpson", dataTableOutput("alpha_3")),
              tabPanel("InvSimpson", dataTableOutput("alpha_4")),
              tabPanel("Fisher", dataTableOutput("alpha_5")),
              tabPanel("Chao1", dataTableOutput("alpha_6")),
              tabPanel("ACE", dataTableOutput("alpha_7")),
              tabPanel("ICE", dataTableOutput("alpha_8")),
              tabPanel("PD", dataTableOutput("alpha_9")),
              position = "static-top"
            )
          })
          
          output$alpha_1  = try(renderDataTable({data.table(dunn.out[["Observed"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_2  = try(renderDataTable({data.table(dunn.out[["Shannon"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_3  = try(renderDataTable({data.table(dunn.out[["Simpson"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_4  = try(renderDataTable({data.table(dunn.out[["InvSimpson"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_5  = try(renderDataTable({data.table(dunn.out[["Fisher"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_6  = try(renderDataTable({data.table(dunn.out[["Chao1"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_7  = try(renderDataTable({data.table(dunn.out[["ACE"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_8  = try(renderDataTable({data.table(dunn.out[["ICE"]])}, options = list(dom = 't')), silent = TRUE)
          output$alpha_9  = try(renderDataTable({data.table(dunn.out[["PD"]])}, options = list(dom = 't')), silent = TRUE)
          
          output$alpha_references <- renderUI({
            tagList(
              p(" ", style = "margin-top: 20px;"),
              box(title = strong("References", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  p(ALPHA_KRUSKAL_REFERENCE, style = "font-size:11pt")
              )
            )
          })
          
          output$alpha_downloadTable <- renderUI({
            tagList(
              p(" ", style = "margin-top: +20px;"),
              box(title = strong("Download Output Table", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  p("You can download the data analysis outputs.",
                    style = "font-size:11pt"), 
                  h5("Data Analysis Outputs"),
                  downloadButton("downloadTabl_alpha", "Download", width = '50%', style = "color:black; background-color: red3")
              )
            )
          })
          
          output$downloadTabl_alpha <- downloadHandler(
            filename = function() {
              paste("Alpha.DA.Output.txt")
            },
            content = function(file) {
              sink(file); print(table.out); sink()    
            }
          )
        }
        
        ## Multinomial Logistic Regression ---------
        
        else if (input$alpha_method == "Multinomial Logistic Regression"){
          
          new.alpha.div <- scale(alpha.div)
          data.out <- try(data_mani(new.alpha.div, sam_dat_alpha, input$alpha_prim, covariate), silent = TRUE) 
          
          shinyjs::hide("barPanel")
          shinyjs::show("alpha_display_pair")
          
          incProgress(3/10, message = "Multinomial Logistic Regression")
          
          # Global Test
          
          alpha.nominal.global.test.summary <- try(alpha.nominal.global.test(data.out, covariate = covariate), silent = TRUE)
          
          # Pairwise comparison
          
          alpha.nominal.logit <- try(vglm(prime_var ~ ., data = data.out, family = multinomial), silent = TRUE)
          
          alpha.ctable.list <- try(alpha.nominal.result.prep(alpha.nominal.logit, num_cat), silent = TRUE)
          
          alpha.list.len <- try(length(alpha.ctable.list), silent = TRUE) # 3 or 4 in example data
          
          output$alpha_display_results <- renderUI({
            tagList(
              box(title = strong("Multinomial Logistic Regression (Global Test)", style = "color:white"),  
                  align = "center", width = NULL, solidHeader = TRUE, status = "info", br(),
                  dataTableOutput("alpha_nominal_global", height = "auto", width = 800)
              )
            )
          })
          
          output$alpha_display_pair = renderUI({
            tagList(
              box(title = strong("Multinomial Logistic Regression (Pairwise Comparison)", style = "color:white"),  
                  align = "center", width = NULL, solidHeader = TRUE, status = "info", 
                  do.call(tabsetPanel, lapply(1:alpha.list.len, function(i){
                    tabPanel(title = paste0("Reference/", names(sort(alpha.cat.order)[-1])[i]), align = "center",
                             plotOutput(paste0("alpha_nominal_forest_plot", i), height = 500, width = 800))
                  }))
              )
            )
          })

          output$alpha_nominal_global <- try(renderDataTable({
            alpha.nominal.global.test.summary
          }), silent = TRUE)
          
          lapply(1:alpha.list.len, function(j) {
            output[[paste0("alpha_nominal_forest_plot", j)]] <- try(renderPlot({
              alpha.nominal.forest.plot(alpha.ctable.list[[j]])
            }), silent = TRUE)
          })
          
          output$alpha_references <- renderUI({
            tagList(
              p(" ", style = "margin-top: 20px;"),
              box(title = strong("References", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  p(ALPHA_MULTINOM_REFERENCE, style = "font-size:11pt")
              )
            )
          })
          
          output$alpha_downloadTable <- renderUI({
            tagList(
              p(" ", style = "margin-top: +20px;"),
              box(title = strong("Download Output Table", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  p("You can download the data analysis outputs.",
                    style = "font-size:11pt"), 
                  h5("Data Analysis Outputs"),
                  downloadButton("downloadTabl_alpha", "Download", width = '50%', style = "color:black; background-color: red3")
              )
            )
          })
          
          output$downloadTabl_alpha <- downloadHandler(
            filename = function() {
              paste("Alpha.DA.Output.txt")
            },
            content = function(file) {
              sink(file); print(alpha.ctable.list); sink()
            }
          )
          
        }
        
        ## Proportional Odds Model ---------
        
        else if (input$alpha_method == "Proportional Odds Model"){
          
          new.alpha.div <- scale(alpha.div)
          data.out <- try(data_mani(new.alpha.div, sam_dat_alpha, input$alpha_prim, covariate), silent = TRUE) 
          
          shinyjs::hide("barPanel")
          shinyjs::hide("alpha_display_pair")
          
          incProgress(3/10, message = "Proportional Odds")
          
          alpha.ordinal.logit <- try(vglm(prime_var ~ ., data = data.out, family = cumulative(parallel = TRUE)), silent = TRUE)
          
          alpha.ctable <- try(alpha.ordinal.result.prep(alpha.ordinal.logit, num_cat), silent = TRUE)
          
          output$alpha_display_results = renderUI({
            tagList(
              box(title = strong("Proportional Odds Model", style = "color:white"),  
                  align = "center", width = NULL, solidHeader = TRUE, status = "info", 
                  plotOutput("alpha_ordinal_forest_plot", height = 500, width = 800)
              )
            )
          })
          
          output$alpha_ordinal_forest_plot <- renderPlot({
            tryCatch(alpha.ordinal.forest.plot(alpha.ctable), error = function(e){
              message("No outcome is available!")
              showModal(modalDialog(div("No outcome is available!")))
              return(NULL)
            })
          })
          
          output$alpha_references <- renderUI({
            tagList(
              p(" ", style = "margin-top: 20px;"),
              box(title = strong("References", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  p(ALPHA_PROPODDS_REFERENCE, style = "font-size:11pt")
              )
            )
          })
          
          output$alpha_downloadTable <- renderUI({
            tagList(
              p(" ", style = "margin-top: +20px;"),
              box(title = strong("Download Output Table", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  p("You can download the data analysis outputs.",
                    style = "font-size:11pt"), 
                  h5("Data Analysis Outputs"),
                  downloadButton("downloadTabl_alpha", "Download", width = '50%', style = "color:black; background-color: red3")
              )
            )
          })
          
          output$downloadTabl_alpha <- downloadHandler(
            filename = function() {
              paste("Alpha.DA.Output.txt")
            },
            content = function(file) {
              sink(file); print(alpha.ctable); sink()
            }
          )
        }
      }) 
    
    
    shinyjs::enable("alpha_primvars")
    shinyjs::enable("alpha_primvars_ref")
    shinyjs::enable("alpha_primvars_rename")
    shinyjs::enable("alpha_var_type")
    shinyjs::enable("alpha_covariate")
    shinyjs::enable("alpha_method")
    shinyjs::enable("alpha_nom_downloadTabUI")
    shinyjs::enable("alpha_references")
    
  })
  
  ## 3. BETA -------------------
  
  observeEvent(input$beta_runButton, {
    shinyjs::disable("beta_primvars")
    shinyjs::disable("beta_primvars_ref")
    shinyjs::disable("beta_primvars_rename")
    shinyjs::disable("beta_var_type")
    shinyjs::disable("beta_covariate")
    shinyjs::disable("beta_method")
    shinyjs::disable("beta_downloadTabUI")
    shinyjs::disable("beta_reference")
    
    withProgress(
      message = "Calculation in progress",
      detail = "This may take a while...", value = 0, {
        
        Ds.Ks <- ds.Ks$res
        
        Ds <- Ds.Ks$Ds
        Ks <- Ds.Ks$Ks
        
        sam_dat <- chooseData$sam.dat
        
        beta.categors <- try(alpha.bin.cat.func(sam_dat, input$beta_prim), silent = TRUE)
        
        beta.bin.categos <- c() 
        beta.cat.order <- c()
        for (i in 1:length(beta.categors)){
          name_ind <- eval(parse(text = paste0("input$betaCat_", i)))
          beta.bin.categos <- c(beta.bin.categos, name_ind)
          
          name_ind <- eval(parse(text = paste0("input$betaCat_ord", i)))
          beta.cat.order <- c(beta.cat.order, name_ind)
        }
        
        sam_dat <- bin.cat.recode.func.mult(sam_dat, input$beta_prim, beta.categors, beta.bin.categos)
        
        sam_dat[[input$beta_prim]] <- factor(sam_dat[[input$beta_prim]], levels = beta.bin.categos) # Multicategory variable
        
        if(input$beta_covariate_yn == "Covariate(s)"){
          
          covariate = input$beta_covariate_options
          na_ind <- sort(unique(which(is.na(sam_dat[,c(input$beta_prim, covariate)]), arr.ind = TRUE)[,1]))
          
          for(cov in covariate){
            if(col.str.check(sam_dat, cov) == "numeric"){
              sam_dat[[cov]] <- as.numeric(sam_dat[[cov]])
            } else if(col.str.check(sam_dat, cov) == "factor"){
              sam_dat[[cov]] <- as.factor(sam_dat[[cov]])
            }
          }
        }
        else{
          covariate = NULL
          na_ind <- which(is.na(sam_dat[[input$beta_prim]]))
        }
        
        if(length(na_ind) != 0){
          sam_dat <- sam_dat[-na_ind,]
          for(i in 1:length(Ds)){
            Ds[[i]] <- Ds[[i]][-na_ind, -na_ind]
          }
          for(i in 1:length(Ks)){
            Ks[[i]] <- Ks[[i]][-na_ind, -na_ind]
          }
        }
        
        re_sam_dat <- sam_dat 
        
        ### MiRKAT-MC -----
        
        if (input$beta_method == "MiRKAT-MC"){
          
          incProgress(3/10, message = "MiRKAT-MC") 
          
          var.type <- ifelse(input$beta_var_type == "Nominal", "nominal", "ordinal")
          
          if(var.type  == "Nominal"){
            names(beta.cat.order) <- beta.bin.categos
            y <- factor(sam_dat[[input$beta_prim]], levels = names(sort(beta.cat.order))) # Multicategory variable
          }
          else{
            names(beta.cat.order) <- beta.bin.categos
            y <- factor(sam_dat[[input$beta_prim]]) # Multicategory variable
            y <- ordered(sam_dat[[input$beta_prim]], levels = names(sort(beta.cat.order)))
          }
          
          # y <- as.factor(sam_dat[[input$beta_prim]]) # Multicategory variable
          
          if(is.null(covariate)){
            dat <- data.frame(y)
          } 
          else{
            dat <- cbind(y, sam_dat[,covariate])
          }
          
          incProgress(2/10, message = "Calculating p-values")
          mirkatmc.pvs.nom <- numeric()
          
          set.seed(521)
          for (i in 1:length(Ks)) {
            mirkatmc.pvs.nom[i] <- try(MiRKATMC(formula = y ~ ., random = NULL, data.type = 'nominal', Ks = Ks[i], data = dat), silent = TRUE)
          }
          
          output$beta_nom_results <- renderUI({
            tagList(
              box(width = NULL, status = "info", solidHeader = TRUE, align = "center",
                  title = strong("PCoA Plot for Beta Diversity", style = "color:white"),
                  plotOutput("beta_nom_PCoA", height = 800, width = 550))
            )
          })
          
          output$beta_nom_PCoA <- renderPlot({
            tryCatch(beta.PCoA.plot(Ds = Ds, pvs = mirkatmc.pvs.nom, y = y), error = function(e){
              message("No outcome is available!")
              showModal(modalDialog(div("No outcome is available!")))
              return(NULL)})
            
          })
          shinyjs::hide("beta_barPanel")
          
          output$beta_downloadTabUI <- renderUI({
            tagList(
              p(" ", style = "margin-top: 20px;"),
              box(title = strong("Download Output Table", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE, 
                  p("You can download the data analysis outputs.",
                    style = "font-size:11pt"),
                  h5("Data Analysis Output"),
                  downloadButton("beta_nom_downloadTable", "Download", width = '50%', style = "background-color: red3")
              )
            )
          })
          
          output$beta_nom_downloadTable <- downloadHandler(
            filename = function() {
              paste("Beta.DA.Output.txt")
            },
            content = function(file) {
              out_temp = data.frame(mirkatmc.pvs.nom)
              rownames(out_temp) = c("Jaccard","Bray.Curtis","U.UniFrac","G.UniFrac","W.UniFrac")
              colnames(out_temp) = "p-value"
              write.table(out_temp, file, sep="\t")
            }
          )
          
          output$beta_reference <- renderUI({
            tagList(
              p(" ", style = "margin-top: 20px;"),
              box(title = strong("References", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  p(BETA_DA_REFERENCE, style = "font-size:11pt")
              )
            )
          })
        }
        
        ## PERMANOVA -----
        
        else if (input$beta_method == "PERMANOVA") {
          incProgress(3/10, message = "PERMANOVA") 
          
          names(beta.cat.order) <- beta.bin.categos
          sam_dat[[input$beta_prim]] <- factor(sam_dat[[input$beta_prim]], levels = names(sort(beta.cat.order))) # Multicategory variable
          
          set.seed(521)
          beta_div <- try(reduced_beta(Ds, sam_dat, input$beta_prim), silent = TRUE)   #요기 
          
          dat_1 <- beta_div 
          dat_2 <- sam_dat 
          dat_3 <- input$beta_prim
          
          # permanova.out <- tryCatch(permanova.mult.united(num_perm = 3000, beta_div, sam_dat, input$beta_prim, method_adj = "BH"), error = function(e) {
          #   message ("No outcome is available!")
          #   showModal(modalDialog(div("No outcome is available!")))
          #   return(NULL)
          # })
          # 
          # permanova.pair.out <- permanova.out$table
          # permanova.pair.out.2 <- permanova.out$download 
          
          set.seed(521)
          permanova.global <- tryCatch(beta.permanova.mult.glob(num_perm = 3000, beta_div, sam_dat, input$beta_prim, download = TRUE), 
                                        error = function(e) {
                                          message ("No outcome is available!")
                                          showModal(modalDialog(div("No outcome is available!")))
                                          return(NULL)})
          
          permanova.global.out <- permanova.global$table 
          table.out <- permanova.global.out
          # table.out <- try(list(Global_PERMANOVA = permanova.global.out, Pairwise_PERMANOVA = permanova.pair.out.2))
          
          # shinyjs::show("beta_barPanel")
          # output$beta_barPanel = renderUI({
          #   navbarPage(
          #     title = "PERMANOVA (Pairwise Comparison)",
          #     tabPanel("Jaccard", dataTableOutput("beta_1")),
          #     tabPanel("Bray.Curtis", dataTableOutput("beta_2")), 
          #     tabPanel("U.UniFrac", dataTableOutput("beta_3")), 
          #     tabPanel("G.UniFrac", dataTableOutput("beta_4")), 
          #     tabPanel("W.UniFrac", dataTableOutput("beta_5")),
          #     position = "static-top"
          #   )
          # })
          # 
          # output$beta_1  = try(renderDataTable({data.table(permanova.pair.out[["Jaccard"]])}, options = list(dom = 't')), silent = TRUE)
          # output$beta_2  = try(renderDataTable({data.table(permanova.pair.out[["Bray.Curtis"]])}, options = list(dom = 't')), silent = TRUE)
          # output$beta_3  = try(renderDataTable({data.table(permanova.pair.out[["U.UniFrac"]])}, options = list(dom = 't')), silent = TRUE)
          # output$beta_4  = try(renderDataTable({data.table(permanova.pair.out[["G.UniFrac"]])}, options = list(dom = 't')), silent = TRUE)
          # output$beta_5  = try(renderDataTable({data.table(permanova.pair.out[["W.UniFrac"]])}, options = list(dom = 't')), silent = TRUE)
          
          output$beta_nom_results <- renderUI({
            tagList(
              box(width = NULL, status = "info", solidHeader = TRUE, align = "center",
                  title = strong("PCoA Plot for Beta Diversity", style = "color:white"),
                  plotOutput("beta_nom_PCoA", height = 800, width = 550))
            )
          })
          
          y <- as.factor(sam_dat[[input$beta_prim]])
          
          output$beta_nom_PCoA <- renderPlot({
            beta.PCoA.plot(Ds = Ds, 
                           pvs = permanova.global.out[, "P.value"],
                           y = y)
          })
          
          
          output$beta_reference <- renderUI({
            tagList(
              p(" ", style = "margin-top: 20px;"),
              box(title = strong("References", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  p(BETA_PERMANOVA_REFERENCE, style = "font-size:11pt")
              )
            )
          })
        }
        
        output$beta_downloadTabUI <- renderUI({
          tagList(
            p(" ", style = "margin-top: +20px;"),
            box(title = strong("Download Output Table", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                p("You can download the data analysis outputs.",
                  style = "font-size:11pt"), 
                h5("Data Analysis Outputs"),
                downloadButton("downloadTabl_beta", "Download", width = '50%', style = "color:black; background-color: red3")
            )
          )
        })
        
        output$downloadTabl_beta <- downloadHandler(
          filename = function() {
            paste("Beta.DA.Output.txt")
          },
          content = function(file) {
            sink(file); print(table.out); sink()    
          }
        )
        
        incProgress(6/10, message = "SAVE")
      }
    )
    shinyjs::enable("beta_primvars")
    shinyjs::enable("beta_primvars_ref")
    shinyjs::enable("beta_primvars_rename")
    shinyjs::enable("beta_var_type")
    shinyjs::enable("beta_covariate")
    shinyjs::enable("beta_method")
    shinyjs::enable("beta_downloadTabUI")
    shinyjs::enable("beta_reference")
  })
  
  ## 4. TAXA -------------------
  
  observeEvent(input$taxa_runButton, {
    shinyjs::disable("taxa_primvars")
    shinyjs::disable("taxa_primvars_ref")
    shinyjs::disable("taxa_primvars_rename")
    shinyjs::disable("taxa_var_type")
    shinyjs::disable("taxa_covariate")
    shinyjs::disable("taxa_method")
    shinyjs::disable("include_species_taxa")
    shinyjs::disable("taxa_downloadTabUI")
    shinyjs::disable("taxa_references")
    
    prim_length = length(names(table(chooseData$sam.dat[,input$taxa_prim])))
    
    for (i in 1:prim_length){
      shinyjs::disable(paste0("taxaCat_", i))
    } 
    
    withProgress(
      message = "Calculation in progress",
      detail = "This may take a while...", value = 0, {
        
        incProgress(2/10, message = "Calculating distance")
        
        sam_dat_taxa <- chooseData$sam.dat
        taxa.out <- chooseData$taxa.out
        
        taxa.names.rank.out <- taxa.names.rank(taxa.out$count)
        
        for(j in 1:6){
          for(k in 1:6){
            names(taxa.out[[j]][[k]]) <- taxa.names.rank.out$names[[k]]
          }
        }
        
        if(input$taxa_covariate_yn == "Covariate(s)"){
          covariate <- input$taxa_covariate_options
        } 
        else{
          covariate <- NULL
        } 
        
        if(input$include_species_taxa == "Phylum - Genus (16S)"){
          include <- FALSE
        }
        else{
          include <- TRUE
        }
        level.names <- get.level.names(include = include)
        
        if (input$dataType_taxa == "Count (Rarefied)") {
          taxa.dataType = "rare.count"
        } else if (input$dataType_taxa == "CLR (Default)") {
          taxa.dataType = "clr"
        } else if (input$dataType_taxa == "Arcsine-root") {
          taxa.dataType = "arcsin"  
        } else if (input$dataType_taxa == "Proportion"){
          taxa.dataType = "imp.prop"
        }
        
        # Reference level
        
        taxa.categors <- try(alpha.bin.cat.func(sam_dat_taxa, input$taxa_prim), silent = TRUE)
        
        taxa.bin.categos <- c() 
        taxa.cat.order <- c()
        
        for (i in 1:prim_length){
          name_ind <- eval(parse(text = paste0("input$taxaCat_", i)))
          taxa.bin.categos <- c(taxa.bin.categos, name_ind)
          
          name_ind <- eval(parse(text = paste0("input$taxaCat_ord", i)))
          taxa.cat.order <- c(taxa.cat.order, name_ind)
        }
        
        sam_dat_taxa <- try(bin.cat.recode.func.mult(sam_dat_taxa, input$taxa_prim, taxa.categors, taxa.bin.categos), silent = TRUE)
        
        taxa.mani <- try(taxa.data.mani(sam_dat_taxa, taxa.out, input$taxa_prim, taxa.dataType, covariate), silent = TRUE) 
        sam_dat_taxa <- try(no.missing.var(sam_dat_taxa, input$taxa_prim), silent = TRUE)
        
        incProgress(3/10, message = "Calculating")
        
        sam_dat_taxa[[input$taxa_prim]] <- factor(sam_dat_taxa[[input$taxa_prim]], levels = taxa.bin.categos)
        num_cat <- length(category.names(sam_dat_taxa, input$taxa_prim))
        
        # covariate as.numeric or as.character
        for(name in level.names){
          
          if(input$taxa_var_type == "Nominal"){
            names(taxa.cat.order) <- taxa.bin.categos
            
            if(input$taxa_method == "Multinomial Logistic Regression"){
              taxa.mani[[name]][["prim.var"]] <- factor(taxa.mani[[name]][["prim.var"]], levels = names(c(sort(taxa.cat.order)[-1], sort(taxa.cat.order)[1]))) # Multicategory variable
            }
            else{
              taxa.mani[[name]][["prim.var"]] <- factor(taxa.mani[[name]][["prim.var"]], levels = names(sort(taxa.cat.order))) # Multicategory variable
            }
          }
          else{
            names(taxa.cat.order) <- taxa.bin.categos
            taxa.mani[[name]][["prim.var"]] <- factor(taxa.mani[[name]][["prim.var"]]) # Multicategory variable
            taxa.mani[[name]][["prim.var"]] <- ordered(taxa.mani[[name]][["prim.var"]], levels = names(sort(taxa.cat.order)))
          }
          
          if(input$taxa_covariate_yn == "Covariate(s)"){
            for(cov in covariate){
              if(col.str.check(taxa.mani[[name]], cov) == "numeric"){
                taxa.mani[[name]][[cov]] <- as.numeric(taxa.mani[[name]][[cov]])
              }
              else if(col.str.check(taxa.mani[[name]], cov) == "factor"){
                taxa.mani[[name]][[cov]] <- as.factor(taxa.mani[[name]][[cov]])
              }
            }
          }
        }

        
        ### ANOVA ---------
        
        if (input$taxa_method == "ANOVA"){
          
          shinyjs::show("taxa_display_area")
          shinyjs::show("taxa_barPanel")
          shinyjs::hide("taxa_display_global")
          shinyjs::show("taxa_display")
          shinyjs::hide("taxa_display_forest")
          
          incProgress(3/10, message = "F test")
          
          f.overall.result <<- tryCatch(taxa.f.pair.mult.overall(taxa.mani, include), 
                                        error = function(e) {
                                          message ("No outcome is available!")
                                          showModal(modalDialog(div("No outcome is available!")))
                                          return(NULL)})
          
          taxa.rmanova.out.ori <<- try(f.overall.result$global_test, silent = TRUE)
          global.p.val.only <<- try(f.overall.result$pval, silent = TRUE)
          
          taxa.rmanova.out <<- tryCatch(q_val_combined_table.anova(taxa.rmanova.out.ori, global.p.val.only), 
                                        error = function(e) {
                                          message ("No outcome is available!")
                                          showModal(modalDialog(div("No outcome is available!")))
                                          return(NULL)})
          
          
          pairwise_result.ori <<- tryCatch(taxa.f.pair.mult.tukey(taxa.mani, include), 
                                           error = function(e) {
                                             message ("No outcome is available!")
                                             showModal(modalDialog(div("No outcome is available!")))
                                             return(NULL)})
          
          tukey.p.list <<- try(make_p_val_list(pairwise_result.ori), silent = TRUE)
          tukey.q.list <<- try(make_q_val_list(tukey.p.list), silent = TRUE)
          pairwise_result <<- try(q_val_dat(tukey.q.list, pairwise_result.ori, FALSE), silent = TRUE)
          
          num_row = numeric()
          
          for (r in 1:(5+include)) {
            row.num <- ceiling(sum(as.numeric(global.p.val.only[[r]])[complete.cases(as.numeric(global.p.val.only[[r]]))] < 0.05)/4)  
            if (row.num > 0) {
              num_row[r] <- row.num
            } else {
              num_row[r] <- 1
            } 
          }
          
          if (prim_length < 4){
            width_boxplot <- 750 
          }else if (4 <= prim_length | prim_length <= 6){
            width_boxplot <- 800 
          }else {
            width_boxplot <- 850 
          }
          
          output$taxa_references <- renderUI({
            tagList(
              p(" ", style = "margin-top: 20px;"),
              box(title = strong("References", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  p(ALPHA_ANOVAF_REFERENCE, style = "font-size:11pt")
              )
            )
          })
          
          
          if (input$include_species_taxa == "Phylum - Genus (16S)"){
            
            incProgress(3/10, message = "Displaying Results in progress")
            
            output$taxa_display = renderUI({
              tagList(
                tabBox(title = strong("ANOVA F-test (Global Test)", style = "color:black", side = "right"), width = NULL,
                       tabPanel("Phylum", align = "center",
                                plotOutput("rank1", height = num_row[1]*250, width = width_boxplot),
                       )
                       ,
                       tabPanel("Class", align = "center",
                                plotOutput("rank2", height = num_row[2]*250, width = width_boxplot),
                       )
                       ,tabPanel("Order", align = "center",
                                 plotOutput("rank3", height = num_row[3]*250, width = width_boxplot),
                       )
                       ,tabPanel("Family", align = "center",
                                 plotOutput("rank4", height = num_row[4]*250, width = width_boxplot),
                       )
                       ,tabPanel("Genus", align = "center",
                                 plotOutput("rank5", height = num_row[5]*250, width = width_boxplot),
                       )
                )
              )
            })
            
            output$rank1 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 1, global.p.val.only)
            }), silent = TRUE)
            
            output$rank2 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 2, global.p.val.only)
            }), silent = TRUE)
            
            output$rank3 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 3, global.p.val.only)
            }), silent = TRUE)
            
            output$rank4 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 4, global.p.val.only)
            }), silent = TRUE)
            
            output$rank5 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 5, global.p.val.only)
            }), silent = TRUE)
            
            sig_global_taxa <- try(sig.taxa(global.p.val.only), silent = TRUE)
            sig_pair <- try(rank_sig(pairwise_result, sig_global_taxa), silent = TRUE)
            sig_pair_dat <- try(list_to_dat(sig_pair), silent = TRUE)
            
            output$taxa_barPanel = renderUI({ 
              navbarPage(
                title = "Tukey's HSD (Pairwise Comparison)",
                tabPanel("Phylum", dataTableOutput("phylum_list")),
                tabPanel("Class", dataTableOutput("class_list")),
                tabPanel("Order", dataTableOutput("order_list")),
                tabPanel("Family", dataTableOutput("family_list")),
                tabPanel("Genus", dataTableOutput("genus_list")),
                position = "static-top"
              )
            })
            
            length_each_page <- nrow(as.matrix(sig_pair_dat[[1]]))/length(sig_pair[[1]])
            
            length_page_1 <- nrow(as.matrix(sig_pair_dat[[1]]))/length(sig_pair[[1]])
            length_page_2 <-  nrow(as.matrix(sig_pair_dat[[2]]))/length(sig_pair[[2]])
            length_page_3 <- nrow(as.matrix(sig_pair_dat[[3]]))/length(sig_pair[[3]])
            length_page_4 <- nrow(as.matrix(sig_pair_dat[[4]]))/length(sig_pair[[4]])
            length_page_5 <- nrow(as.matrix(sig_pair_dat[[5]]))/length(sig_pair[[5]])
            
            for (i in 1:length(sig_pair_dat)){
              colnames(sig_pair_dat[[i]])[colnames(sig_pair_dat[[i]]) == "Adj.P.value"] = "Adj. P.value"
              if("P.value" %in% colnames(sig_pair_dat[[i]])){
                ind <- which(colnames(sig_pair_dat[[i]]) == "P.value")
                sig_pair_dat[[i]] <- sig_pair_dat[[i]][, -c(ind)]
              }
            }
            
            output$phylum_list <- try(renderDataTable({datatable(sig_pair_dat[[1]], rownames = FALSE, options = list(pageLength = length_page_1, dom = '<"top" p>'))}), silent = TRUE)
            output$class_list <- try(renderDataTable({datatable(sig_pair_dat[[2]], rownames = FALSE, options = list(pageLength =length_page_2, dom = '<"top" p>'))}), silent = TRUE)
            output$order_list <- try(renderDataTable({datatable(sig_pair_dat[[3]], rownames = FALSE, options = list(pageLength = length_page_3, dom = '<"top" p>'))}), silent = TRUE)
            output$family_list <- try(renderDataTable({datatable(sig_pair_dat[[4]], rownames = FALSE, options = list(pageLength = length_page_4, dom = '<"top" p>'))}), silent = TRUE)
            output$genus_list <- try(renderDataTable({datatable(sig_pair_dat[[5]], rownames = FALSE, options = list(pageLength = length_page_5, dom = '<"top" p>'))}), silent = TRUE)
          
            output$taxa_downloadTable = renderUI({
              tagList(
                p(" ", style = "margin-bottom: +20px;"),
                box(title = strong("Download Output Table", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                    p("You can download the data analysis outputs.",
                      style = "font-size:11pt"),
                    h5("Data Analysis Outputs"),
                    downloadButton("tdownloadTabl", "Download", width = '50%', style = "color:black; background-color: red3")
                )
              )
            })
            
            output$tdownloadTabl <- downloadHandler(
              filename = function() {
                paste("Taxa.Analysis.Output.zip")
              }, 
              content = function(DA.file) {
                temp <- setwd(tempdir())
                on.exit(setwd(temp))
                dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt")
                
                capture.output(list(Global = as.data.frame(taxa.rmanova.out[[1]]), Pairwise = pairwise_result[[1]]), file = "Phylum.txt")
                capture.output(list(Global = as.data.frame(taxa.rmanova.out[[2]]), Pairwise = pairwise_result[[2]]), file = "Class.txt")
                capture.output(list(Global = as.data.frame(taxa.rmanova.out[[3]]), Pairwise = pairwise_result[[3]]), file = "Order.txt")
                capture.output(list(Global = as.data.frame(taxa.rmanova.out[[4]]), Pairwise = pairwise_result[[4]]), file = "Family.txt")
                capture.output(list(Global = as.data.frame(taxa.rmanova.out[[5]]), Pairwise = pairwise_result[[5]]), file = "Genus.txt")
                zip(zipfile=DA.file, files=dataFiles)
              }
            )
            
            }else{
            incProgress(3/10, message = "Displaying Results in progress")
            
            output$taxa_display = renderUI({
              tagList(
                tabBox(title = strong("ANOVA F-test (Global Test)", style = "color:black", side = "right"), width = NULL,
                       tabPanel("Phylum", align = "center",
                                plotOutput("rank1", height = num_row[1]*250, width = width_boxplot),
                       )
                       ,
                       tabPanel("Class", align = "center",
                                plotOutput("rank2", height = num_row[2]*250, width = width_boxplot),
                       )
                       ,tabPanel("Order", align = "center",
                                 plotOutput("rank3", height = num_row[3]*250, width = width_boxplot),
                       )
                       ,tabPanel("Family", align = "center",
                                 plotOutput("rank4", height = num_row[4]*250, width = width_boxplot),
                       )
                       ,tabPanel("Genus", align = "center",
                                 plotOutput("rank5", height = num_row[5]*250, width = width_boxplot),
                       )
                       ,tabPanel("Species", align = "center",
                                 plotOutput("rank6", height = num_row[6]*250, width = width_boxplot),
                       )
                )
              )
            })
            
            output$rank1 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 1, global.p.val.only)
            }), silent = TRUE)
            
            output$rank2 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 2, global.p.val.only)
            }), silent = TRUE)
            
            output$rank3 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 3, global.p.val.only)
            }), silent = TRUE)
            
            output$rank4 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 4, global.p.val.only)
            }), silent = TRUE)
            
            output$rank5 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 5, global.p.val.only)
            }), silent = TRUE)
            
            output$rank6 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 6, global.p.val.only)
            }), silent = TRUE)
            
            sig_global_taxa <<- try(sig.taxa(global.p.val.only), silent = TRUE)
            sig_pair <<- try(rank_sig(pairwise_result, sig_global_taxa), silent = TRUE)
            sig_pair_dat <<- try(list_to_dat(sig_pair), silent = TRUE)
            
            output$taxa_barPanel = renderUI({ 
              navbarPage(
                title = "Tukey's HSD (Pairwise Comparison)",
                tabPanel("Phylum", dataTableOutput("phylum_list")),
                tabPanel("Class", dataTableOutput("class_list")),
                tabPanel("Order", dataTableOutput("order_list")),
                tabPanel("Family", dataTableOutput("family_list")),
                tabPanel("Genus", dataTableOutput("genus_list")),
                tabPanel("Species", dataTableOutput("species_list")),
                position = "static-top"
              )
            })
            
            length_each_page <- nrow(as.matrix(sig_pair_dat[[1]]))/length(sig_pair[[1]])
            
            length_page_1 <- nrow(as.matrix(sig_pair_dat[[1]]))/length(sig_pair[[1]])
            length_page_2 <-  nrow(as.matrix(sig_pair_dat[[2]]))/length(sig_pair[[2]])
            length_page_3 <- nrow(as.matrix(sig_pair_dat[[3]]))/length(sig_pair[[3]])
            length_page_4 <- nrow(as.matrix(sig_pair_dat[[4]]))/length(sig_pair[[4]])
            length_page_5 <- nrow(as.matrix(sig_pair_dat[[5]]))/length(sig_pair[[5]])
            length_page_6 <- nrow(as.matrix(sig_pair_dat[[6]]))/length(sig_pair[[6]])
            
            for (i in 1:length(sig_pair_dat)){
              colnames(sig_pair_dat[[i]])[colnames(sig_pair_dat[[i]]) == "Adj.P.value"] = "Adj. P.value"
              if("P.value" %in% colnames(sig_pair_dat[[i]])){
                ind <- which(colnames(sig_pair_dat[[i]]) == "P.value")
                sig_pair_dat[[i]] <- sig_pair_dat[[i]][, -c(ind)]
              }
            }
            
            output$phylum_list <- try(renderDataTable({datatable(sig_pair_dat[[1]], rownames = FALSE, options = list(pageLength = length_page_1, dom = '<"top" p>'))}), silent = TRUE)
            output$class_list <- try(renderDataTable({datatable(sig_pair_dat[[2]], rownames = FALSE, options = list(pageLength =length_page_2, dom = '<"top" p>'))}), silent = TRUE)
            output$order_list <- try(renderDataTable({datatable(sig_pair_dat[[3]], rownames = FALSE, options = list(pageLength = length_page_3, dom = '<"top" p>'))}), silent = TRUE)
            output$family_list <- try(renderDataTable({datatable(sig_pair_dat[[4]], rownames = FALSE, options = list(pageLength = length_page_4, dom = '<"top" p>'))}), silent = TRUE)
            output$genus_list <- try(renderDataTable({datatable(sig_pair_dat[[5]], rownames = FALSE, options = list(pageLength = length_page_5, dom = '<"top" p>'))}), silent = TRUE)
            output$species_list <- try(renderDataTable({datatable(sig_pair_dat[[6]], rownames = FALSE, options = list(pageLength = length_page_6, dom = '<"top" p>'))}), silent = TRUE)
          
            output$taxa_downloadTable = renderUI({
              tagList(
                p(" ", style = "margin-bottom: +20px;"),
                box(title = strong("Download Output Table", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                    p("You can download the data analysis outputs.",
                      style = "font-size:11pt"),
                    h5("Data Analysis Outputs"),
                    downloadButton("tdownloadTabl", "Download", width = '50%', style = "color:black; background-color: red3")
                )
              )
            })
            
            output$tdownloadTabl <- downloadHandler(
              filename = function() {
                paste("Taxa.Analysis.Output.zip")
              }, 
              content = function(DA.file) {
                temp <- setwd(tempdir())
                on.exit(setwd(temp))
                dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
                
                capture.output(list(Global = as.data.frame(taxa.rmanova.out[[1]]), Pairwise = pairwise_result[[1]]), file = "Phylum.txt")
                capture.output(list(Global = as.data.frame(taxa.rmanova.out[[2]]), Pairwise = pairwise_result[[2]]), file = "Class.txt")
                capture.output(list(Global = as.data.frame(taxa.rmanova.out[[3]]), Pairwise = pairwise_result[[3]]), file = "Order.txt")
                capture.output(list(Global = as.data.frame(taxa.rmanova.out[[4]]), Pairwise = pairwise_result[[4]]), file = "Family.txt")
                capture.output(list(Global = as.data.frame(taxa.rmanova.out[[5]]), Pairwise = pairwise_result[[5]]), file = "Genus.txt")
                capture.output(list(Global = as.data.frame(taxa.rmanova.out[[6]]), Pairwise = pairwise_result[[6]]), file = "Species.txt")
                zip(zipfile=DA.file, files=dataFiles)
              }
            )
            
          }
          
        }
        
        ## Kruskal-Wallis ---------
        
        else if (input$taxa_method == "Kruskal-Wallis"){
          
          shinyjs::show("taxa_display_area")
          shinyjs::show("taxa_barPanel")
          shinyjs::hide("taxa_display_global")
          shinyjs::show("taxa_display")
          shinyjs::hide("taxa_display_forest")
          
          incProgress(3/10, message = "Kruskal-Wallis test")
          
          test.out <<- tryCatch(taxa.kruskal.mult.overall(taxa.mani, include), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          
          taxa.kruskal.out.ori <<- try(test.out $global_test, silent = TRUE)
          global.p.val.only <<- try(test.out$pval, silent = TRUE)
          
          taxa.kruskal.out <<- tryCatch(q_val_combined_table.kruskal(taxa.rmanova.out.ori, global.p.val.only), 
                                        error = function(e) {
                                          message ("No outcome is available!")
                                          showModal(modalDialog(div("No outcome is available!")))
                                          return(NULL)})
          
          
          pairwise_result.ori <<- tryCatch(taxa.pair.mult.dunn(taxa.mani, include), error = function(e) {
            message ("No outcome is available!")
            showModal(modalDialog(div("No outcome is available!")))
            return(NULL)
          })
          
          dunn.p.list <<- try(make_p_val_list(pairwise_result.ori), silent = TRUE)
          dunn.q.list <<- try(make_q_val_list(dunn.p.list), silent = TRUE)
          pairwise_result <<- try(q_val_dat(dunn.q.list, pairwise_result.ori, FALSE), silent = TRUE)
          
          num_row = numeric()
          
          for (r in 1:(5+include)) {
            row.num <- ceiling(sum(as.numeric(global.p.val.only[[r]])[complete.cases(as.numeric(global.p.val.only[[r]]))] < 0.05)/4)  
            if (row.num > 0) {
              num_row[r] <- row.num
            } else {
              num_row[r] <- 1
            } 
          }
          
          num_row <<- num_row  
          if (prim_length < 4){
            width_boxplot <- 750 
          }else if (4 <= prim_length | prim_length <= 6){
            width_boxplot <- 800 
          }else {
            width_boxplot <- 850 
          }
          
          
          output$taxa_references <- renderUI({
            tagList(
              p(" ", style = "margin-top: 20px;"),
              box(title = strong("References", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  p(ALPHA_KRUSKAL_REFERENCE, style = "font-size:11pt")
              ) 
            )
          })
          
          
          if (input$include_species_taxa == "Phylum - Genus (16S)"){
            
            incProgress(3/10, message = "Displaying Results in progress")
            
            
            output$taxa_display = renderUI({
              tagList(
                
                tabBox(title = strong("Kruskal-Wallis test (Global Test)", style = "color:black", side = "right"), width = NULL,
                       tabPanel("Phylum", align = "center",
                                plotOutput("rank1", height = num_row[1]*250, width = width_boxplot),
                       )
                       ,
                       tabPanel("Class", align = "center",
                                plotOutput("rank2", height = num_row[2]*250, width = width_boxplot),
                       )
                       ,tabPanel("Order", align = "center",
                                 plotOutput("rank3", height = num_row[3]*250, width = width_boxplot),
                       )
                       ,tabPanel("Family", align = "center",
                                 plotOutput("rank4", height = num_row[4]*250, width = width_boxplot),
                       )
                       ,tabPanel("Genus", align = "center",
                                 plotOutput("rank5", height = num_row[5]*250, width = width_boxplot),
                       )
                )
              )
            })
            
            output$rank1 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 1, global.p.val.only)
            }), silent = TRUE)
            
            output$rank2 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 2, global.p.val.only)
            }), silent = TRUE)
            
            output$rank3 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 3, global.p.val.only)
            }), silent = TRUE)
            
            output$rank4 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 4, global.p.val.only)
            }), silent = TRUE)
            
            output$rank5 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 5, global.p.val.only)
            }), silent = TRUE)
            
            sig_global_taxa <- try(sig.taxa(global.p.val.only), silent = TRUE)
            sig_pair <- try(rank_sig(pairwise_result, sig_global_taxa), silent = TRUE)
            sig_pair_dat <- try(list_to_dat(sig_pair), silent = TRUE)
            
            output$taxa_barPanel = renderUI({ 
              navbarPage(
                title = "Dunn’s test (Pairwise Comparison)",
                tabPanel("Phylum", dataTableOutput("phylum_list")),
                tabPanel("Class", dataTableOutput("class_list")),
                tabPanel("Order", dataTableOutput("order_list")),
                tabPanel("Family", dataTableOutput("family_list")),
                tabPanel("Genus", dataTableOutput("genus_list")),
                position = "static-top"
              )
            })
            
            length_each_page <- nrow(as.matrix(sig_pair_dat[[1]]))/length(sig_pair[[1]])
            
            length_page_1 <- nrow(as.matrix(sig_pair_dat[[1]]))/length(sig_pair[[1]])
            length_page_2 <-  nrow(as.matrix(sig_pair_dat[[2]]))/length(sig_pair[[2]])
            length_page_3 <- nrow(as.matrix(sig_pair_dat[[3]]))/length(sig_pair[[3]])
            length_page_4 <- nrow(as.matrix(sig_pair_dat[[4]]))/length(sig_pair[[4]])
            length_page_5 <- nrow(as.matrix(sig_pair_dat[[5]]))/length(sig_pair[[5]])
            
            for (i in 1:length(sig_pair_dat)){
              colnames(sig_pair_dat[[i]])[colnames(sig_pair_dat[[i]]) == "Adj.P.value"] = "Adj. P.value"
              if("P.value" %in% colnames(sig_pair_dat[[i]])){
                ind <- which(colnames(sig_pair_dat[[i]]) == "P.value")
                sig_pair_dat[[i]] <- sig_pair_dat[[i]][, -c(ind)]
              }
            }
            
            output$phylum_list <- try(renderDataTable({datatable(sig_pair_dat[[1]], rownames = FALSE, options = list(pageLength = length_page_1, dom = '<"top" p>'))}), silent = TRUE)
            output$class_list <- try(renderDataTable({datatable(sig_pair_dat[[2]], rownames = FALSE, options = list(pageLength =length_page_2, dom = '<"top" p>'))}), silent = TRUE)
            output$order_list <- try(renderDataTable({datatable(sig_pair_dat[[3]], rownames = FALSE, options = list(pageLength = length_page_3, dom = '<"top" p>'))}), silent = TRUE)
            output$family_list <- try(renderDataTable({datatable(sig_pair_dat[[4]], rownames = FALSE, options = list(pageLength = length_page_4, dom = '<"top" p>'))}), silent = TRUE)
            output$genus_list <- try(renderDataTable({datatable(sig_pair_dat[[5]], rownames = FALSE, options = list(pageLength = length_page_5, dom = '<"top" p>'))}), silent = TRUE)
          
            output$taxa_downloadTable = renderUI({
              tagList(
                p(" ", style = "margin-bottom: +20px;"),
                box(title = strong("Download Output Table", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                    p("You can download the data analysis outputs.",
                      style = "font-size:11pt"),
                    h5("Data Analysis Outputs"),
                    downloadButton("tdownloadTabl", "Download", width = '50%', style = "color:black; background-color: red3")
                )
              )
            })
            
            output$tdownloadTabl <- downloadHandler(
              filename = function() {
                paste("Taxa.Analysis.Output.zip")
              }, 
              content = function(DA.file) {
                temp <- setwd(tempdir())
                on.exit(setwd(temp))
                dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt")
                
                capture.output(list(Global = as.data.frame(taxa.kruskal.out[[1]]), Pairwise = pairwise_result[[1]]), file = "Phylum.txt")
                capture.output(list(Global = as.data.frame(taxa.kruskal.out[[2]]), Pairwise = pairwise_result[[2]]), file = "Class.txt")
                capture.output(list(Global = as.data.frame(taxa.kruskal.out[[3]]), Pairwise = pairwise_result[[3]]), file = "Order.txt")
                capture.output(list(Global = as.data.frame(taxa.kruskal.out[[4]]), Pairwise = pairwise_result[[4]]), file = "Family.txt")
                capture.output(list(Global = as.data.frame(taxa.kruskal.out[[5]]), Pairwise = pairwise_result[[5]]), file = "Genus.txt")
                zip(zipfile=DA.file, files=dataFiles)
              }
            )
            
            }else{
            incProgress(3/10, message = "Displaying Results in progress")
            
            output$taxa_display = renderUI({
              tagList(
                tabBox(title = strong("Kruskal-Wallis test (Global Test)", style = "color:black", side = "right"), width = NULL,
                       tabPanel("Phylum", align = "center",
                                plotOutput("rank1", height = num_row[1]*250, width = width_boxplot),
                       )
                       ,
                       tabPanel("Class", align = "center",
                                plotOutput("rank2", height = num_row[2]*250, width = width_boxplot),
                       )
                       ,tabPanel("Order", align = "center",
                                 plotOutput("rank3", height = num_row[3]*250, width = width_boxplot),
                       )
                       ,tabPanel("Family", align = "center",
                                 plotOutput("rank4", height = num_row[4]*250, width = width_boxplot),
                       )
                       ,tabPanel("Genus", align = "center",
                                 plotOutput("rank5", height = num_row[5]*250, width = width_boxplot),
                       )
                       ,tabPanel("Species", align = "center",
                                 plotOutput("rank6", height = num_row[6]*250, width = width_boxplot),
                       )
                )
              )
            })
            
            output$rank1 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 1, global.p.val.only)  #
            }), silent = TRUE)
            
            output$rank2 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 2, global.p.val.only)
            }), silent = TRUE)
            
            output$rank3 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 3, global.p.val.only)
            }), silent = TRUE)
            
            output$rank4 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 4, global.p.val.only)
            }), silent = TRUE)
            
            output$rank5 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 5, global.p.val.only)
            }), silent = TRUE)
            
            output$rank6 = try(renderPlot({ 
              taxa.bin.mult (taxa.mani, sam_dat_taxa, input$taxa_prim, 6, global.p.val.only)
            }), silent = TRUE)
            
            sig_global_taxa <<- try(sig.taxa(global.p.val.only), silent = TRUE)
            sig_pair <<- try(rank_sig(pairwise_result, sig_global_taxa), silent = TRUE)
            sig_pair_dat <<- try(list_to_dat(sig_pair), silent = TRUE)
            
            output$taxa_barPanel = renderUI({ 
              navbarPage(
                title = "Dunn’s test (Pairwise Comparison)",
                tabPanel("Phylum", dataTableOutput("phylum_list")),
                tabPanel("Class", dataTableOutput("class_list")),
                tabPanel("Order", dataTableOutput("order_list")),
                tabPanel("Family", dataTableOutput("family_list")),
                tabPanel("Genus", dataTableOutput("genus_list")),
                tabPanel("Species", dataTableOutput("species_list")),
                position = "static-top"
              )
            })
            
            # print(nrow(as.matrix(sig_pair_dat[[1]])))
            length_each_page <- nrow(as.matrix(sig_pair_dat[[1]]))/length(sig_pair[[1]])
            
            length_page_1 <- nrow(as.matrix(sig_pair_dat[[1]]))/length(sig_pair[[1]])
            length_page_2 <-  nrow(as.matrix(sig_pair_dat[[2]]))/length(sig_pair[[2]])
            length_page_3 <- nrow(as.matrix(sig_pair_dat[[3]]))/length(sig_pair[[3]])
            length_page_4 <- nrow(as.matrix(sig_pair_dat[[4]]))/length(sig_pair[[4]])
            length_page_5 <- nrow(as.matrix(sig_pair_dat[[5]]))/length(sig_pair[[5]])
            length_page_6 <- nrow(as.matrix(sig_pair_dat[[6]]))/length(sig_pair[[6]])
            
            for (i in 1:length(sig_pair_dat)){
              colnames(sig_pair_dat[[i]])[colnames(sig_pair_dat[[i]]) == "Adj.P.value"] = "Adj. P.value"
              if("P.value" %in% colnames(sig_pair_dat[[i]])){
                ind <- which(colnames(sig_pair_dat[[i]]) == "P.value")
                sig_pair_dat[[i]] <- sig_pair_dat[[i]][, -c(ind)]
              }
            }
            
            output$phylum_list <- try(renderDataTable({datatable(sig_pair_dat[[1]], rownames = FALSE, options = list(pageLength = length_page_1, dom = '<"top" p>'))}), silent = TRUE)
            output$class_list <- try(renderDataTable({datatable(sig_pair_dat[[2]], rownames = FALSE, options = list(pageLength =length_page_2, dom = '<"top" p>'))}), silent = TRUE)
            output$order_list <- try(renderDataTable({datatable(sig_pair_dat[[3]], rownames = FALSE, options = list(pageLength = length_page_3, dom = '<"top" p>'))}), silent = TRUE)
            output$family_list <- try(renderDataTable({datatable(sig_pair_dat[[4]], rownames = FALSE, options = list(pageLength = length_page_4, dom = '<"top" p>'))}), silent = TRUE)
            output$genus_list <- try(renderDataTable({datatable(sig_pair_dat[[5]], rownames = FALSE, options = list(pageLength = length_page_5, dom = '<"top" p>'))}), silent = TRUE)
            output$species_list <- try(renderDataTable({datatable(sig_pair_dat[[6]], rownames = FALSE, options = list(pageLength = length_page_6, dom = '<"top" p>'))}), silent = TRUE)
            
            output$taxa_downloadTable = renderUI({
              tagList(
                p(" ", style = "margin-bottom: +20px;"),
                box(title = strong("Download Output Table", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                    p("You can download the data analysis outputs.",
                      style = "font-size:11pt"),
                    h5("Data Analysis Outputs"),
                    downloadButton("tdownloadTabl", "Download", width = '50%', style = "color:black; background-color: red3")
                )
              )
            })
            
            output$tdownloadTabl <- downloadHandler(
              filename = function() {
                paste("Taxa.Analysis.Output.zip")
              }, 
              content = function(DA.file) {
                temp <- setwd(tempdir())
                on.exit(setwd(temp))
                dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
                
                capture.output(list(Global = as.data.frame(taxa.kruskal.out[[1]]), Pairwise = pairwise_result[[1]]), file = "Phylum.txt")
                capture.output(list(Global = as.data.frame(taxa.kruskal.out[[2]]), Pairwise = pairwise_result[[2]]), file = "Class.txt")
                capture.output(list(Global = as.data.frame(taxa.kruskal.out[[3]]), Pairwise = pairwise_result[[3]]), file = "Order.txt")
                capture.output(list(Global = as.data.frame(taxa.kruskal.out[[4]]), Pairwise = pairwise_result[[4]]), file = "Family.txt")
                capture.output(list(Global = as.data.frame(taxa.kruskal.out[[5]]), Pairwise = pairwise_result[[5]]), file = "Genus.txt")
                capture.output(list(Global = as.data.frame(taxa.kruskal.out[[6]]), Pairwise = pairwise_result[[6]]), file = "Species.txt")
                zip(zipfile=DA.file, files=dataFiles)
              }
            )
            
          }
        }
        
        ## Multinomial Logistic Regression ---------
        
        else if (input$taxa_method == "Multinomial Logistic Regression"){
          
          shinyjs::hide("taxa_display_area")
          shinyjs::hide("taxa_barPanel")
          shinyjs::show("taxa_display_global")
          shinyjs::hide("taxa_display")
          shinyjs::show("taxa_display_forest")
          
          # Global Test
          taxa.nominal.global.summary.list <- try(taxa.nominal.global(taxa.mani, level.names, covariate = covariate), silent = TRUE)
          taxa.nominal.global.table <- try(taxa.nominal.global.outcome(taxa.nominal.global.summary.list, level.names), silent = TRUE)
          
          # Pairwise comparison
          
          ctable.list <- try(taxa.nominal.result.merge(taxa.mani, level.names, num_cat, covariate), silent = TRUE)
          nlr <- num_cat - 1
          taxa_thres <- 0.05
          
          forest_nrow <- c()
          height_forest <- c()
          num_comp <- c()
          rclm <- list()
          
          for(i in 1:nlr){
            ctable <- try(taxa.forest.plot.pages1(ctable.list[[i]], thres = taxa_thres, mult.test.cor = TRUE), silent = TRUE)
            if(is.null(ctable)){
              rclm[i] <- list(NULL)
            }
            else{
              rclm[[i]] <- ctable
            }
            
            forest_nrow <- try(c(forest_nrow, taxa.forest.plot.pages(ctable.list[[i]], thres = taxa_thres, mult.test.cor = TRUE)), silent = TRUE)
            
            if(is.null(rclm[[i]])){ # No significant taxon
              height_forest <- c(height_forest, 200)
            }
            else if(nrow(rclm[[i]][[1]]) == 1){
              height_forest <- c(height_forest, 150)
            }
            else if(nrow(rclm[[i]][[1]]) > 1 & nrow(rclm[[i]][[1]]) < 5){
              height_forest <- c(height_forest, 300)
            }
            else if(nrow(rclm[[i]][[1]]) >= 5 & nrow(rclm[[i]][[1]]) < 20){
              height_forest <- c(height_forest, 500)
            }
            else {
              height_forest <- c(height_forest, 800)
            }
            
            num_comp <- try(c(num_comp, ifelse(is.null(rclm[[i]]), 1, length(rclm[[i]]))), silent = TRUE)
          }
          
          output$taxa_display_global <- renderUI({
            tagList(
              box(title = strong("Multinomial Logistic Regression (Global Test)", style = "color:white"),  
                  align = "center", width = NULL, solidHeader = TRUE, status = "info", br(),
                  dataTableOutput("taxa_nominal_global_table", height = "auto", width = 750)
              )
            )
          })
          
          output$taxa_nominal_global_table <- try(renderDataTable({
            taxa.nominal.global.table
          }), silent = TRUE)
          
          if(num_cat == 3){
            output$taxa_display_forest = renderUI({
              tagList(
                box(title = strong("Multinomial Logistic Regression (Pairwise Comparison)", style = "color:white"),  
                    align = "center", width = NULL, solidHeader = TRUE, status = "info", 
                    tabsetPanel(tabPanel(title = paste0("Reference/", names(sort(taxa.cat.order)[-1])[1]), align = "center", 
                                         do.call(tabsetPanel, lapply(1:forest_nrow[1], function(i){
                                           tabPanel(title = paste0("Page", i), align = "center", br(),
                                                    plotOutput(paste0("taxa_forest_plot_cat3_comp1", i), height = height_forest[1], width = 750))
                                         })))),
                    tabsetPanel(tabPanel(title = paste0("Reference/", names(sort(taxa.cat.order)[-1])[2]), align = "center", 
                                         do.call(tabsetPanel, lapply(1:forest_nrow[2], function(i){
                                           tabPanel(title = paste0("Page", i), align = "center", br(),
                                                    plotOutput(paste0("taxa_forest_plot_cat3_comp2", i), height = height_forest[2], width = 750))
                                         }))))
                )
              )
            })
            
            lapply(1:forest_nrow[1], function(i) {
              output[[paste0("taxa_forest_plot_cat3_comp1", i)]] <- renderPlot({
                tryCatch(taxa.forest.plot(rclm[[1]], page = i), error = function(e){
                  message("No outcome is available!")
                  showModal(modalDialog(div("No outcome is available!")))
                  return(NULL)})
              })
            })
            
            lapply(1:forest_nrow[2], function(i) {
              output[[paste0("taxa_forest_plot_cat3_comp2", i)]] <- renderPlot({
                tryCatch(taxa.forest.plot(rclm[[2]], page = i), error = function(e){
                  message("No outcome is available!")
                  showModal(modalDialog(div("No outcome is available!")))
                  return(NULL)})
              })
            })
          } else if(num_cat == 4){
            output$taxa_display_forest = renderUI({
              tagList(
                box(title = strong("Multinomial Logistic Regression (Pairwise Comparison)", style = "color:white"),  
                    align = "center", width = NULL, solidHeader = TRUE, status = "info", 
                    tabsetPanel(tabPanel(title = paste0("Reference/", names(sort(taxa.cat.order)[-1])[1]), align = "center", 
                                         do.call(tabsetPanel, lapply(1:forest_nrow[1], function(i){
                                           tabPanel(title = paste0("Page", i), align = "center", br(),
                                                    plotOutput(paste0("taxa_forest_plot_cat4_comp1", i), height = height_forest[1], width = 750))
                                         })))),
                    tabsetPanel(tabPanel(title = paste0("Reference/", names(sort(taxa.cat.order)[-1])[2]), align = "center", 
                                         do.call(tabsetPanel, lapply(1:forest_nrow[2], function(i){
                                           tabPanel(title = paste0("Page", i), align = "center", br(),
                                                    plotOutput(paste0("taxa_forest_plot_cat4_comp2", i), height = height_forest[2], width = 750))
                                         })))),
                    tabsetPanel(tabPanel(title = paste0("Reference/", names(sort(taxa.cat.order)[-1])[3]), align = "center", 
                                         do.call(tabsetPanel, lapply(1:forest_nrow[3], function(i){
                                           tabPanel(title = paste0("Page", i), align = "center", br(),
                                                    plotOutput(paste0("taxa_forest_plot_cat4_comp3", i), height = height_forest[3], width = 750))
                                         }))))
                )
              )
            })
            
            lapply(1:forest_nrow[1], function(i) {
              output[[paste0("taxa_forest_plot_cat4_comp1", i)]] <- renderPlot({
                tryCatch(taxa.forest.plot(rclm[[1]], page = i), error = function(e){
                  message("No outcome is available!")
                  showModal(modalDialog(div("No outcome is available!")))
                  return(NULL)})
              })
            })
            
            lapply(1:forest_nrow[2], function(i) {
              output[[paste0("taxa_forest_plot_cat4_comp2", i)]] <- renderPlot({
                tryCatch(taxa.forest.plot(rclm[[2]], page = i), error = function(e){
                  message("No outcome is available!")
                  showModal(modalDialog(div("No outcome is available!")))
                  return(NULL)})
              })
            })
            
            lapply(1:forest_nrow[3], function(i) {
              output[[paste0("taxa_forest_plot_cat4_comp3", i)]] <- renderPlot({
                tryCatch(taxa.forest.plot(rclm[[3]], page = i), error = function(e){
                  message("No outcome is available!")
                  showModal(modalDialog(div("No outcome is available!")))
                  return(NULL)})
              })
            })
          } else if(num_cat == 5){
            output$taxa_display_forest = renderUI({
              tagList(
                box(title = strong("Multinomial Logistic Regression (Pairwise Comparison)", style = "color:white"),  
                    align = "center", width = NULL, solidHeader = TRUE, status = "info", 
                    tabsetPanel(tabPanel(title = paste0("Reference/", names(sort(taxa.cat.order)[-1])[1]), align = "center", 
                                         do.call(tabsetPanel, lapply(1:forest_nrow[1], function(i){
                                           tabPanel(title = paste0("Page", i), align = "center", br(),
                                                    plotOutput(paste0("taxa_forest_plot_cat5_comp1", i), height = height_forest[1], width = 750))
                                         })))),
                    tabsetPanel(tabPanel(title = paste0("Reference/", names(sort(taxa.cat.order)[-1])[2]), align = "center", 
                                         do.call(tabsetPanel, lapply(1:forest_nrow[2], function(i){
                                           tabPanel(title = paste0("Page", i), align = "center", br(),
                                                    plotOutput(paste0("taxa_forest_plot_cat5_comp2", i), height = height_forest[2], width = 750))
                                         })))),
                    tabsetPanel(tabPanel(title = paste0("Reference/", names(sort(taxa.cat.order)[-1])[3]), align = "center", 
                                         do.call(tabsetPanel, lapply(1:forest_nrow[3], function(i){
                                           tabPanel(title = paste0("Page", i), align = "center", br(),
                                                    plotOutput(paste0("taxa_forest_plot_cat5_comp3", i), height = height_forest[3], width = 750))
                                         })))),
                    tabsetPanel(tabPanel(title = paste0("Reference/", names(sort(taxa.cat.order)[-1])[4]), align = "center", 
                                         do.call(tabsetPanel, lapply(1:forest_nrow[4], function(i){
                                           tabPanel(title = paste0("Page", i), align = "center", br(),
                                                    plotOutput(paste0("taxa_forest_plot_cat5_comp4", i), height = height_forest[4], width = 750))
                                         }))))
                )
              )
            })
            
            lapply(1:forest_nrow[1], function(i) {
              output[[paste0("taxa_forest_plot_cat5_comp1", i)]] <- renderPlot({
                tryCatch(taxa.forest.plot(rclm[[1]], page = i), error = function(e){
                  message("No outcome is available!")
                  showModal(modalDialog(div("No outcome is available!")))
                  return(NULL)})
              })
            })
            
            lapply(1:forest_nrow[2], function(i) {
              output[[paste0("taxa_forest_plot_cat5_comp2", i)]] <- renderPlot({
                tryCatch(taxa.forest.plot(rclm[[2]], page = i), error = function(e){
                  message("No outcome is available!")
                  showModal(modalDialog(div("No outcome is available!")))
                  return(NULL)})
              })
            })
            
            lapply(1:forest_nrow[3], function(i) {
              output[[paste0("taxa_forest_plot_cat5_comp3", i)]] <- renderPlot({
                tryCatch(taxa.forest.plot(rclm[[3]], page = i), error = function(e){
                  message("No outcome is available!")
                  showModal(modalDialog(div("No outcome is available!")))
                  return(NULL)})
              })
            })
            
            lapply(1:forest_nrow[4], function(i) {
              output[[paste0("taxa_forest_plot_cat5_comp4", i)]] <- renderPlot({
                tryCatch(taxa.forest.plot(rclm[[4]], page = i), error = function(e){
                  message("No outcome is available!")
                  showModal(modalDialog(div("No outcome is available!")))
                  return(NULL)})
              })
            })
          } else if(num_cat == 6){
            output$taxa_display_forest = renderUI({
              tagList(
                box(title = strong("Multinomial Logistic Regression (Pairwise Comparison)", style = "color:white"),  
                    align = "center", width = NULL, solidHeader = TRUE, status = "info", 
                    tabsetPanel(tabPanel(title = paste0("Reference/", names(sort(taxa.cat.order)[-1])[1]), align = "center", 
                                         do.call(tabsetPanel, lapply(1:forest_nrow[1], function(i){
                                           tabPanel(title = paste0("Page", i), align = "center", br(),
                                                    plotOutput(paste0("taxa_forest_plot_cat6_comp1", i), height = height_forest[1], width = 750))
                                         })))),
                    tabsetPanel(tabPanel(title = paste0("Reference/", names(sort(taxa.cat.order)[-1])[2]), align = "center", 
                                         do.call(tabsetPanel, lapply(1:forest_nrow[2], function(i){
                                           tabPanel(title = paste0("Page", i), align = "center", br(),
                                                    plotOutput(paste0("taxa_forest_plot_cat6_comp2", i), height = height_forest[2], width = 750))
                                         })))),
                    tabsetPanel(tabPanel(title = paste0("Reference/", names(sort(taxa.cat.order)[-1])[3]), align = "center", 
                                         do.call(tabsetPanel, lapply(1:forest_nrow[3], function(i){
                                           tabPanel(title = paste0("Page", i), align = "center", br(),
                                                    plotOutput(paste0("taxa_forest_plot_cat6_comp3", i), height = height_forest[3], width = 750))
                                         })))),
                    tabsetPanel(tabPanel(title = paste0("Reference/", names(sort(taxa.cat.order)[-1])[4]), align = "center", 
                                         do.call(tabsetPanel, lapply(1:forest_nrow[4], function(i){
                                           tabPanel(title = paste0("Page", i), align = "center", br(),
                                                    plotOutput(paste0("taxa_forest_plot_cat6_comp4", i), height = height_forest[4], width = 750))
                                         })))),
                    tabsetPanel(tabPanel(title = paste0("Reference/", names(sort(taxa.cat.order)[-1])[5]), align = "center", 
                                         do.call(tabsetPanel, lapply(1:forest_nrow[5], function(i){
                                           tabPanel(title = paste0("Page", i), align = "center", br(),
                                                    plotOutput(paste0("taxa_forest_plot_cat6_comp5", i), height = height_forest[5], width = 750))
                                         }))))
                )
              )
            })
            
            lapply(1:forest_nrow[1], function(i) {
              output[[paste0("taxa_forest_plot_cat6_comp1", i)]] <- renderPlot({
                tryCatch(taxa.forest.plot(rclm[[1]], page = i), error = function(e){
                  message("No outcome is available!")
                  showModal(modalDialog(div("No outcome is available!")))
                  return(NULL)})
              })
            })
            
            lapply(1:forest_nrow[2], function(i) {
              output[[paste0("taxa_forest_plot_cat6_comp2", i)]] <- renderPlot({
                tryCatch(taxa.forest.plot(rclm[[2]], page = i), error = function(e){
                  message("No outcome is available!")
                  showModal(modalDialog(div("No outcome is available!")))
                  return(NULL)})
              })
            })
            
            lapply(1:forest_nrow[3], function(i) {
              output[[paste0("taxa_forest_plot_cat6_comp3", i)]] <- renderPlot({
                tryCatch(taxa.forest.plot(rclm[[3]], page = i), error = function(e){
                  message("No outcome is available!")
                  showModal(modalDialog(div("No outcome is available!")))
                  return(NULL)})
              })
            })
            
            lapply(1:forest_nrow[4], function(i) {
              output[[paste0("taxa_forest_plot_cat6_comp4", i)]] <- renderPlot({
                tryCatch(taxa.forest.plot(rclm[[4]], page = i), error = function(e){
                  message("No outcome is available!")
                  showModal(modalDialog(div("No outcome is available!")))
                  return(NULL)})
              })
            })
            
            lapply(1:forest_nrow[5], function(i) {
              output[[paste0("taxa_forest_plot_cat6_comp5", i)]] <- renderPlot({
                tryCatch(taxa.forest.plot(rclm[[5]], page = i), error = function(e){
                  message("No outcome is available!")
                  showModal(modalDialog(div("No outcome is available!")))
                  return(NULL)})
              })
            })
          } else if(num_cat == 7){
            output$taxa_display_forest = renderUI({
              tagList(
                box(title = strong("Multinomial Logistic Regression (Pairwise Comparison)", style = "color:white"),  
                    align = "center", width = NULL, solidHeader = TRUE, status = "info", 
                    tabsetPanel(tabPanel(title = paste0("Reference/", names(sort(taxa.cat.order)[-1])[1]), align = "center", 
                                         do.call(tabsetPanel, lapply(1:forest_nrow[1], function(i){
                                           tabPanel(title = paste0("Page", i), align = "center", br(),
                                                    plotOutput(paste0("taxa_forest_plot_cat7_comp1", i), height = height_forest[1], width = 750))
                                         })))),
                    tabsetPanel(tabPanel(title = paste0("Reference/", names(sort(taxa.cat.order)[-1])[2]), align = "center", 
                                         do.call(tabsetPanel, lapply(1:forest_nrow[2], function(i){
                                           tabPanel(title = paste0("Page", i), align = "center", br(),
                                                    plotOutput(paste0("taxa_forest_plot_cat7_comp2", i), height = height_forest[2], width = 750))
                                         })))),
                    tabsetPanel(tabPanel(title = paste0("Reference/", names(sort(taxa.cat.order)[-1])[3]), align = "center", 
                                         do.call(tabsetPanel, lapply(1:forest_nrow[3], function(i){
                                           tabPanel(title = paste0("Page", i), align = "center", br(),
                                                    plotOutput(paste0("taxa_forest_plot_cat7_comp3", i), height = height_forest[3], width = 750))
                                         })))),
                    tabsetPanel(tabPanel(title = paste0("Reference/", names(sort(taxa.cat.order)[-1])[4]), align = "center", 
                                         do.call(tabsetPanel, lapply(1:forest_nrow[4], function(i){
                                           tabPanel(title = paste0("Page", i), align = "center", br(),
                                                    plotOutput(paste0("taxa_forest_plot_cat7_comp4", i), height = height_forest[4], width = 750))
                                         })))),
                    tabsetPanel(tabPanel(title = paste0("Reference/", names(sort(taxa.cat.order)[-1])[5]), align = "center", 
                                         do.call(tabsetPanel, lapply(1:forest_nrow[5], function(i){
                                           tabPanel(title = paste0("Page", i), align = "center", br(),
                                                    plotOutput(paste0("taxa_forest_plot_cat7_comp5", i), height = height_forest[5], width = 750))
                                         })))),
                    tabsetPanel(tabPanel(title = paste0("Reference/", names(sort(taxa.cat.order)[-1])[6]), align = "center", 
                                         do.call(tabsetPanel, lapply(1:forest_nrow[6], function(i){
                                           tabPanel(title = paste0("Page", i), align = "center", br(),
                                                    plotOutput(paste0("taxa_forest_plot_cat7_comp6", i), height = height_forest[6], width = 750))
                                         }))))
                )
              )
            })
            
            lapply(1:forest_nrow[1], function(i) {
              output[[paste0("taxa_forest_plot_cat7_comp1", i)]] <- renderPlot({
                tryCatch(taxa.forest.plot(rclm[[1]], page = i), error = function(e){
                  message("No outcome is available!")
                  showModal(modalDialog(div("No outcome is available!")))
                  return(NULL)})
              })
            })
            
            lapply(1:forest_nrow[2], function(i) {
              output[[paste0("taxa_forest_plot_cat7_comp2", i)]] <- renderPlot({
                tryCatch(taxa.forest.plot(rclm[[2]], page = i), error = function(e){
                  message("No outcome is available!")
                  showModal(modalDialog(div("No outcome is available!")))
                  return(NULL)})
              })
            })
            
            lapply(1:forest_nrow[3], function(i) {
              output[[paste0("taxa_forest_plot_cat7_comp3", i)]] <- renderPlot({
                tryCatch(taxa.forest.plot(rclm[[3]], page = i), error = function(e){
                  message("No outcome is available!")
                  showModal(modalDialog(div("No outcome is available!")))
                  return(NULL)})
              })
            })
            
            lapply(1:forest_nrow[4], function(i) {
              output[[paste0("taxa_forest_plot_cat7_comp4", i)]] <- renderPlot({
                tryCatch(taxa.forest.plot(rclm[[4]], page = i), error = function(e){
                  message("No outcome is available!")
                  showModal(modalDialog(div("No outcome is available!")))
                  return(NULL)})
              })
            })
            
            lapply(1:forest_nrow[5], function(i) {
              output[[paste0("taxa_forest_plot_cat7_comp5", i)]] <- renderPlot({
                tryCatch(taxa.forest.plot(rclm[[5]], page = i), error = function(e){
                  message("No outcome is available!")
                  showModal(modalDialog(div("No outcome is available!")))
                  return(NULL)})
              })
            })
            
            lapply(1:forest_nrow[6], function(i) {
              output[[paste0("taxa_forest_plot_cat7_comp6", i)]] <- renderPlot({
                tryCatch(taxa.forest.plot(rclm[[6]], page = i), error = function(e){
                  message("No outcome is available!")
                  showModal(modalDialog(div("No outcome is available!")))
                  return(NULL)})
              })
            })
          }
          
          output$taxa_references <- renderUI({
            tagList(
              p(" ", style = "margin-top: 20px;"),
              box(title = strong("References", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  p(ALPHA_MULTINOM_REFERENCE, style = "font-size:11pt")
              )
            )
          })
          
          output$taxa_downloadTable = renderUI({
            tagList(
              p(" ", style = "margin-bottom: +20px;"),
              box(title = strong("Download Output Table", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  p("You can download the data analysis outputs.",
                    style = "font-size:11pt"),
                  h5("Data Analysis Outputs"),
                  downloadButton("tdownloadTabl", "Download", width = '50%', style = "color:black; background-color: red3")
              )
            )
          })
          
          output$tdownloadTabl <- downloadHandler(
            filename = function() {
              paste("Taxa.Analysis.Output.txt")
            },
            content = function(file) {
              sink(file); print(ctable.list); sink()
            }
          )
        }
        
        ## Proportional Odds Model ---------
        
        else if (input$taxa_method == "Proportional Odds Model"){
          
          shinyjs::hide("taxa_display_area")
          shinyjs::hide("taxa_barPanel")
          shinyjs::hide("taxa_display_global")
          shinyjs::hide("taxa_display")
          shinyjs::show("taxa_display_forest")
          
          taxa_thres <- 0.05
          
          result.ctable <- try(taxa.ordinal.result.merge(taxa.mani, level.names, covariate = covariate), silent = TRUE)
          
          forest_nrow <- try(taxa.forest.plot.pages(result.ctable, thres = taxa_thres, mult.test.cor = TRUE), silent = TRUE)
          
          result.ctable.list <- try(taxa.forest.plot.pages1(result.ctable, thres = taxa_thres, mult.test.cor = TRUE), silent = TRUE)
          
          height_forest <- 0
          
          if(is.null(result.ctable.list[[1]])){ # No significant taxon
            height_forest <- 200
          }
          else if(nrow(result.ctable.list[[1]]) == 1){
            height_forest <- 150
          }
          else if(nrow(result.ctable.list[[1]]) > 1 & nrow(result.ctable.list[[1]]) < 5){
            height_forest <- 300
          }
          else if(nrow(result.ctable.list[[1]]) >= 5 & nrow(result.ctable.list[[1]]) < 20){
            height_forest <- 500
          }
          else {
            height_forest <- 800
          }
          
          output$taxa_display_forest = renderUI({
            tagList(
              box(title = strong("Proportional Odds Model", style = "color:white"),  
                  align = "center", width = NULL, solidHeader = TRUE, status = "info", 
                  do.call(tabsetPanel, lapply(1:forest_nrow, function(i){
                    tabPanel(title = paste0("Page", i), align = "center",
                             plotOutput(paste0("taxa_ordinal_forest_plot", i), height = height_forest, width = 750))
                  }))
              )
            )
          })
          
          lapply(1:forest_nrow, function(j) {
            output[[paste0("taxa_ordinal_forest_plot", j)]] <- try(renderPlot({
              taxa.forest.plot(result.ctable.list, page = j)
            }), silent = TRUE)
          })
          
          output$taxa_references <- renderUI({
            tagList(
              p(" ", style = "margin-top: 20px;"),
              box(title = strong("References", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  p(ALPHA_PROPODDS_REFERENCE, style = "font-size:11pt")
              )
            )
          })
          
          output$taxa_downloadTable = renderUI({
            tagList(
              p(" ", style = "margin-bottom: +20px;"),
              box(title = strong("Download Output Table", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  p("You can download the data analysis outputs.",
                    style = "font-size:11pt"),
                  h5("Data Analysis Outputs"),
                  downloadButton("tdownloadTabl", "Download", width = '50%', style = "color:black; background-color: red3")
              )
            )
          })
          
          output$tdownloadTabl <- downloadHandler(
            filename = function() {
              paste("Taxa.Analysis.Output.txt")
            },
            content = function(file) {
              sink(file); print(result.ctable); sink()
            }
          )
        }
      })
    
    for (i in 1:prim_length){
      shinyjs::enable(paste0("taxaCat_", i))
    }
    shinyjs::enable("taxa_primvars")
    shinyjs::enable("taxa_primvars_ref")
    shinyjs::enable("taxa_primvars_rename")
    shinyjs::enable("taxa_var_type")
    shinyjs::enable("taxa_covariate")
    shinyjs::enable("taxa_method")
    shinyjs::enable("include_species_taxa")
    shinyjs::enable("taxa_downloadTabUI")
    shinyjs::enable("taxa_references")
    
  })
  
  ## 5-1. RF -------------------
  
  observeEvent(input$rf_nom_runButton, {
    shinyjs::disable("rf_nom_dataType")
    shinyjs::disable("rf_nom_response")
    shinyjs::disable("rf_nom_data_input_opt")
    shinyjs::disable("rf_nom_nfold")
    shinyjs::disable("rf_nom_ntree")
    shinyjs::disable("rf_nom_var_num")
    shinyjs::disable("rf_nom_include_species")
    shinyjs::disable("rf_nom_runButton")
    
    withProgress(
      message = "Calculation in progress",
      detail = "This may take a while...", value = 0, {
        if (input$rf_nom_dataType == "Count (Rarefied)") {
          type = "rare.count"
          ref = RF_REFERENCE_RC
        }
        else if (input$rf_nom_dataType == "Proportion") {
          type = "prop"
          ref = RF_REFERENCE
        }
        else if (input$rf_nom_dataType == "CLR (Default)") {
          type = "clr"
          ref = RF_REFERENCE_CLR
        }
        else if(input$rf_nom_dataType == "Arcsine-root"){
          type = "arcsin"
          ref = RF_REFERENCE
        }

        if(input$rf_nom_include_species == "Phylum - Genus (16S)"){
          level.names = get.level.names(include = FALSE)
        }
        else if(input$rf_nom_include_species == "Phylum - Species (Metagenomics)"){
          level.names = get.level.names(include = TRUE)
        }
        
        taxa.out <- chooseData$taxa.out

        for(j in 1:6){
          for(name in level.names){
            names(taxa.out[[j]][[name]]) <- chooseData$taxa.names.out$names[[name]]
          }
        }
        data <- taxa.out[[type]]
        rf.colnames.list.mult <- colnames.to.ind(data)
        data <- change.colnames(data, rf.colnames.list.mult$new)

        input.data <- remove.na(data = data, sam.dat = chooseData$sam.dat, y.name = input$rf_nom_response, level.names = level.names)
        data <- input.data[[1]]
        sam.dat.na <- input.data[[2]]
        
        y.name <- input$rf_nom_response
        nfold <- as.numeric(input$rf_nom_nfold)
        ntree <- as.numeric(input$rf_nom_ntree)
        
        rf.list <- list()
        for(name in level.names){
          incProgress(1/10, message = sprintf("Random Forest: %s in progress", str_to_title(name)))
          set.seed(578)
          rf.list[[name]] <- try(rf.cla.rev(data = data,
                                            sam.dat.na = sam.dat.na,
                                            y.name = y.name,
                                            nfold = nfold,
                                            ntree = ntree,
                                            stratified = TRUE,
                                            name = name,
                                            p = 1), 
                                 silent = TRUE)
        }

        incProgress(2/10, message = "Visualizations in progress")

        rf.imp.plot.list <- list()
        pd.plot.list <- list()
        
        rf_nom_length <- length(category.names(sam.dat.na, y.name))
        rf_nom_name <- c()
        for(i in 1:rf_nom_length){
          rf_nom_name <- c(rf_nom_name, input[[paste0("rf_nom_label", i)]])
        }
        
        for(name in level.names){
          rf.imp.plot.list[[name]] <- try(rf.imp.plot(rf.list, 
                                                      name, 
                                                      n = as.numeric(input$rf_nom_var_num), 
                                                      type = 2, 
                                                      is.cat = TRUE, 
                                                      data = data, 
                                                      data.type = type), 
                                          silent = TRUE)
          pd.plot.list[[name]] <- try(rf.pdp.mult(rf.list = rf.list, 
                                                  n = as.numeric(input$rf_nom_var_num), 
                                                  name = name, 
                                                  data.type = type, 
                                                  level.names = level.names,
                                                  label = rf_nom_name), 
                                      silent = TRUE)
        }

        rf.width.list <- list()
        for(name in level.names){
          if(as.numeric(input$rf_nom_var_num) > ncol(data[[name]])){
            n <- ncol(data[[name]])
            if(n %% 5 == 0){
              n <- n / 5
            }
            else{
              n <- (n / 5) + 1
            }
          }
          else{
            n <- as.numeric(input$rf_nom_var_num)
            n <- n / 5
          }
          rf.width.list[[name]] <- n * 160
        }

        output$rf_nom_results <- renderUI({
          tagList(
            do.call(tabsetPanel, lapply(1:length(level.names), function(i) {
              tabPanel(title = str_to_title(level.names[i]), align = "center",
                       tabsetPanel(tabPanel(title = "Importance", align = "center", br(),
                                            plotOutput(paste0("rf_nom_imp", i), height = 700, width = 510),
                                            dataTableOutput(paste0("rf_nom_column_table1", i), height = "auto", width = 510)),
                                   tabPanel(title = "Partial Dependence", align = "center", br(),
                                            plotOutput(paste0("rf_nom_pd", i), height = 750, width = rf.width.list[[i]]), br(),
                                            dataTableOutput(paste0("rf_nom_column_table2", i), height = "auto", width = 510)),
                                   tabPanel(title = "CV Error", align = "center", br(),
                                            plotOutput(paste0("rf_nom_cv", i), height = 700, width = 700)),
                                   tabPanel(title = "OOB Error", align = "center", br(),
                                            plotOutput(paste0("rf_nom_oob", i), height = 700, width = 700))))}))
          )
        })

        lapply(1:length(level.names), function(j) {
          output[[paste0("rf_nom_imp", j)]] <- renderPlot({
            tryCatch(rf.imp.plot.list[[level.names[j]]], error = function(e){
              message("Visualization not available! Check the input.")
              showModal(modalDialog(div("Visualization not available! Check the input.")))
              return(NULL)
            })
          })

          output[[paste0("rf_nom_pd", j)]] <- renderPlot({
            tryCatch(pd.plot.list[[level.names[j]]], error = function(e){
              message("Visualization not available! Check the input.")
              showModal(modalDialog(div("Visualization not available! Check the input.")))
              return(NULL)
            })
          })

          output[[paste0("rf_nom_column_table1", j)]] <- renderDataTable({
            tryCatch(rf.pd.var.used(rf.list, level.names[[j]], rf.colnames.list.mult, n = as.numeric(input$rf_nom_var_num), is.cat = TRUE), error = function(e){
              message("Visualization not available! Check the input.")
              showModal(modalDialog(div("Visualization not available! Check the input.")))
              return(NULL)
            })
          })
          
          output[[paste0("rf_nom_column_table2", j)]] <- renderDataTable({
            tryCatch(rf.pd.var.used(rf.list, level.names[[j]], rf.colnames.list.mult, n = as.numeric(input$rf_nom_var_num), is.cat = TRUE), error = function(e){
              message("Visualization not available! Check the input.")
              showModal(modalDialog(div("Visualization not available! Check the input.")))
              return(NULL)
            })
          })

          output[[paste0("rf_nom_cv", j)]] <- renderPlot({
            tryCatch(cv.mtry(rf.list, level.names[j], is.cat = TRUE), error = function(e){
              message("Visualization not available! Check the input.")
              showModal(modalDialog(div("Visualization not available! Check the input.")))
              return(NULL)
            })
          })

          output[[paste0("rf_nom_oob", j)]] <- renderPlot({
            tryCatch(error.plot(rf.list, level.names[j], is.cat = TRUE), error = function(e){
              message("Visualization not available! Check the input.")
              showModal(modalDialog(div("Visualization not available! Check the input.")))
              return(NULL)
            })
          })

          output$rf_nom_downloadTabUI <- renderUI({
            tagList(
              p(" ", style = "margin-top: 20px;"),
              box(title = strong("Download Output Table", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  p("You can download the data analysis outputs.",
                    style = "font-size:11pt"),
                  h5("Feature Importance"),
                  downloadButton("rf_nom_downloadTable", "Download", width = '50%', style = "background-color: red3")
              )
            )
          })

          output$rf_nom_downloadTable <- downloadHandler(
            filename = function() {
              paste("RF_Importance.zip")
            },
            content = function(DA.file) {
              temp <- setwd(tempdir())
              on.exit(setwd(temp))
              if (length(level.names) == 5) {
                dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt")
                for(i in 1:length(dataFiles)){
                  write.table(rf.imp.df(rf.list, level.names[i], type = 0), file = dataFiles[i], sep = "\t")
                }
              }
              else if(length(level.names) == 6) {
                dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
                for(i in 1:length(dataFiles)){
                  write.table(rf.imp.df(rf.list, level.names[i], type = 0), file = dataFiles[i], sep = "\t")
                }
              }
              zip(zipfile=DA.file, files=dataFiles)
            })

          output$rf_nom_reference <- renderUI({
            tagList(
              p(" ", style = "margin-top: 20px;"),
              box(title = strong("References", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  p(ref, style = "font-size:11pt")
              )
            )
          })
        })
      }
    )
    shinyjs::enable("rf_nom_dataType")
    shinyjs::enable("rf_nom_response")
    shinyjs::enable("rf_nom_data_input_opt")
    shinyjs::enable("rf_nom_nfold")
    shinyjs::enable("rf_nom_ntree")
    shinyjs::enable("rf_nom_var_num")
    shinyjs::enable("rf_nom_include_species")
    shinyjs::enable("rf_nom_runButton")
  })
  
  ## 5-2. XGB -------------------
  
  observeEvent(input$xgb_nom_runButton, {
    shinyjs::disable("xgb_nom_dataType")
    shinyjs::disable("xgb_nom_response")
    shinyjs::disable("xgb_nom_data_input_opt")
    shinyjs::disable("xgb_nom_nfold")
    shinyjs::disable("xgb_nom_eta")
    shinyjs::disable("xgb_nom_penalty")
    shinyjs::disable("xgb_nom_nrounds")
    shinyjs::disable("xgb_nom_var_num")
    shinyjs::disable("xgb_nom_include_species")
    shinyjs::disable("xgb_nom_runButton")

    withProgress(
      message = "Calculation in progress",
      detail = "This may take a while...", value = 0, {
        if (input$xgb_nom_dataType == "Count (Rarefied)") {
          type = "rare.count"
          ref = XGB_REFERENCE_RC
        }
        else if (input$xgb_nom_dataType == "Proportion") {
          type = "prop"
          ref = XGB_REFERENCE
        }
        else if (input$xgb_nom_dataType == "CLR (Default)") {
          type = "clr"
          ref = XGB_REFERENCE_CLR
        }
        else if(input$xgb_nom_dataType == "Arcsine-root"){
          type = "arcsin"
          ref = XGB_REFERENCE
        }

        if(input$xgb_nom_penalty == "Yes (Default)"){
          alpha = 0
          lambda = 1
        }
        else if(input$xgb_nom_penalty == "No"){
          alpha = 0
          lambda = 0
        }

        if(input$xgb_nom_include_species == "Phylum - Genus (16S)"){
          level.names = get.level.names(include = FALSE)
        }
        else if(input$xgb_nom_include_species == "Phylum - Species (Metagenomics)"){
          level.names = get.level.names(include = TRUE)
        }

        taxa.out <- chooseData$taxa.out

        for(j in 1:6){
          for(name in level.names){
            names(taxa.out[[j]][[name]]) <- chooseData$taxa.names.out$names[[name]]
          }
        }
        data <- taxa.out[[type]]
        xgb.colnames.list.mult <- colnames.to.ind(data)
        data <- change.colnames(data, xgb.colnames.list.mult$new)
        
        input.data <- remove.na(data = data, sam.dat = chooseData$sam.dat, y.name = input$xgb_nom_response, level.names = level.names)
        data <- try(input.data[[1]], silent = TRUE)
        sam.dat.na <- try(input.data[[2]], silent = TRUE)
        
        y.name <- input$xgb_nom_response
        eval <- xgb.model.input.cla$eval
        eta <- as.numeric(input$xgb_nom_eta)
        nfold <- as.numeric(input$xgb_nom_nfold)
        nrounds <- as.numeric(input$xgb_nom_nrounds)
        
        if(!is.numeric(sam.dat.na[[y.name]])){
          sam.dat.na[[y.name]] <- as.factor(sam.dat.na[[y.name]])
          levels(sam.dat.na[[y.name]]) <- 0:(length(category.names(sam.dat.na, y.name))-1)
          sam.dat.na[[y.name]] <- as.numeric(sam.dat.na[[y.name]])-1
        }

        xgb.list <- list()
        for(name in level.names){
          incProgress(1/10, message = sprintf("Extreme Gradient Boosting: %s in progress", str_to_title(name)))
          set.seed(578)
          xgb.list[[name]] <- try(xgb.mult(data = data,
                                          sam.dat.na = sam.dat.na,
                                          y.name = y.name,
                                          eta = eta,
                                          nrounds = nrounds,
                                          nfold = nfold,
                                          alpha = alpha,
                                          lambda = lambda,
                                          stratified = TRUE,
                                          name = name),
                                  silent = TRUE)
        }


        incProgress(2/10, message = "Visualizations in progress")
        xgb.importance <- try(xgb.imp.list(xgb.list, level.names), silent = TRUE)

        xgb.loss <- list()
        xgb.imp.list <- list()
        xgb.shap <- list()
        xgb.shap.dep <- list()
        
        xgb_nom_length <- length(category.names(sam.dat.na, y.name))
        xgb_nom_name <- c()
        for(i in 1:xgb_nom_length){
          xgb_nom_name <- c(xgb_nom_name, input[[paste0("xgb_nom_label", i)]])
        }
        
        for(name in level.names){
          xgb.loss[[name]] <- try(xgb.error.plot.2(xgb.list, 
                                                   name), 
                                  silent = TRUE)
          xgb.shap[[name]] <- try(xgb.shap.summary.2(data = xgb.list[[name]]$data$x, 
                                                     top_n = as.numeric(input$xgb_nom_var_num), 
                                                     model = xgb.list[[name]]$model), 
                                  silent = TRUE)
          xgb.shap.dep[[name]] <- try(xgb.pdp.mult(xgb.list = xgb.list, 
                                                   n = as.numeric(input$xgb_nom_var_num), 
                                                   name = name, 
                                                   data.type = type, 
                                                   level.names = level.names,
                                                   label = xgb_nom_name), 
                                      silent = TRUE)
        }

        xgb.width.list.nom <- list()
        for(name in level.names){
          if(as.numeric(input$xgb_nom_var_num) > ncol(data[[name]])){
            n <- ncol(data[[name]])
            if(n %% 5 == 0){
              n <- n / 5
            }
            else{
              n <- (n / 5) + 1
            }
          }
          else{
            n <- as.numeric(input$xgb_nom_var_num)
            n <- n / 5
          }
          xgb.width.list.nom[[name]] <- n * 160
        }
        
        output$xgb_nom_results <- renderUI({
          tagList(
            do.call(tabsetPanel, lapply(1:length(level.names), function(i) {
              tabPanel(title = str_to_title(level.names[i]), align = "center",
                       tabsetPanel(tabPanel(title = "Importance", align = "center", br(),
                                            plotOutput(paste0("xgb_nom_SHAP", i), height = 700, width = 510), br(),
                                            dataTableOutput(paste0("xgb_nom_column_table1", i), height = "auto", width = 510)),
                                   tabPanel(title = "Partial Dependence", align = "center", br(),
                                            plotOutput(paste0("xgb_nom_SHAP_dep", i), height = 750, width = xgb.width.list.nom[[i]]), br(),
                                            dataTableOutput(paste0("xgb_nom_column_table2", i), height = "auto", width = 510)),
                                   tabPanel(title = "Error Plot", align = "center", br(),
                                            plotOutput(paste0("xgb_nom_loss", i), height = 700, width = 700))))}))
          )
        })

        lapply(1:length(level.names), function(j) {
          output[[paste0("xgb_nom_SHAP", j)]] <- renderPlot({
            tryCatch(xgb.shap[[level.names[j]]], error = function(e){
              message("Visualization not available! Check the input.")
              showModal(modalDialog(div("Visualization not available! Check the input.")))
              return(NULL)
            })
          })

          output[[paste0("xgb_nom_SHAP_dep", j)]] <- renderPlot({
            tryCatch(xgb.shap.dep[[level.names[j]]], error = function(e){
              message("Visualization not available! Check the input.")
              showModal(modalDialog(div("Visualization not available! Check the input.")))
              return(NULL)
            })
          })
          
          output[[paste0("xgb_nom_column_table1", j)]] <- renderDataTable({
            tryCatch(xgb.shap.imp.var(xgb.list, level.names[[j]], xgb.colnames.list.mult, n = as.numeric(input$xgb_nom_var_num)), error = function(e){
              message("Visualization not available! Check the input.")
              showModal(modalDialog(div("Visualization not available! Check the input.")))
              return(NULL)
            })
          })
          
          output[[paste0("xgb_nom_column_table2", j)]] <- renderDataTable({
            tryCatch(xgb.shap.imp.var(xgb.list, level.names[[j]], xgb.colnames.list.mult, n = as.numeric(input$xgb_nom_var_num)), error = function(e){
              message("Visualization not available! Check the input.")
              showModal(modalDialog(div("Visualization not available! Check the input.")))
              return(NULL)
            })
          })
          
          output[[paste0("xgb_nom_loss", j)]] <- renderPlot({
            tryCatch(xgb.loss[[level.names[j]]], error = function(e){
              message("Visualization not available! Check the input.")
              showModal(modalDialog(div("Visualization not available! Check the input.")))
              return(NULL)
            })
          })
          
          output$xgb_nom_downloadTabUI <- renderUI({
            tagList(
              p(" ", style = "margin-top: 20px;"),
              box(title = strong("Download Output Table", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  p("You can download the data analysis outputs.",
                    style = "font-size:11pt"),
                  h5("Feature Importance"),
                  downloadButton("xgb_nom_downloadTable1", "Download", width = '50%', style = "background-color: red3")
              )
            )
          })

          output$xgb_nom_downloadTable1 <- downloadHandler(
            filename = function() {
              paste("XGB_Importance.zip")
            },
            content = function(DA.file) {
              temp <- setwd(tempdir())
              on.exit(setwd(temp))
              if (length(level.names) == 5) {
                dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt")
                for(i in 1:length(dataFiles)){
                  write.table(as.data.frame(xgb.importance[[i]]), file = dataFiles[i], sep = "\t")
                }
              }
              else if(length(level.names) == 6) {
                dataFiles = c("Phylum.txt", "Class.txt", "Order.txt" ,"Family.txt", "Genus.txt", "Species.txt")
                for(i in 1:length(dataFiles)){
                  write.table(as.data.frame(xgb.importance[[i]]), file = dataFiles[i], sep = "\t")
                }
              }
              zip(zipfile=DA.file, files=dataFiles)
            })

          output$xgb_nom_reference <- renderUI({
            tagList(
              p(" ", style = "margin-top: 20px;"),
              box(title = strong("References", style = "color:white"), width = NULL, status = "info", solidHeader = TRUE,
                  p(ref, style = "font-size:11pt")
              )
            )
          })
        })
      }
    )

    shinyjs::enable("xgb_nom_dataType")
    shinyjs::enable("xgb_nom_response")
    shinyjs::enable("xgb_nom_data_input_opt")
    shinyjs::enable("xgb_nom_nfold")
    shinyjs::enable("xgb_nom_eta")
    shinyjs::enable("xgb_nom_penalty")
    shinyjs::enable("xgb_nom_nrounds")
    shinyjs::enable("xgb_nom_var_num")
    shinyjs::enable("xgb_nom_include_species")
    shinyjs::enable("rf_nom_runButton")
    shinyjs::enable("xgb_nom_runButton")
  })
}

# RUN ------
shinyApp(ui = ui, server = server)