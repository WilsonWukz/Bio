# Install CRAN packages
install.packages("Seurat")
package.version("Seurat")
# install.packages("torch")
install.packages("sigclust")
install.packages("glmpca")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("scCATCH")
install.packages("dplyr")
install.packages("cowplot")
install.packages("patchwork")
install.packages("Matrix")
install.packages("RCurl")
install.packages("devtools")
install.packages("plotrix")
install.packages("magrittr")
install.packages("ggrepel")
install.packages("ggpubr")
install.packages("ggthemes")
BiocManager::install("multtest") #metap包需要这个包的依赖，而直接装会报错版本不和
install.packages("metap")
install.packages("plotly")
install.packages("readxl")
install.packages("ggnewscale")
install.packages("data.table")
install.packages("igraph")
install.packages("ggupset")
install.packages("europepmc")
install.packages("rgl")
install.packages("NMF")
install.packages("venn")
install.packages("ggalluvial")
install.packages("AnnoProbe")
install.packages("Cairo")


# Instal Bioconductor packages
BiocManager::install('tximport')
BiocManager::install('scran')
BiocManager::install('CoGAPS')
BiocManager::install('Nebulosa')
BiocManager::install('schex')
BiocManager::install('biomaRt', force = TRUE) #
BiocManager::install('AnnotationHub')
BiocManager::install('BiocGenerics')
BiocManager::install('SingleR')
BiocManager::install('celldex', force = TRUE) #
BiocManager::install('scRNAseq', force = TRUE) #
BiocManager::install('KEGGREST')
BiocManager::install('KEGGgraph')
BiocManager::install('mygene', force = TRUE) #
BiocManager::install('EnhancedVolcano')
BiocManager::install('DESeq2')
BiocManager::install('SingleCellExperiment',force = TRUE)
BiocManager::install('multtest', force = TRUE)
BiocManager::install('clusterProfiler', force = TRUE) #
BiocManager::install('DOSE', force = TRUE)
BiocManager::install('enrichplot', force = TRUE)
BiocManager::install('pathview', force = TRUE)
BiocManager::install("slingshot", force = TRUE)
BiocManager::install("monocle", force = TRUE)



options(timeout = 600)  # Increase to 600 seconds or more if needed
BiocManager::install('org.Hs.eg.db', force = TRUE)#Dependence for pathview
BiocManager::install('VennDetail', force = TRUE)

## load Github packages ()
library(usethis) #
library(gitcreds) #

usethis::use_git_config(user.name = "AccountName", user.email = "Email")
# usethis::create_github_token()

credentials::set_github_pat("YOUR_APT_KEY")
# Save your PAT in environ
usethis::edit_r_environ()

# Install Github packages
# 安装并加载 devtools 包
install.packages("devtools")
library(devtools)

# 使用 devtools::install_github() 从 GitHub 安装所需的 R 包
devtools::install_github("satijalab/seurat-wrappers")       # SeuratWrappers
devtools::install_github("perou-lab/MultiK")                   # MultiK
devtools::install_github("jinworks/CellChat")
BiocManager::install("ComplexHeatmap") #dependencies for cellchat

devtools::install_github("navinlabcode/copykat")            # copykat
devtools::install_github("cole-trapnell-lab/monocle3")      # monocle3
devtools::install_github("immunogenomics/harmony")          # harmony
devtools::install_github("JEFworks/liger")                  # liger
devtools::install_github("welch-lab/liger")                 # liger
devtools::install_github("dieterich-lab/CellPlot")          # CellPlot
devtools::install_github("r-lib/scales")                    # scales
install.packages("scales")                                  # github cant installed, then try CRAN

devtools::install_github("rstudio/reticulate")              # reticulate
devtools::install_github("immunogenomics/presto")           # presto
install.packages("magick")                                  #Dependence for scRNAtoolVis
devtools::install_github("junjunlab/jjAnno")                #Dependence for scRNAtoolVis
devtools::install_github("junjunlab/scRNAtoolVis")

devtools::install_github("TheHumphreysLab/plot1cell", force = T)
#Install dependencies for plot1cell
BiocManager::install('DoubletFinder', force = TRUE)
devtools::install_github("chris-mcginnis-ucsf/DoubletFinder")
BiocManager::install('loomR', force = TRUE)

devtools::install_github("cole-trapnell-lab/monocle3", force = TRUE)
devtools::install_github("mojaveazure/loomR")

BiocManager::install('EnsDb.Hsapiens.v86')
BiocManager::install('GEOquery')
BiocManager::install('simplifyEnrichment')

devtools::install_github("cole-trapnell-lab/monocle-release@develop")

#Install Dr.Meng's Package
devtools::install_github("LingzhangMeng/TidyGenePlot")
devtools::install_github("LingzhangMeng/scTidyGene")
devtools::install_github("LingzhangMeng/OptiRes")
# Install My Github packages
devtools::install_github("WilsonWukz/EasyCellMarker2")
devtools::install_github("WilsonWukz/SimilarityDistanceCalculate")
