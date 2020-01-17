###
### Expression Analysis, Configuration Commands
###


## libraries


## core utilities

# install.packages("BiocManager")
# install.packages("devtools")
# install.packages("https://cran.r-project.org/src/contrib/Archive/curl/curl_4.0.tar.gz",repo=NULL,type="source")
BiocManager::install()

## plotting

# install.packages("ggplot2")
# install.packages("ggfortify")
# install.packages("pheatmap")
# install.packages("RColorBrewer")
# install.packages("Rtsne")
# install.packages("pROC")
# install.packages("VennDiagram")
# install.packages("cowplot")
# install.packages("gridGraphics")
library("ggplot2")
# library("ggfortify")
library("pheatmap")
library("RColorBrewer")
# library("Rtsne")
# library("pROC")
# library("VennDiagram")
library("cowplot")
library("gridGraphics")

## data parsing

# BiocManager::install("statmod")
# BiocManager::install("apeglm")
# BiocManager::install("BiocParallel")
# install.packages("tidyverse")
# install.packages("dplyr")
# install.packages("readr")
# install.packages("readxl")
library("statmod")
library("apeglm")
library("BiocParallel")
# library("tidyverse")
library("readr")
# library("dplyr")
library("readxl")

## annotation

# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("ensembldb")
# BiocManager::install("EnsDb.Hsapiens.v86")
# BiocManager::install("EnsDb.Mmusculus.v79")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("ensembldb")
library("EnsDb.Hsapiens.v86")

## gene/txp expression

# BiocManager::install("EDASeq")
# BiocManager::install("limma")
# BiocManager::install("edgeR")
# BiocManager::install("DESeq2")
# BiocManager::install("sva")
# BiocManager::install("RUVSeq")
# BiocManager::install("WGCNA")
# BiocManager::install("tximport")
# library("EDASeq")
# library("limma")
# library("edgeR")
library("DESeq2")
library("sva")
# library("RUVSeq")
# library("WGCNA"); allowWGCNAThreads()
# library("tximport")

## gene set enrichment

# BiocManager::install("goseq")
# BiocManager::install("GO.db")
# BiocManager::install("fgsea")
# BiocManager::install("reactome.db")
# BiocManager::install("GSEABase")
# BiocManager::install("GSVA")
# devtools::install_github('dviraran/xCell')
library("goseq")
library("GO.db")
library("fgsea")
# library("reactome.db")
library("GSEABase")            # must be loaded for xCell
library("GSVA")                # must be loaded for xCell
library("xCell")

## survival

# install.packages("survival")
# library("survival")

## data

# BiocManager::install("TCGAbiolinks") #for 3.5.1 required archived version of cmprsk `install_version("cmprsk", version="2.2-7", repos="http://cran.us.r-project.org")`
# library('TCGAbiolinks')
# BiocManager::install("curatedCRCData")
# library("curatedCRCData")

## methylation
# BiocManager::install("ChAMP")
# library("ChAMP")


## functions

set_group_colors <- function(annotation_file_name, groups_of_interest) {
  ## define group colors
  setwd(anno_dir)
  pheno <- read.table(file=annotation_file_name, sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="--")
  colors <- list()
  for (i in 1:length(groups_of_interest)) {
    if (class(pheno[, groups_of_interest[i]]) == "character") {
      gn <- length(levels(as.factor(pheno[, groups_of_interest[i]])))
      ## get hcl color wheel coordinates for various cohorts (ggplot starts at 15 and increases to 360 by 360/gn)
      inc <- 360/gn
      cords <- vector()
      for (j in seq(15, 360, inc)) {
        if (j==15) {
          cords <- j
        } else {
          cords <- append(cords, j)
        }
      }
      ## get hcl hexidecimal color codes
      hex <- hcl(h=cords, c=100, l=65)
      names(hex) <- levels(as.factor(pheno[,groups_of_interest[i]]))
      colors[[i]] <- hex
    } else {
      # continuous_values <- pheno[, groups_of_interest[i]]
      # continuous_values <- continuous_values[!is.na(continuous_values)]
      colors[[i]] <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
      # colors[[i]] <- continuous_values
    }
  }
  names(colors) <- groups_of_interest
  return(colors)
}


## project variables

## datasets
base_sets <- c("simple","mss","msi_mss_upper","msi_mss_lower")
names(base_sets) <- base_sets

## paths
anno_dir <- "/scratch/chd5n/aneuploidy/raw-data/annotations/"
count_dir <- "/scratch/chd5n/aneuploidy/raw-data/counts/"
r_dir <- "/scratch/chd5n/aneuploidy/results/rdata/"
plot_dir <- "/scratch/chd5n/aneuploidy/results/plots/"
table_dir <- "/scratch/chd5n/aneuploidy/results/tables/"
script_dir <- "/home/chd5n/projects/aneuploidy/scripts/"
#wgcna_dir <- paste0(work_dir, "wgcna/")
#enrich_dir <- paste0(work_dir, "enrich/")

## colors
filename <- paste0(anno_dir, paste0("rna_set_",set_name,".tsv"))
groups <- c("category","disease_type","primary_site","tissue_type","stage","msi_status","vital_status","race","sex","capture_kit_name")
names(groups) <- groups
group_colors <- set_group_colors(filename,groups)

cat("Aneuploidy configuration loaded. Ready to analyze!\n\n")
