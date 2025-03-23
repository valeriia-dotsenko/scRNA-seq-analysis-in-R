# Load libraries
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(Seurat)

setwd("C:/Users/mevado/OneDrive - TUNI.fi/ISE/virus infection/papers/reovirus")

meta <- read_xlsx(path="./Stem+TA_metada.xlsx")

# Read in the dataset
sce <- readRDS(file = "ileum_infection_Stem+TA.rds")
Idents(sce)
head(sce@meta.data)

# run all DE methods
methods <- c("wilcox")
cells.1	<- meta[  meta$label == "M",]$Column #meta$cell_type=="Stem"&
cells.2 <- meta[meta$label == "V1",]$Column #meta$cell_type=="Stem" & 

DE <- list()
for (m in methods){
  outfile <- paste("./",m,"_M_vs_V1.tab", sep='')
  if(!file.exists(outfile)){
    DE[[m]]<- FindMarkers(object = sce, ident.1 = cells.1, ident.2 = cells.2,test.use = m)
    write.table(DE[[m]],file=outfile,sep="\t",quote=F)
  }
}


methods <- c("wilcox")
cells.1	<- meta[meta$cell_type=="Stem"& meta$label == "M",]$Column #meta$cell_type=="Stem"&
cells.2 <- meta[meta$cell_type=="Stem"& meta$label == "V1",]$Column #meta$cell_type=="Stem" & 

DE <- list()
for (m in methods){
  outfile <- paste("./",m,"_M_vs_V1_Stem cells only.tab", sep='')
  if(!file.exists(outfile)){
    DE[[m]]<- FindMarkers(object = sce, ident.1 = cells.1, ident.2 = cells.2,test.use = m)
    write.table(DE[[m]],file=outfile,sep="\t",quote=F)
  }
}
