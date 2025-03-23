library(readxl)
library(ggpubr)
library(ggplot2)
library(ggbeeswarm)

setwd("C:/Users/mevado/OneDrive - TUNI.fi/ISE/virus infection/papers/reovirus")
#setwd("D:/OneDrive/OneDrive - TUNI.fi/ISE/virus infection/papers/reovirus")

meta <- read_xlsx(path="./metadata_ileum.xlsx")
#abbr <- data.frame(`Cell.Type` = c("Enterocytes","Fibroblasts","Endothelial cells", "T cells", "Stem/TA",
#                                  "Macrophages", "Axial mural cells", "B cells", "Dendritic cells",
#                                  "Smooth-muscle cells", "Lymphatic cells", "Enteric neurons", "Goblet cells",
#                                  "Entero-endochrine cells"),
#                   `Abbr` = c("E","F", "Endo", "T","S/TA", "M",
#                             "AMC", "B", "DC", "SMC", "LC", "EN", "G", "EEC" ))
#meta <- merge(meta,abbr, by="Cell.Type")
#meta$Sample2 <- meta$Sample
#meta[meta$Sample %in% c("M1","M4"),]$Sample2 <- "M"
#writexl::write_xlsx(meta, path="./metadata_ileum.xlsx")

#counts <- read_xlsx(path="./Gene_Count_per_Cell.xlsx")
#counts <- read.table(file="./Gene_Count_per_Cell.txt", header=T)
ileum_infection_sc <- readRDS(file = "ileum_infection_sc.rds")

#cell percentage
pct_cells <- table(meta$Cell.Type,meta$Sample2)
pct_cells <- (t(pct_cells)/Matrix::colSums(pct_cells))*100
pct_cells
Matrix::rowSums(pct_cells)

gene <- "Hmgcs2"
TG2_counts <- as.data.frame(ileum_infection_sc@assays$RNA@data[paste0("GRCm38-",gene),]) #data- normalized, counts - raw counts
colnames(TG2_counts) <- gene
TG2_counts$Column <- rownames(TG2_counts)
TG2_counts <- merge(TG2_counts, meta, by = "Column")

p1 <- ggviolin(TG2_counts, x = "Abbr", y = gene,
         add = c("jitter"),facet.by = c("Sample2"), 
         order = c("E","S/TA","G", "EEC", "F", "Endo", "SMC","AMC", "T", "M", "B", "DC", "LC", "EN" ))

p1

TG2_counts2 <- TG2_counts[TG2_counts$Cell.Type %in% c("Enterocytes", "Stem/TA", "Goblet cells","Entero-endochrine cells"),]
p2 <- ggviolin(TG2_counts2, x = "Abbr", y = gene,
               add = c("jitter"),facet.by = c("Sample2"), 
               order = c("E","S/TA","G", "EEC"))

p2

ggexport(p2, filename = paste0(gene, "_expression_by_type.png"),
         width = 1800, length = 1800)
#-------------------------------------------------------------------Viral infection------------------------------------------------------

reoviral_genes <- c("ReoT1L-T1LReoS1", "ReoT1L-T1LReoS2", "ReoT1L-T1LReoS3", "ReoT1L-T1LReoS4", "ReoT1L-T1LReoM1", "ReoT1L-T1LReoM2", "ReoT1L-T1LReoM3", "ReoT1L-T1LReoL1", "ReoT1L-T1LReoL2", "ReoT1L-T1LReoL3")
viral_expression <- as.data.frame(ileum_infection_sc@assays$RNA@counts[reoviral_genes,])

vir <- as.data.frame(Matrix::colSums(viral_expression))
colnames(vir) <- "vir.counts"
vir$Column <- row.names(vir)
vir <- merge(vir,meta, by = "Column")

highlight <- vir[vir$vir.counts>0,]$Column
highlight_df <- TG2_counts[TG2_counts$Column %in% highlight,]

#p1 <- ggviolin(TG2_counts, x = "Abbr", y = gene,
#               add = c("jitter"),facet.by = c("Sample2"), order = c("E","S/TA","G", "EEC", "F", "Endo", "SMC","AMC", "T", "M", "B", "DC", "LC", "EN" ))+
#  geom_point(data=highlight_df, 
#             aes(x=Abbr,y=Muc1), 
#             color='red',
#             size=3)
#p1
#ggexport(p1, filename = "TGM2_expression_by_type+infection.png",
#         width = 1800, length = 1800)


vir <- merge(vir, TG2_counts[,c("Column",gene)], by="Column")
ggscatter(vir, x="vir.counts", y=gene)
#---------------------------------------------------------------------Stem cell population----------------------------------------------
meta <- read_xlsx(path="./Stem+TA_metada.xlsx")

counts <- read_xlsx(path="./Stem+TA_Gene_Count_per_Cell.xlsx")
row.names(counts) <- counts$Gene
counts[1:10,1:10]


TG2_counts <- as.numeric(t(counts["GRCm38-Tgm2",])[-647,])
correlation <- data.frame()
for(i in 1:nrow(counts)){
  y <- as.numeric(t(counts[counts$Gene[i],])[-647,])
  df <- data.frame(gene=counts$Gene[i], cor=cor(TG2_counts, y))
  correlation <- rbind(correlation,df)
}

y <- as.numeric(t(counts["GRCm38-Bola3",])[-647,])

plot(x=TG2_counts, y=y)

#Enterocytes progenitor
meta2 <- meta[meta$cell_type=="Enterocytes progenitor",]$Column
counts2 <- as.data.frame(counts[, c(meta2, "Gene")])

TG2_counts2 <- as.numeric(t(counts2[counts2$Gene == "GRCm38-Tgm2",])[-324,])
correlation2 <- data.frame()
for(i in 1:nrow(counts2)){
  y <- as.numeric(t(counts2[counts2$Gene == counts2$Gene[i],])[-324,])
  df <- data.frame(gene=counts2$Gene[i], cor=cor(TG2_counts2, y))
  correlation2 <- rbind(correlation2,df)
}

y <- as.numeric(t(counts2[counts2$Gene == "GRCm38-Mgat4b",])[-324,]) 

plot(x=TG2_counts2, y=y)


TG2_counts <- as.data.frame(TG2_counts)
TG2_counts$Column <- rownames(TG2_counts)
colnames(TG2_counts)[1] <- "Tgm2"
TG2_counts <- merge(TG2_counts, meta, by = "Column")
TG2_counts$Tgm2 <- as.numeric(TG2_counts$Tgm2)



p1 <- ggviolin(TG2_counts, x = "Abbr", y = "Tgm2",
               add = c("jitter", "mean_sd"),facet.by = c("Sample2"), order = c("S","P","TA","Epr"))
p1
ggexport(p1, filename = "TGM2_expression_by_Stem_type.png",
         width = 1800, length = 1800)

gene <- "Hmgcs2"
TG2_counts <- as.data.frame(t(counts[paste0("GRCm38-",gene),])) #data- normalized, counts - raw counts
colnames(TG2_counts) <- gene
TG2_counts$Column <- rownames(TG2_counts)
TG2_counts <- merge(TG2_counts, meta, by = "Column")
TG2_counts$Hmgcs2 <- as.numeric(TG2_counts$Hmgcs2)
p1 <- ggviolin(TG2_counts, x = "Abbr", y = gene,
               add = c("jitter", "mean_sd"),facet.by = c("label"), order = c("S","P","TA","Epr"))
p1

ggexport(p1, filename = paste0(gene, "_expression_by_type Stem cells subpopulation.png"),
         width = 1800, length = 1800)

#-------------------------------------------------------------------Viral infection------------------------------------------------------

reoviral_genes <- c("ReoT1L-T1LReoS1", "ReoT1L-T1LReoS2", "ReoT1L-T1LReoS3", "ReoT1L-T1LReoS4", "ReoT1L-T1LReoM1", "ReoT1L-T1LReoM2", "ReoT1L-T1LReoM3", "ReoT1L-T1LReoL1", "ReoT1L-T1LReoL2", "ReoT1L-T1LReoL3")

viral_expression <- counts[reoviral_genes,]

vir <- as.data.frame(Matrix::colSums(viral_expression[,1:646]))
colnames(vir) <- "vir.counts"
vir$Column <- row.names(vir)
vir <- merge(vir,meta, by = "Column")

highlight <- vir[vir$vir.counts>0,]$Column
highlight_df <- TG2_counts[TG2_counts$Column %in% highlight,]

p1 <- ggviolin(TG2_counts, x = "Abbr", y = "Tgm2",
               add = c("jitter", "mean_sd"),facet.by = c("Sample2"), order = c("S","P","TA","Epr"))+
  geom_point(data=highlight_df, 
             aes(x=Abbr,y=Tgm2), 
             color='red',
             size=3)
p1

#--------------------------------------------------------------------------MHCII------------------------------------------------------------
MHCII <- paste0("GRCm38-",c("H2-Ab1", "H2-Aa", "Ciita", "Cd74", "H2-DMa", "H2-DMb1"))
MHCII_expression <- counts[counts$Gene %in% MHCII,]

library(reshape2)
MHCII_expression_melt <- melt(MHCII_expression)
MHCII_expression_melt$log2 <- log2(MHCII_expression_melt$value+1)
MHCII_expression_melt <- merge(MHCII_expression_melt, meta, by.x = "variable", by.y = "Column")

p1 <- ggviolin(MHCII_expression_melt, x = "Abbr", y = "log2",
               add = c("jitter", "mean_sd"),facet.by = c("label"), order = c("S","P","TA","Epr"))
p1
ggexport(p1, filename = "MHCII_genes_expression_by_Stem_type.png",
         width = 1800, length = 1800)
