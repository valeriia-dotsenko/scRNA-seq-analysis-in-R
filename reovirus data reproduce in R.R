library(Seurat)
library(hdf5r) #devtools::install_github("hhoeflin/hdf5r")
library(readxl)
library(dplyr)
library(scclusteval) #devtools::install_github("crazyhottommy/scclusteval")
library(ggpubr)

setwd("C:/Users/mevado/OneDrive - TUNI.fi/ISE/virus infection/papers/reovirus")
#setwd("~/Documents/reovirus")
Mock_Ileum_D1PI_sc <-  Read10X_h5("./GSM5705262_scRNAseq_Mock_Ileum_D1PI_filtered_feature_bc_matrix.h5")
T1L_Ileum_D1PI_sc <-  Read10X_h5("./GSM5705263_scRNAseq_T1L_Ileum_D1PI_filtered_feature_bc_matrix.h5")
Mock_Ileum_D4PI_sc <-  Read10X_h5('./GSM5705265_scRNAseq_Mock_Ileum_D4PI_filtered_feature_bc_matrix.h5')
T1L_Ileum_D4PI_sc <-  Read10X_h5('./GSM5705266_scRNAseq_T1L_Ileum_D4PI_filtered_feature_bc_matrix.h5')


celseq1 <- CreateSeuratObject(counts = Mock_Ileum_D1PI_sc, project = "M1")
celseq2 <- CreateSeuratObject(counts = T1L_Ileum_D1PI_sc, project = "V1")
celseq3 <- CreateSeuratObject(counts = Mock_Ileum_D4PI_sc, project = "M4")
celseq4 <- CreateSeuratObject(counts = T1L_Ileum_D4PI_sc, project = "V4")

ileum_infection_sc <- merge(celseq1, y = c(celseq2, celseq3,celseq4), 
                            add.cell.ids = c("M1", "V1", "M4","V4"), project = "ileum_infection_sc")

dim(ileum_infection_sc) #[1] 54848 24901
ileum_infection_sc
unique(sapply(X = strsplit(colnames(ileum_infection_sc), split = "_"), FUN = "[", 1))
table(ileum_infection_sc$orig.ident)

#saveRDS(ileum_infection_sc, file = "ileum_infection_sc.rds")

ileum_infection_sc <- readRDS(file = "ileum_infection_sc.rds")
reoviral_genes <- c("ReoT1L-T1LReoS1", "ReoT1L-T1LReoS2", "ReoT1L-T1LReoS3", "ReoT1L-T1LReoS4", "ReoT1L-T1LReoM1", "ReoT1L-T1LReoM2", "ReoT1L-T1LReoM3", "ReoT1L-T1LReoL1", "ReoT1L-T1LReoL2", "ReoT1L-T1LReoL3")

viral_expression <- ileum_infection_sc@assays$RNA@counts[reoviral_genes,]
dim(viral_expression) #[1]    10 24901

viral_expression.columns = dimnames(viral_expression)[2]

# filter cells and genes
counts_per_cell <- Matrix::colSums(ileum_infection_sc@assays$RNA@counts)
counts_per_gene <- Matrix::rowSums(ileum_infection_sc@assays$RNA@counts)

a <- counts_per_cell[counts_per_cell > 1] #154, not 178
b <- counts_per_gene[counts_per_gene != 0] #25313

filtered_matrix <- ileum_infection_sc@assays$RNA@counts[names(b),names(a)]
dim(filtered_matrix) #[1] 29535 24747 
ileum_infection_sc <- CreateSeuratObject(counts = filtered_matrix)
head(ileum_infection_sc@meta.data)

#mt genes percent
mt_genes <- grep(pattern = "^GRCm38-mt-", x = rownames(x = ileum_infection_sc), value = TRUE)
pct_mito <- (Matrix::colSums(ileum_infection_sc[mt_genes, ]))/(Matrix::colSums(ileum_infection_sc))*100
#pct_mito[1:5]
head(ileum_infection_sc@meta.data)  # Before adding
ileum_infection_sc <- AddMetaData(object = ileum_infection_sc, metadata = pct_mito, col.name = "percent.mito")
head(ileum_infection_sc@meta.data)  # After adding

#viral genes total + percent
total_viral  <-  Matrix::colSums(viral_expression)
#total_viral[total_viral > 10]
ileum_infection_sc <- AddMetaData(object = ileum_infection_sc, metadata = total_viral, col.name = "total_viral")
head(ileum_infection_sc@meta.data)  # After adding

pct_viral = (Matrix::colSums(ileum_infection_sc[reoviral_genes, ]))/(Matrix::colSums(ileum_infection_sc))*100
#pct_viral[pct_viral > 1]
ileum_infection_sc <- AddMetaData(object = ileum_infection_sc, metadata = pct_viral, col.name = "pct_viral")
head(ileum_infection_sc@meta.data)  # After adding

#violin plot
qc <- VlnPlot(object = ileum_infection_sc, features =c('nCount_RNA',"nFeature_RNA","percent.mito","pct_viral"))
#qc
ggexport(qc, filename ="violin_ileum_qc.png")

# filter data 
ileum_infection_sc <- subset(ileum_infection_sc, subset = percent.mito < 30)
#summary(ileum_infection_sc@meta.data[,"percent.mito"])
dim(ileum_infection_sc) #[1] 29535  8597

ileum_infection_sc <- subset(ileum_infection_sc, subset = nCount_RNA >= 200)
dim(ileum_infection_sc) #[1] 29535  8083  514 cells filtered

counts_per_gene <- Matrix::rowSums(ileum_infection_sc@assays$RNA@counts)
b <- counts_per_gene[counts_per_gene > 10] #10326

filtered_matrix <- ileum_infection_sc@assays$RNA@counts[names(b),]
dim(filtered_matrix) #[1] 18885  8083
ileum_infection_sc <- CreateSeuratObject(counts = filtered_matrix)
head(ileum_infection_sc@meta.data)
table(ileum_infection_sc$orig.ident)

#Normalize data
ileum_infection_sc <- NormalizeData(ileum_infection_sc, normalization.method = "LogNormalize", scale.factor = 1e4)

#Identification of highly variable features (feature selection)
ileum_infection_sc <- FindVariableFeatures(ileum_infection_sc)
plot1 <- VariableFeaturePlot(ileum_infection_sc)
plot1


#Cell cycle
s_genes <- c('Mcm4', 'Exo1', 'Slbp', 'Gmnn', 'Cdc45', 'Msh2', 'Mcm6', 'Rrm2', 'Pold3', 'Blm', 'Ubr7', 'Mcm5',
             'Clspn', 'Hells', 'Nasp', 'Rpa2', 'Rad51ap1', 'Tyms', 'Rrm1', 'Rfc2', 'Prim1', 'Brip1', 'Usp1', 
             'Ung', 'Pola1', 'Mcm2', 'Fen1', 'Tipin', 'Pcna', 'Cdca7', 'Uhrf1', 'Casp8ap2', 'Cdc6', 'Dscc1', 
             'Wdr76', 'E2f8', 'Dtl', 'Ccne2', 'Atad2', 'Gins2', 'Chaf1b', 'Pcna-ps2')
s_genes <- paste0("GRCm38-",s_genes)
g2m_genes = c('Nuf2', 'Psrc1', 'Ncapd2', 'Ccnb2', 'Smc4', 'Lbr', 'Tacc3', 'Cenpa', 'Kif23', 'Cdca2', 'Anp32e',
              'G2e3', 'Cdca3', 'Anln', 'Cenpe', 'Gas2l3', 'Tubb4b', 'Cenpf', 'Dlgap5', 'Hjurp',  'Gtse1',
              'Bub1', 'Birc5', 'Ube2c', 'Rangap1', 'Hmmr', 'Ect2', 'Tpx2', 'Ckap5', 'Cbx5', 'Nek2', 'Ttk', 'Cdca8',
              'Nusap1', 'Ctcf', 'Cdc20', 'Cks2', 'Mki67', 'Tmpo', 'Ckap2l', 'Aurkb', 'Kif2c', 'Cdk1', 'Kif20b', 
              'Top2a', 'Aurka', 'Ckap2', 'Hmgb2', 'Cdc25c', 'Ndc80', 'Kif11') #'Cks1brt',
g2m_genes <- paste0("GRCm38-",g2m_genes)

ileum_infection_sc <- CellCycleScoring(ileum_infection_sc, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)
# view cell cycle scores and phase assignments
head(ileum_infection_sc@meta.data)
# Visualize the distribution of cell cycle markers across
RidgePlot(ileum_infection_sc, features = c("GRCm38-Pcna", "GRCm38-Top2A", "GRCm38-Mcm6", "GRCm38-Mki67"), ncol = 2)


#Scaling the data
mt_genes <- grep(pattern = "^GRCm38-mt-", x = rownames(x = ileum_infection_sc), value = TRUE)
pct_mito <- (Matrix::colSums(ileum_infection_sc[mt_genes, ]))/(Matrix::colSums(ileum_infection_sc))*100
#pct_mito[1:5]
head(ileum_infection_sc@meta.data)  # Before adding
ileum_infection_sc <- AddMetaData(object = ileum_infection_sc, metadata = pct_mito, col.name = "percent.mito")
head(ileum_infection_sc@meta.data)  # After adding

all.genes <- rownames(ileum_infection_sc)
ileum_infection_sc <- ScaleData(ileum_infection_sc,vars.to.regress = c("nCount_RNA", "percent.mito", "S.Score", "G2M.Score"), features = all.genes)

saveRDS(ileum_infection_sc, file = "ileum_infection_sc.rds")
ileum_infection_sc <- readRDS(file = "ileum_infection_sc.rds")
head(ileum_infection_sc@meta.data)  # After adding

#Perform linear dimensional reduction
ileum_infection_sc <- RunPCA(ileum_infection_sc, features = VariableFeatures(object = ileum_infection_sc))
# Examine and visualize PCA results a few different ways
DimPlot(ileum_infection_sc, reduction = "pca", group.by = c("orig.ident","nCount_RNA","nFeature_RNA","percent.mito","Phase"))

#Cluster the cells
ileum_infection_sc <- FindNeighbors(ileum_infection_sc, dims = 1:20)
ileum_infection_sc <- RunUMAP(ileum_infection_sc, dims = 1:20)
# note that you can set `label = TRUE` or use the LabelClusters function to help label  individual clusters
DimPlot(ileum_infection_sc, reduction = "umap", group.by = "orig.ident")

#Cluster cells
ileum_infection_sc <- FindClusters(ileum_infection_sc, resolution = 0.5)
head(Idents(ileum_infection_sc), 5)

ileum_infection_sc <- FindClusters(ileum_infection_sc, resolution = 0.3)

head(ileum_infection_sc@meta.data)
head(Idents(ileum_infection_sc), 5)

DimPlot(ileum_infection_sc, reduction = "umap", group.by = "RNA_snn_res.0.5")
DimPlot(ileum_infection_sc, reduction = "umap", group.by = "RNA_snn_res.0.3")

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
seurat.markers <- FindAllMarkers(ileum_infection_sc,  test.use = "wilcox", logfc.threshold = 0.25)
seurat.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

new.cluster.ids <- c("Enterocytes",
                     "Fibroblasts",
                     "Endothelial cells",
                     "T cells",
                     "Stem/TA",
                     "Enterocytes",
                     "Macrophages",
                     "Fibroblasts",
                     "Axial mural cells",
                     "B cells",
                     "Dendritic cells",
                     "Fibroblasts",
                     "Smooth-muscle cells",
                     "Lymphatic cells",
                     "Enteric neurons",
                     "Goblet cells",
                     "Entero-endochrine cells")
names(new.cluster.ids) <- levels(ileum_infection_sc)
ileum_infection_sc <- RenameIdents(ileum_infection_sc, new.cluster.ids)
UMAP <- DimPlot(ileum_infection_sc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
UMAP
ggexport(UMAP, filename ="UMAP_ileum.png")
head(ileum_infection_sc@meta.data)

saveRDS(ileum_infection_sc, file = "ileum_infection_sc.rds")


markers <- VlnPlot(ileum_infection_sc, features = paste0("GRCm38-",c("Lyz1", "Olfm4", "Apoa4","Mki67",'Muc2',"Chga","Dclk1","Tgm2")),slot = "counts", log = TRUE)
ggexport(markers, filename ="Markers_ileum.png")
markers <- VlnPlot(ileum_infection_sc, features = paste0("GRCm38-",c("Tgm2")),slot = "counts", log = TRUE)
ggexport(markers, filename ="Markers_ileum.png")

p1 <- FeaturePlot(ileum_infection_sc, features = paste0("GRCm38-",c("Lyz1", "Olfm4", "Apoa4","Mki67",'Muc2',"Chga","Dclk1","Tgm2")))
p1
ggexport(p1, filename ="Markers2_ileum.png",
         height = 2000, width = 2000, res = 300)

p2 <- FeaturePlot(ileum_infection_sc, features = paste0("GRCm38-",c("Hmgcs2")), split.by = "orig.ident")
p2
ggexport(p2, filename ="TGM2_ileum.png",
         height = 2000, width = 2000, res = 300)


#install.packages("clustermole")
#BiocManager::install("GSVA")
#BiocManager::install("singscore")

#library(clustermole)
levels(Idents(ileum_infection_sc))
p2 <- FeaturePlot(ileum_infection_sc, features = "percent.mito")
p2
ggexport(p2, filename ="TGM2_ileum.png",
         height = 2000, width = 2000, res = 300)

write.table(ileum_infection_sc@active.ident, file='Convert_UMI_Label.tsv', quote=FALSE, sep='\t', col.names = TRUE)

write.table(ileum_infection_sc@assays[["RNA"]]@counts, file = 'Gene_Count_per_Cell.tsv')


# calculate each cell type numbers
df <- ileum_infection_sc@meta.data
new.cluster.ids <- as.data.frame(c("0"="Enterocytes",
                                   "1"="Fibroblasts",
                                   "2"="Endothelial cells",
                                   "3"="T cells",
                                   "4"="Stem/TA",
                                   "5"="Enterocytes",
                                   "6"="Macrophages",
                                   "7"="Fibroblasts",
                                   "8"="Axial mural cells",
                                   "9"="B cells",
                                   "10"="Dendritic cells",
                                   "11"="Fibroblasts",
                                   "12"="Smooth-muscle cells",
                                   "13"="Lymphatic cells",
                                   "14"="Enteric neurons",
                                   "15"="Goblet cells",
                                   "16"="Entero-endochrine cells"))
new.cluster.ids$seurat_clusters <- row.names(new.cluster.ids)
colnames(new.cluster.ids) <- c( "cluster_name","seurat_clusters")

df <- merge(df, new.cluster.ids, by = "seurat_clusters")

df %>%
  group_by(orig.ident, cluster_name) %>%
  summarise(n = n()) -> cellnumber

ggbarplot(cellnumber[cellnumber$cluster_name %in% "Stem/TA" & cellnumber$orig.ident != "M1",], x=colnames(cellnumber)[1], y=colnames(cellnumber)[3],label = TRUE, label.pos = "out", fill = colnames(cellnumber)[1])+
  scale_x_discrete(labels = c("M","Vd1", "Vd4"))+
  ylab("Number Stem/TA cells")+
  scale_fill_manual(values=c("#00AFBB", "#E7B800", "#FC4E07"))+
  ggtitle("Reovirus infected mouse ileum")+
  theme( axis.title.x = element_blank(),
         legend.position = "none")

#some additional plots
# Rename identity classes, so M1 and M4 are the same

orig.ident <- as.vector(ileum_infection_sc$orig.ident)

orig.ident <- replace(orig.ident, orig.ident == "M1", "M")
orig.ident <- replace(orig.ident, orig.ident == "M4", "M")

ileum_infection_sc$orig.ident <- NULL
head(ileum_infection_sc@meta.data)

ileum_infection_sc <- AddMetaData(
  object = ileum_infection_sc,
  metadata = orig.ident,
  col.name = 'study.rep'
)

UMAP <- DimPlot(ileum_infection_sc, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "study.rep" , ncol = 3) + NoLegend()
UMAP

p2 <- FeaturePlot(ileum_infection_sc, features = paste0("GRCm38-",c("Hmgcs2")), split.by = "study.rep")
p2

p1_3 <- ggarrange(UMAP,p2,nrow = 2,heights = c(1,1))
p1_3
ggexport(p1_3, filename = "UMAP + Hmgcs2 by sample.tiff",
         width = 1200, # 7 inch
         height = 700)
