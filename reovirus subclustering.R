library(readxl)
library(ggpubr)
library(Seurat)

setwd("C:/Users/mevado/OneDrive - TUNI.fi/ISE/virus infection/papers/reovirus")
#setwd("D:/OneDrive/OneDrive - TUNI.fi/ISE/virus infection/papers/reovirus")

meta <- read_xlsx(path="./metadata_ileum.xlsx")

stem_cells <- meta[meta$Cell.Type == "Stem/TA",]
pct_stem <- table(stem_cells$Sample2)/table(meta$Sample2)*100

ileum_infection_sc <- readRDS(file = "ileum_infection_sc.rds")
head(ileum_infection_sc@meta.data)

filtered_matrix <- ileum_infection_sc@assays$RNA@counts[,stem_cells$Column2]
dim(filtered_matrix) #[1] 18885   646
ileum_infection_STA <- CreateSeuratObject(counts = filtered_matrix)
ileum_infection_STA@meta.data$sub.cluster


names(ileum_infection_sc@graphs)
#Subclustering
ileum_infection_STA <- FindSubCluster(ileum_infection_sc, "Stem/TA",  graph.name="RNA_snn", resolution = 0.5,   algorithm = 1 )

head(ileum_infection_STA@meta.data)
ileum_infection_STA <- subset(x = ileum_infection_STA, idents = "Stem/TA")
saveRDS(ileum_infection_STA, file = "ileum_infection_Stem+TA.rds")
DimPlot(ileum_infection_STA, reduction = "umap",group.by = "sub.cluster", label = TRUE, label.size = 6)
ileum_infection_STA <- SetIdent(ileum_infection_STA, value = ileum_infection_STA@meta.data$sub.cluster)

seurat.markers <- FindAllMarkers(ileum_infection_STA,  test.use = "wilcox", logfc.threshold = 0.25)
top10 <- seurat.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
writexl::write_xlsx(seurat.markers, path ="Markers_Stem+TA_subclusters_ileum.xlsx")
writexl::write_xlsx(top10, path ="Markers_Stem+TA_subclusters_ileum_top10.xlsx")

ileum_infection_STA <- readRDS(file = "ileum_infection_Stem+TA.rds")
head(ileum_infection_STA@meta.data)
ileum_infection_STA <- SetIdent(ileum_infection_STA, value = ileum_infection_STA@meta.data$sub.cluster)
head(ileum_infection_STA@meta.data)

new.cluster.ids <- c("TA",
                     "Stem",
                     "Enterocytes progenitor",
                     "Enterocytes progenitor",
                     "Paneth",
                     "Stem")
names(new.cluster.ids) <- levels(ileum_infection_STA)
ileum_infection_STA <- RenameIdents(ileum_infection_STA, new.cluster.ids)
UMAP <- DimPlot(ileum_infection_STA, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
UMAP
ggexport(UMAP, filename ="UMAP_Stem+TA.png")

p2 <- FeaturePlot(ileum_infection_STA, features = paste0("GRCm38-",c("Tgm2")))
p2
ggexport(p2, filename ="TGM2_ileum_Stem+TA.png",
         height = 2000, width = 2000, res = 300)

metadata <- as.data.frame(ileum_infection_STA@active.ident)
metadata$cell <- rownames(metadata)
writexl::write_xlsx(metadata, path = 'Stem+TA_metada.xlsx')
counts <- as.data.frame(ileum_infection_STA@assays[["RNA"]]@counts)
counts$Gene <- rownames(counts)
writexl::write_xlsx(counts, path = 'Stem+TA_Gene_Count_per_Cell.xlsx')
   