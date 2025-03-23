# scRNA-seq-analysis-in-R
This is analysis of scRNA seq GSE189636 dataset in R.

Data is published in paper: 	
> Mantri M, Hinchman MM, McKellar DW, Wang MFZ et al. Spatiotemporal transcriptomics reveals pathogenesis of viral myocarditis. Nat Cardiovasc Res 2022 Oct;1(10):946-960. PMID: 36970396

The initial analysis by authors performed in Phyton - https://github.com/madhavmantri/reovirus_induced_myocarditis/blob/main/notebooks/ileum_scRNAseq_processing.ipynb

## reovirus data reproduce in R.R
I reproduced ileum analysis in R:

UMAP plot:

![UMAP_ileum](https://github.com/user-attachments/assets/0e55534d-c05c-4eb4-b3c2-4b5704fb42ff)

FeaturePlot of cell type gene markers expression:

![Markers2_ileum](https://github.com/user-attachments/assets/2ad75ea1-d7dd-4881-8ee1-4672612d5a95)

UMAP + FeaturePlot by sample Types: 

![UMAP + Hmgcs2 by sample](https://github.com/user-attachments/assets/dd1e88e6-1b21-4f70-b2f9-d6d98f7ba64e)

## reovirus subclustering file.R
Stem/TA cluster used for subsclustering into "TA", "Stem", "Enterocytes progenitor", and "Paneth" cell types. 

UMAP plot for subclustered cells:

![UMAP_Stem+TA](https://github.com/user-attachments/assets/10366beb-a368-4fb6-b772-d872346885ad)

## Checking gene expression in reovirus ileum.R
Tgm2 expression by cell and sample type:

![TGM2_expression_by_type001](https://github.com/user-attachments/assets/6096bd54-7230-449e-bddc-a18af83b7b72)

## Differential expression analysis for scRNAseq.R
Wilcox test for DEG expression. Head of the results for M vs V1 comparisons when Stem cells compartment is subsetted. 

|         | p_val      | avg_log2FC     | pct.1 | pct.2 | p_val_adj   |
|----------------|------------|----------------|-------|-------|-------------|
| GRCm38-Bst2    | 1.48E-09   | -3.673358687   | 0.312 | 1     | 2.79E-05    |
| GRCm38-Isg15   | 1.44E-09   | -3.667286995   | 0.062 | 1     | 2.72E-05    |
| GRCm38-Ifitm3  | 5.07E-08   | -3.624499485   | 0     | 0.889 | 0.000957401 |
| GRCm38-Xist    | 2.43E-09   | -3.573824539   | 0     | 0.981 | 4.59E-05    |
| GRCm38-Ly6e    | 2.65E-09   | -3.066591269   | 0.062 | 0.981 | 5.01E-05    |



