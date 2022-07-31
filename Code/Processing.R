library(Seurat)
library(biomaRt)
library(dplyr)
library(ggplot2)

#### gene homology (modified from https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/) ####

# Basic function to convert mouse to human gene names

Seu_Mouse_data <- readRDS("/external/rprshnas01/kcni/ychen/git/Ex_Env_Storage/Nr1_celltypes/DGE - Seurat AIBS Object.rds")
Seu_Human_data <- readRDS("/external/rprshnas01/netdata_kcni/stlab/Public/Seurat_objects/Seu_AIBS_obj_update_07JUN21.rds")

testgenes <- row.names(Seu_Mouse_data)
testgenes <- row.names(Seu_Human_data)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

MGI_HGNC = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = testgenes, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
HGNC_MGI = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = testgenes , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)

# dup check and append
dup_hum = MGI_HGNC[duplicated(MGI_HGNC$HGNC.symbol), "HGNC.symbol"]
MGI_HGNC$Duplicate_HGNC <- MGI_HGNC$HGNC.symbol %in% dup_hum
dup_mou = MGI_HGNC[duplicated(MGI_HGNC$MGI.symbol), "MGI.symbol"]
MGI_HGNC$Duplicate_MGI <- MGI_HGNC$MGI.symbol %in% dup_mou
MGI_HGNC$One_to_one <- MGI_HGNC$Duplicate_HGNC == F & MGI_HGNC$Duplicate_MGI == F
write.csv(MGI_HGNC, "./MGI_HGNC.csv")

dup_hum = HGNC_MGI[duplicated(HGNC_MGI$HGNC.symbol), "HGNC.symbol"]
HGNC_MGI$Duplicate_HGNC <- HGNC_MGI$HGNC.symbol %in% dup_hum
dup_mou = HGNC_MGI[duplicated(HGNC_MGI$MGI.symbol), "MGI.symbol"]
HGNC_MGI$Duplicate_MGI <- HGNC_MGI$MGI.symbol %in% dup_mou
HGNC_MGI$One_to_one <- HGNC_MGI$Duplicate_HGNC == F & HGNC_MGI$Duplicate_MGI == F
write.csv(HGNC_MGI, "./HGNC_MGI.csv")

# consolidate

#MGI_HGNC <- read.csv("./MGI_HGNC.csv", row.names = 1)
#HGNC_MGI <- read.csv("./HGNC_MGI.csv", row.names = 1)

HGNC_MGI <- HGNC_MGI[,c(2,1,3:5)]
colnames(MGI_HGNC) == colnames(HGNC_MGI)
Homologous_genes_mouhum <- rbind(HGNC_MGI, MGI_HGNC)
Homologous_genes_mouhum <- unique(Homologous_genes_mouhum)
write.csv(Homologous_genes_mouhum, "./Homologous_genes_mouhum.csv")

# filter

Hom_genes_filtered = Homologous_genes_mouhum[Homologous_genes_mouhum$One_to_one == T,] 

hum_counts <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_human_counts_mini.csv", row.names = 1)
mou_counts <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_mouse_counts_mini.csv", row.names = 1)

Hom_genes_filtered <- Hom_genes_filtered[Hom_genes_filtered$MGI.symbol %in% colnames(mou_counts),]
Hom_genes_filtered <- Hom_genes_filtered[Hom_genes_filtered$HGNC.symbol %in% colnames(hum_counts),]

write.csv(Hom_genes_filtered, "./Homologous_genes_mouhum_filtered.csv")
Hom_genes_filtered <- read.csv("./Homologous_genes_mouhum_filtered.csv", row.names = 1)

#
#### Make Seurat objects ####

# human
hum_counts <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_human_counts_mini.csv", row.names = 1)
hum_meta <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_human_meta_mini.csv", row.names = 1)
row.names(hum_meta) <- hum_meta$sample_name
Seu_hum <- CreateSeuratObject(counts = t(hum_counts), meta.data = hum_meta)
Seu_hum <- readRDS("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/Seu_hum_mini.rds")

# mouse (with gene name conversion for homologous genes)
mou_counts <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_mouse_counts_mini.csv", row.names = 1)
mou_meta <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_mouse_meta_mini.csv", row.names = 1)
row.names(mou_meta) <- mou_meta$sample_name

Hom_genes_filtered <- read.csv("./Homologous_genes_mouhum_filtered.csv", row.names = 1)
test <- merge(t(mou_counts), Hom_genes_filtered[,1:2], by.x = "row.names", by.y = "MGI.symbol", all.x = T, all.y = F)
test <- test %>% mutate(HGNC.symbol = coalesce(HGNC.symbol, Row.names))
row.names(test) <- test$HGNC.symbol
test <- test[,-1]
test <- test[,-2870]
mou_counts <- test

Seu_mou <- CreateSeuratObject(counts = mou_counts, meta.data = mou_meta)
Seu_mou <- readRDS("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/Seu_mou_mini.rds")

#
#### Integration + initial analysis ####

Seu.list <- c(Seu_hum, Seu_mou)

# normalize and identify variable features for each dataset independently
Seu.list <- lapply(X = Seu.list, FUN = function(x) {
  #x <- NormalizeData(x)
  x <- NormalizeData(x , normalization.method = "LogNormalize", scale.factor = 1000000)
  x <- x[Hom_genes_filtered$HGNC.symbol,]
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = Seu.list)

anchors <- FindIntegrationAnchors(object.list = Seu.list, anchor.features = features)

# this command creates an 'integrated' data assay
Seu_intd_obj <- IntegrateData(anchorset = anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(Seu_intd_obj) <- "integrated"

# Run the standard workflow for visualization and clustering
Seu_intd_obj <- ScaleData(Seu_intd_obj, verbose = FALSE)
Seu_intd_obj <- RunPCA(Seu_intd_obj, npcs = 30, verbose = FALSE)
Seu_intd_obj <- RunUMAP(Seu_intd_obj, reduction = "pca", dims = 1:30)
Seu_intd_obj <- FindNeighbors(Seu_intd_obj, reduction = "pca", dims = 1:30)
Seu_intd_obj <- FindClusters(Seu_intd_obj, resolution = 0.5)

# view our integration
plot2 <- DimPlot(Seu_intd_obj, reduction = "umap", group.by = "class_label")
plot1 <- DimPlot(Seu_intd_obj, reduction = "umap", group.by = "species")
plot1 + plot2

saveRDS(Seu_intd_obj, "/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/Second_Integration.rds")
Seu_intd_obj <- readRDS("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/First_Integration.rds")

#try to make sense of our clusters
unique(Seu_intd_obj$seurat_clusters)
Seu_intd_obj@meta.data$species_subclass <- paste0(Seu_intd_obj@meta.data$species, "_", Seu_intd_obj@meta.data$subclass_label)
conf_mtx <- as.data.frame.matrix(table(Seu_intd_obj$species_subclass, Seu_intd_obj$seurat_clusters))
DimPlot(Seu_intd_obj, reduction = "umap", group.by = "seurat_clusters")
FeaturePlot(Seu_intd_obj, features = c("SST", "VIP", "PVALB"))

conf_mtx <- read.csv("./conf_mtx.csv", row.names = 1)

heatmap <- heatmap(as.matrix(conf_mtx))
orig_id_order <- row.names(conf_mtx)[heatmap$rowInd]
denovo_id_order <- names(conf_mtx)[heatmap$colInd]

conf_mtx <- tibble::rownames_to_column(conf_mtx,var = "original_group")
conf_mtx_gathered <- tidyr::gather(conf_mtx, value = "cells", key = "denovo_cluster", 2:28)
conf_mtx_gathered$original_group <- factor(conf_mtx_gathered$original_group, levels = orig_id_order)
conf_mtx_gathered$denovo_cluster <- factor(conf_mtx_gathered$denovo_cluster, levels = denovo_id_order)

ggplot(conf_mtx_gathered) +
  geom_tile(aes(y = original_group, x = denovo_cluster, fill = cells)) +
  scale_fill_gradientn(limits = c(0,150),
                       colours = c("white", "red")) +
  labs(y = "Original cell type", x = "De Novo cluster", fill = "Number of cells")

#
#### More (test) analyses ####

Seu_intd_obj <- readRDS("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/First_Integration.rds")

### species specific DE

Seu_intd_obj@meta.data$species_denovoClus <- paste0(Seu_intd_obj@meta.data$species, "_", Seu_intd_obj@meta.data$seurat_clusters)
Idents(Seu_intd_obj) <- "species_denovoClus"
#DefaultAssay(Seu_intd_obj) <- "RNA"
DefaultAssay(Seu_intd_obj) <- "integrated"
DE_Genes2 <- FindMarkers(Seu_intd_obj, 
                        ident.1 = "human_8", 
                        ident.2 = "mouse_8", 
                        min.diff.pct = 0.25)

# plot
Idents(Seu_intd_obj) <- "seurat_clusters"
cluster8_pvalb <- subset(Seu_intd_obj, idents = "8")
DefaultAssay(cluster8_pvalb) <- "integrated"
Idents(cluster8_pvalb) <- "species"
avg.pvalb <- as.data.frame(log1p(AverageExpression(cluster8_pvalb)$integrated))
avg.pvalb$gene <- rownames(avg.pvalb)

genes.to.label = c("TRPM3", "SST", "TH", "COX3")
genes.to.label = rownames(DE_Genes)
p1 <- ggplot(avg.pvalb, aes(human, mouse)) + geom_point(size = 1, alpha = 0.1) + ggtitle("pvalb")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = 10)
p1

