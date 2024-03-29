library(Seurat)
library(biomaRt)
library(dplyr)

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

Hom_genes_filtered = Homologous_genes_mouhum[Homologous_genes_mouhum$One_to_one == T,] 

#
#### raw data ####

# metadata
hum_meta <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/human/metadata.csv")
mou_meta <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/mouse/2020_10_24/metadata.csv")

hum_meta <- hum_meta[, intersect(colnames(hum_meta), colnames(mou_meta))]
mou_meta <- mou_meta[, intersect(colnames(hum_meta), colnames(mou_meta))]

hum_meta$species <- "human"
mou_meta$species <- "mouse"

hum_meta <- hum_meta[hum_meta$outlier_call == "False",]
mou_meta <- mou_meta[mou_meta$outlier_call == "False",]

write.csv(hum_meta, "/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_human_meta.csv")
write.csv(mou_meta, "/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_mouse_meta.csv")

# count matrices
hum_counts <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/human/matrix.csv", row.names = 1)
mou_counts <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/mouse/2020_10_24/matrix.csv", row.names = 1)

hum_counts <- hum_counts[hum_counts$sample_name %in% hum_meta$sample_name,]
mou_counts <- mou_counts[mou_counts$sample_name %in% mou_meta$sample_name,]

Hom_genes_filtered <- Hom_genes_filtered[Hom_genes_filtered$MGI.symbol %in% colnames(mou_counts),]
Hom_genes_filtered <- Hom_genes_filtered[Hom_genes_filtered$HGNC.symbol %in% colnames(hum_counts),]
write.csv(Hom_genes_filtered, "./Homologous_genes_mouhum_filtered.csv")

hum_counts <- hum_counts[,intersect(colnames(hum_counts), Hom_genes_filtered$HGNC.symbol)]
write.csv(hum_counts, "/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_human_counts.csv")

mou_counts <- mou_counts[,intersect(colnames(mou_counts), Hom_genes_filtered$MGI.symbol)]
write.csv(mou_counts, "/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_mouse_counts.csv")

#### Seurat object ####

row.names(hum_meta) <- hum_meta$sample_name
Seu_hum <- CreateSeuratObject(counts = t(hum_counts), meta.data = hum_meta)
saveRDS(Seu_hum, "/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/Seu_hum.rds")

row.names(mou_meta) <- mou_meta$sample_name
Seu_mou <- CreateSeuratObject(counts = t(mou_counts), meta.data = mou_meta)
saveRDS(Seu_mou, "/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/Seu_mou.rds")

#### Subsample dataset (raw data V2) ####

#human
Seu_hum <- readRDS("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/Seu_hum.rds")
Idents(Seu_hum) <- "subclass_label"
Seu_hum_mini <- subset(Seu_hum, downsample = 150)
saveRDS(Seu_hum_mini, "/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/Seu_hum_mini.rds")
sample_list <- unname(Seu_hum_mini$sample_name)

hum_counts <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/human/matrix.csv", row.names = 1)
hum_counts <- hum_counts[sample_list,]
write.csv(hum_counts, "/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_human_counts_mini.csv")

hum_meta <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_human_meta.csv", row.names = 1)
hum_meta <- hum_meta[hum_meta$sample_name %in% sample_list,]
write.csv(hum_meta, "/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_human_meta_mini.csv")
row.names(hum_meta) <- hum_meta$sample_name
Seu_hum_mini <- CreateSeuratObject(counts = t(hum_counts), meta.data = hum_meta)
saveRDS(Seu_hum_mini, "/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/Seu_hum_mini.rds")

#mouse
Seu_mou <- readRDS("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/Seu_mou.rds")
Idents(Seu_mou) <- "subclass_label"
Seu_mou_mini <- subset(Seu_mou, downsample = 150)
Seu_mou_mini <- subset(Seu_mou_mini, downsample = 80)
saveRDS(Seu_mou_mini, "/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/Seu_mou_mini.rds")
sample_list <- unname(Seu_mou_mini$sample_name)

mou_counts <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/mouse/2020_10_24/matrix.csv", row.names = 1)
mou_counts <- mou_counts[sample_list,]
write.csv(mou_counts, "/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_mouse_counts_mini.csv")

mou_meta <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_mouse_meta.csv", row.names = 1)
mou_meta <- mou_meta[mou_meta$sample_name %in% sample_list,]
write.csv(mou_meta, "/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_mouse_meta_mini.csv")
row.names(mou_meta) <- mou_meta$sample_name
Seu_mou_mini <- CreateSeuratObject(counts = t(mou_counts), meta.data = mou_meta)
saveRDS(Seu_mou_mini, "/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/Seu_mou_mini.rds")

#genesum filter

hum_counts <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_human_counts_mini.csv", row.names = 1)
hum_counts <- hum_counts[,colSums(hum_counts) != 0]
write.csv(hum_counts, "/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_human_counts_mini.csv")
hum_meta <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_human_meta_mini.csv", row.names = 1)
row.names(hum_meta) <- hum_meta$sample_name
Seu_hum_mini <- CreateSeuratObject(counts = t(hum_counts), meta.data = hum_meta)
saveRDS(Seu_hum_mini, "/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/Seu_hum_mini.rds")

mou_counts <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_mouse_counts_mini.csv")
mou_counts <- mou_counts[!is.na(mou_counts$X),] 
row.names(mou_counts) <- mou_counts$X
mou_counts <- mou_counts[,-1]
mou_counts <- mou_counts[,colSums(mou_counts) != 0]
table(colSums(mou_counts) != 0)
write.csv(mou_counts, "/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_mouse_counts_mini.csv")
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
  
Seu_mou_mini <- CreateSeuratObject(counts = mou_counts, meta.data = mou_meta)
saveRDS(Seu_mou_mini, "/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/Seu_mou_mini.rds")
