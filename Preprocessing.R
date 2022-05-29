library(Seurat)
library(biomaRt)

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





