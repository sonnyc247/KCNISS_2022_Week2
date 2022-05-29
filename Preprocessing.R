#### Hom genes ####

counts_matrix <- counts_matrix[-grep('[[:digit:]]+[[:digit:]]+[[:digit:]]+[[:digit:]]+', rownames(counts_matrix)),]

#### gene homology (https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/) ####

# Basic function to convert mouse to human gene names

Seu_Mouse_data <- readRDS("/external/rprshnas01/kcni/ychen/git/Ex_Env_Storage/Nr1_celltypes/DGE - Seurat AIBS Object.rds")
testgenes <- row.names(Seu_Mouse_data)

require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
MGI_HGNC = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = testgenes, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
dup_hum = MGI_HGNC[duplicated(MGI_HGNC$HGNC.symbol), "HGNC.symbol"]
MGI_HGNC$Duplicate_MGI_to_HGNC <- MGI_HGNC$HGNC.symbol %in% dup_hum
MGI_HGNC_filtered = MGI_HGNC[!(MGI_HGNC$HGNC.symbol %in% dup_hum),1:2] 

write.csv(MGI_HGNC, "./MGI_HGNC.csv")
