

hum_counts <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_human_counts_mini.csv", row.names = 1)
hum_meta <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Intralab_collab_scc_projects/KCNISS2002_week2proj/AIBS_human_meta_mini.csv", row.names = 1)

hum_meta <- hum_meta[sample(100),]
hum_counts <- hum_counts[hum_meta$sample_name,]

library(limma)
library(edgeR)
test <- DGEList(t(hum_counts))
