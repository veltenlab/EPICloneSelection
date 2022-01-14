library(WSH)
rnb.set <- load.rnb.set('/users/lvelten/project/Methylome/analysis/selection_pipeline/RnBeads/rnb_report_20210713/rnbSet_preprocessed/')
bam.file <- '/users/lvelten/project/Methylome/data/external/Cabezas/raw/GSM1274427/mapped/GSM1274427.bam/GSM1274427_1_bismark_bt2_pe.bam/GSM1274427_1_bismark_bt2_pe_sorted.bam'
parallel.setup(16)
res <- rnb.calculate.qfdrp(rnb.set,bam.file,cores=16)
write.csv(res,'/users/lvelten/project/Methylome/analysis/selection_pipeline/WSH/GSM1274427_qFDRP.csv')
