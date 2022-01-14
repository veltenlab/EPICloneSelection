##################################### find_IMS.R #####################################
#' This script is used to deteremine intermediately methylated sites from the RnBeads
#' set, focusing only on the HSCs

.libPaths(c('/users/mscherer/conda/envs/rnbeads/lib/R/library/', .libPaths()))
library(yaml)
library(data.table)
library(RnBeads)
library(RnBeads.mm10)

rnb_set_path <- '/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/RnBeads/rnb_report_20211004_reduced/cluster_run/preprocessing_RnBSet/'
out_folder <- '/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/IMS/'
all_dmrs <- c('/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/RnBeads/DMRs/non_pairwise/high_HSCs.csv',
              '/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/RnBeads/DMRs/non_pairwise/high_MPP.csv',
              '/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/RnBeads/DMRs/non_pairwise/high_MPP1.csv',
              '/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/RnBeads/DMRs/non_pairwise/high_MPP2.csv')
config_file <- '/users/mscherer/cluster/project/Methylome/src/selection_pipeline/config.yaml'
config <- yaml.load_file(config_file)

source('/users/mscherer/cluster/project/Methylome/src/selection_pipeline/checkForCutSite.R')

rnb_set <- load.rnb.set(rnb_set_path)
anno_fr <- annotation(rnb_set)
anno_gr <- makeGRangesFromDataFrame(anno_fr)
dmrs <- do.call(rbind.fill, lapply(all_dmrs, read.csv))
dmrs_gr <- makeGRangesFromDataFrame(dmrs, end='Start')
op <- findOverlaps(anno_gr, dmrs_gr)
anno_gr <- anno_gr[-queryHits(op)]
meth_data <- meth(rnb_set)[, c("GSM1274424",
                               "GSM1274425",
                               "GSM1274426")]
meth_data <- meth_data[-queryHits(op), ]
is_intermediate <- apply(meth_data, 1, function(x){
  all(x>0.25&x<0.75)
})
meth_data <- meth_data[is_intermediate, ]
anno_gr <- anno_gr[is_intermediate]
covg_data <- covg(rnb_set)[, c("GSM1274424",
                               "GSM1274425",
                               "GSM1274426")]
covg_data <- covg_data[-queryHits(op), ][is_intermediate, ]
mean_covg <- rowMeans(covg_data)
too_high <- mean_covg>quantile(mean_covg, .95)
mean_covg <- mean_covg[!too_high]
meth_data <- meth_data[!too_high, ]
anno_gr <- anno_gr[!too_high]
meth_data_fr <- data.frame(Chromosome=seqnames(anno_gr),
                    Start=start(anno_gr),
                    End=end(anno_gr),
                    MeanMeth=rowMeans(meth_data),
                    MeanCovg=mean_covg)
res <- checkForCutSite(meth_data_fr,
                       number=1000,
                       config=config_file, 
                       sort.col='MeanCovg')
if(!dir.exists(out_folder)){
  system(paste('mkdir', out_folder))
}
tfbs_sites <- colnames(res)[(which(colnames(res)=='GCContent')+1):ncol(res)]
tfbs_frame <- res[, tfbs_sites]
all_nas <- apply(tfbs_frame, 1, function(x)all(is.na(x)))
res <- res[!all_nas, ]
write.csv(res, paste0(out_folder, '/IMS_annotated.csv'))
