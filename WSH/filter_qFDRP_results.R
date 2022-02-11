##################################### annotate_WSH.R #####################################
#' With this script, we filter for potentially interesting variable regions and annotate them
#' according to further epigenomic annotations such as TFBSs or ChIP-seq data

.libPaths(c('/users/mscherer/conda/envs/rnbeads/lib/R/library/', .libPaths()))
library(yaml)
library(data.table)
library(RnBeads)
library(RnBeads.mm10)

res_folder <- '/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/WSH_subsampled/GSM1274424/'
sample_name <- 'qFDRP_GSM1274424'
out_folder <- '/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/WSH/'
all_dmrs <- c('/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/RnBeads/DMRs/non_pairwise/high_HSCs.csv',
              '/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/RnBeads/DMRs/non_pairwise/high_MPP.csv',
              '/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/RnBeads/DMRs/non_pairwise/high_MPP1.csv',
              '/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/RnBeads/DMRs/non_pairwise/high_MPP2.csv')
config_file <- '/users/mscherer/cluster/project/Methylome/src/selection_pipeline/config.yaml'
config <- yaml.load_file(config_file)

source('/users/mscherer/cluster/project/Methylome/src/selection_pipeline/checkForCutSite.R')

qfdrp <- read.csv(paste0(res_folder, sample_name, '.csv'))
load(paste0(res_folder, 'annotation.RData'))
nas <- is.na(qfdrp$x)
qfdrp <- qfdrp[!nas, ]
annotation <- annotation[!nas]
is_1 <- qfdrp==1
annotation <- annotation[!is_1]
qfdrp <- qfdrp[!is_1]
dmrs <- do.call(rbind.fill, lapply(all_dmrs, read.csv))
dmrs_gr <- makeGRangesFromDataFrame(dmrs, end='Start')
dmrs_gr <- resize(dmrs_gr, width = 500, fix = 'center')
op <- findOverlaps(annotation,dmrs_gr)
qfdrp <- qfdrp[-queryHits(op)]
annotation <- annotation[-queryHits(op)]
qfdrp <- data.frame(Chromosome=seqnames(annotation),
                    Start=start(annotation),
                    End=end(annotation),
                    qFDRP=qfdrp)
res <- checkForCutSite(qfdrp,
                       number=230,
                       config=config_file, 
                       sort.col='qFDRP')
out_folder <- file.path(out_folder, paste0('filtered_', sample_name))
if(!dir.exists(out_folder)){
  system(paste('mkdir', out_folder))
}
tfbs_sites <- colnames(res)[(which(colnames(res)=='GCContent')+1):ncol(res)]
tfbs_frame <- res[, tfbs_sites]
all_nas <- apply(tfbs_frame, 1, function(x)all(is.na(x)))
res <- res[!all_nas, ]
write.csv(res, paste0(out_folder, '/filtered_', sample_name, '.csv'))
