##################################### find_IMS.R #####################################
#' This script is used to deteremine intermediately methylated sites from the RnBeads
#' set, focusing only on the HSCs

library(yaml)
library(data.table)
library(RnBeads)
library(RnBeads.mm10)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-r", "--rnbset", type="character")
parser$add_argument("-o", "--output", type="character")
parser$add_argument("-d", "--dmrs", type="character")
parser$add_argument("-p", "--pdrs", type="character")
parser$add_argument("-c", "--config", type="character", default="config.yaml")
parser$add_argument("-a", "--digest", type="character")
args <- parser$parse_args()

rnb_set_path <- args$rnbset
out_folder <- args$output
all_dmrs <- list.files(args$dmrs, pattern='high_filtered_', full.names=TRUE)
pdrs <- sapply(list.dirs(args$pdrs, recursive=FALSE), function(x){
  list.files(file.path(x, 'PDR'), pattern='.csv', full.names=TRUE)
})
pdr_annotations <- sapply(list.dirs(args$pdrs, recursive=FALSE), function(x){
  list.files(file.path(x, 'PDR'), pattern='.RData', full.names=TRUE)
})
config_file <- args$config
config <- yaml.load_file(config_file)
load(pdr_annotations[1])
for(i in 2:length(pdr_annotations)){
  last_anno <- annotation
  load(pdr_annotations[i])
  op <- findOverlaps(last_anno, annotation)
  last_anno <- last_anno[queryHits(op)]
}
all_pdr <- matrix(nrow=length(last_anno), ncol=length(pdrs))
for(i in 1:length(pdr_annotations)){
  load(pdr_annotations[i])
  op <- findOverlaps(last_anno, annotation)
  pdr <- read.csv(pdrs[i])
  all_pdr[queryHits(op), i] <- pdr[subjectHits(op), ]
}
annotation <- last_anno
cut_file <- args$digest

source(cut_file)

rnb_set <- load.rnb.set(rnb_set_path)
anno_fr <- annotation(rnb_set)
anno_gr <- makeGRangesFromDataFrame(anno_fr)
dmrs <- do.call(rbind.fill, lapply(all_dmrs, read.csv))
dmrs_gr <- makeGRangesFromDataFrame(dmrs, end='Start')
dmrs_gr <- resize(dmrs_gr, width = 500, fix = 'center')
op <- findOverlaps(anno_gr, dmrs_gr)
anno_gr <- anno_gr[-queryHits(op)]
meth_data <- meth(rnb_set)[, c("GSM1274424"), drop=FALSE]#,
#                               "GSM1274425",
#                               "GSM1274426"), drop=FALSE]
meth_data <- meth_data[-queryHits(op), , drop=FALSE]
is_intermediate <- apply(meth_data, 1, function(x){
  all(x>0.25&x<0.75)
})
meth_data <- meth_data[is_intermediate, , drop=FALSE]
anno_gr <- anno_gr[is_intermediate]
covg_data <- covg(rnb_set)[, c("GSM1274424"), drop=FALSE]#,
#                               "GSM1274425",
#                              "GSM1274426"), drop=FALSE]
covg_data <- covg_data[-queryHits(op), , drop=FALSE][is_intermediate, , drop=FALSE]
mean_covg <- rowMeans(covg_data)
too_high <- mean_covg>quantile(mean_covg, .95)
mean_covg <- mean_covg[!too_high]
meth_data <- meth_data[!too_high, , drop=FALSE]
anno_gr <- anno_gr[!too_high]
op <- findOverlaps(anno_gr, annotation)
pdrs <- rowMeans(all_pdr)
pdrs <- pdrs[subjectHits(op)]
meth_data_fr <- data.frame(Chromosome=seqnames(anno_gr),
                    Start=start(anno_gr),
                    End=end(anno_gr),
                    MeanMeth=rowMeans(meth_data),
                    MeanCovg=mean_covg)
meth_data_fr[queryHits(op), 'PDR'] <- pdrs
res <- checkForCutSite(na.omit(meth_data_fr),
                       number=config[['general']][['n_imcs']],
                       config=config_file, 
                       sort.col='PDR',
                       decreasing=FALSE,
                       use.extended=TRUE)
if(!dir.exists(out_folder)){
  system(paste('mkdir', out_folder))
}
write.csv(res, paste0(out_folder, '/IMC_annotated_all.csv'))
tfbs_sites <- colnames(res)[(which(colnames(res)=='GCContent')+1):ncol(res)]
tfbs_frame <- res[, tfbs_sites]
all_nas <- apply(tfbs_frame, 1, function(x)all(is.na(x)))
res <- res[!all_nas, ]
write.csv(res, paste0(out_folder, '/IMC_annotated.csv'))
