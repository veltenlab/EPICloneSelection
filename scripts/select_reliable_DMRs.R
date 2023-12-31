##################################### select_reliable_DMRs.R #####################################
#' This file computes reliable DMRs from and RnBeads output directory. We define reliable as thos
#' that are consistently methylated/unmethylated per cell type in all pairwise comparison.

#.libPaths(c('/users/mscherer/R/R-4.1.2',.libPaths()))

library(yaml)
library(data.table)
library(RnBeads)
library(RnBeads.mm10)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-r", "--report", type="character")
parser$add_argument("-o", "--output", type="character")
parser$add_argument("-d", "--digest", type="character")
parser$add_argument("-c", "--config", type="character", default="config.yaml")
args <- parser$parse_args()

report <- args$report
output <- args$output
config_file <- args$config
config <- yaml.load_file(config_file)
cut_file <- args$digest

source(cut_file)

all.comparisons <- list.files(file.path(report,'differential_methylation_data'), full.names=TRUE, pattern = 'diffMethTable_site')
system(paste0('rm -rf ', output,'/high_*.csv'))
system(paste0('rm -rf ', output,'/low_*.csv'))
dmrs <- lapply(all.comparisons,function(comp){
    print(comp)
    diff.table <- as.data.frame(fread(comp))
    group.names <- colnames(diff.table)[grepl('sd.',colnames(diff.table))]
    group.names <- gsub('sd.','', group.names)
    print(paste('Processing diffMethTable for ', group.names[1], 'and', group.names[2]))
    diff.table <- diff.table[abs(diff.table$mean.diff)>config[['dmcs']][['min_diff']],]
    row.names(diff.table) <- paste0(diff.table$Chromosome, '_', diff.table$Start)
    diff.table[['comparison']] <- paste(group.names[1], 'vs.', group.names[2])
    diff.positive <- diff.table[diff.table$mean.diff>0,]
    diff.positive$type <- paste0('high_', group.names[1])
    row.names(diff.positive) <- paste0(diff.positive$Chromosome, '_', diff.positive$Start)
    file.positive <- file.path(output,paste0('high_', group.names[1], '.csv'))
    if(file.exists(file.positive)){
        all.positive <- read.csv(file.positive,row.names = 1)
        if(nrow(all.positive)<4){
            stop('Please delete the high_* files in the output directory')
        }
        both <- intersect(row.names(diff.positive),row.names(all.positive))
        all.positive <- data.frame(all.positive[both,],diff.positive[both,])
        write.csv(all.positive,file.positive)
    }else{
        write.csv(diff.positive,file.positive)
    }
    return(diff.positive)
})

all.files <- list.files(output, full.names = TRUE, pattern = 'high')
dmrs.all <- lapply(all.files, function(dmr){
    print(dmr)
    dmr <- read.csv(dmr)
    dmr <- checkForCutSite(dmr,
                    number=min(nrow(dmr), 500),
                    config=config_file,
                    sort.col=c('mean.diff'),
                    use.extended=TRUE)
    print(nrow(dmr))
    tfbs_sites <- colnames(dmr)[(which(colnames(dmr)=='GCContent')+1):ncol(dmr)]
    tfbs_frame <- dmr[, tfbs_sites]
    all_nas <- apply(tfbs_frame, 1, function(x)all(is.na(x)))
    print(sum(!all_nas))
    dmr <- dmr[!all_nas, ]
})
all.names <- list.files(output, pattern = 'high')
all.names <- gsub('high_','high_filtered_',all.names)
for(i in 1:length(dmrs.all)){
    write.csv(dmrs.all[[i]],file.path(output,all.names[i]))
}

all.dmrs <- data.frame(MeanDiff=abs(c(dmrs.all[[1]]$mean.diff,dmrs.all[[2]]$mean.diff,dmrs.all[[3]]$mean.diff,dmrs.all[[4]]$mean.diff)),
                       Type=c(dmrs.all[[1]]$type,dmrs.all[[2]]$type,dmrs.all[[3]]$type,dmrs.all[[4]]$type))
theme <- theme(panel.background = element_rect(color='black',fill='white'),
               panel.grid=element_blank(),
               text=element_text(color='black',size=15))
plot <- ggplot(all.dmrs,aes(x=Type,y=MeanDiff))+geom_boxplot()+theme
ggsave(file.path(output, 'DMR_mean_diff.pdf'), plot)
