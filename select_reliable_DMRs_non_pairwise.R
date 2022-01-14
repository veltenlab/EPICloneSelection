##################################### select_reliable_DMRs.R #####################################
#' This file computes reliable DMRs from and RnBeads output directory. We define reliable as thos
#' that are consistently methylated/unmethylated per cell type in all pairwise comparison.

library(yaml)
library(data.table)
library(RnBeads)
library(RnBeads.mm10)

report <- '/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/RnBeads/rnb_report_20211019_differential_non_pairwise/'
output <- '/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/RnBeads/DMRs/non_pariwise_reduced/'
pu1.chip <- '/users/mscherer/cluster/project/Methylome/data/external/ZhangChIP/GSM774291_FLDN1PU.1I_mm10.bed'
config_file <- '/users/mscherer/cluster/project/Methylome/src/selection_pipeline/config.yaml'
config <- yaml.load_file(config_file)

source('/users/mscherer/cluster/project/Methylome/src/selection_pipeline/checkForCutSite.R')

all.comparisons <- list.files(file.path(report,'differential_methylation_data'), full.names=TRUE, pattern = 'diffMethTable_site')
system(paste0('rm -rf ', output,'/high_*.csv'))
dmrs <- lapply(all.comparisons,function(comp){
    print(comp)
    comp <- read.csv(comp)
    fr <- checkForCutSite(comp,
                          number = 500,
                          config = config_file,
                          sort.col='mean.diff')    
    fr$X <- paste0(fr$Chromosome,'_',fr$Start)
    return(fr)
}) 
names(dmrs) <- c('HSCs', 'MPP', 'MPP1', 'MPP2')
for(n in names(dmrs)){
    write.csv(dmrs[[n]],file.path(output, paste0('high_', n, '.csv')))
}
