####################### generate_final_panel.R ####################### 
#' With this file, we merge together all the selected CpGs and construct
#' the final table for upload into the MissionBio Designer

library(GenomicRanges)
n_hscs <- 105
# Overthink this!
n_mpp1 <- 70
n_mpp2 <- 70
n_mpp <- 75
n_imrs <- 210
n_wsh <- 80
n_non_cut <- 50
n_meth <- 20
n_unmeth <- 20
n_primers <- 10

hsc_dmrs <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/RnBeads/DMRs/high_filtered_extended_HSC.csv')
mpp1_dmrs <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/RnBeads/DMRs/high_filtered_extended_MPP1.csv')
mpp2_dmrs <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/RnBeads/DMRs/high_filtered_extended_MPP2.csv')
mpp_dmrs <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/RnBeads/DMRs/high_filtered_extended_MPP.csv')
imrs <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/IMS/IMS_annotated.csv')
wsh <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/WSH/filtered_qFDRP/filtered_qFDRP.csv')
non_cut <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/non_cut_amplicons.csv')
always_meth <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/RnBeads/meth_unmeth/always_meth_filtered.csv')
always_unmeth <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/RnBeads/meth_unmeth/always_unmeth_filtered.csv')

all_frames <- list(HSC=hsc_dmrs,
                   MPP1=mpp1_dmrs,
                   MPP2=mpp2_dmrs,
                   MPP=mpp_dmrs,
                   IMR=imrs,
                   WSH=wsh,
                   NonCut=non_cut,
                   AlwaysMeth=always_meth,
                   AlwaysUnmeth=always_unmeth)
all_frames <- lapply(all_frames, function(x){
  x <- x[x$GCContent<0.6, ]
  if(!('End'%in%colnames(x))){
    x[, 'End'] <- x[, 'Start']+1
  }
  x
})
for(type in c('HSC', 'MPP', 'IMR', 'WSH')){
  fr <- all_frames[[type]]
  is_non_na <- !is.na(fr[, 'enhancer_annotation']) | !is.na(fr[, 'promoter']) | !is.na(fr[, 'AgingDMC'])
  is_non_na[!is_non_na][sample(1:sum(!is_non_na), round(length(is_non_na)*(sum(is_non_na)*0.25/length(is_non_na))))] <- TRUE
  fr <- fr[is_non_na, ]
  all_frames[[type]] <- fr
}
final_frame <- data.frame(rbind(all_frames$HSC[1:n_hscs, c('Chromosome', 'Start', 'End')],
                                all_frames$MPP1[1:n_mpp1, c('Chromosome', 'Start', 'End')],
                                all_frames$MPP2[1:n_mpp2, c('Chromosome', 'Start', 'End')],
                                all_frames$MPP[1:n_mpp, c('Chromosome', 'Start', 'End')],
                                all_frames$IMR[1:n_imrs, c('Chromosome', 'Start', 'End')],
                                all_frames$WSH[1:n_wsh, c('Chromosome', 'Start', 'End')],
                                all_frames$NonCut[1:n_non_cut, c('Chromosome', 'Start', 'End')],
                                all_frames$AlwaysMeth[1:n_meth, c('Chromosome', 'Start', 'End')],
                                all_frames$AlwaysUnmeth[1:n_unmeth, c('Chromosome', 'Start', 'End')]
                                ),
                          Name=c(paste0(rep('HSC_high', n_hscs), 1:n_hscs),
                                 paste0(rep('MPP1_high', n_mpp1), 1:n_mpp1),
                                 paste0(rep('MPP2_high', n_mpp2), 1:n_mpp2),
                                 paste0(rep('MPP_high', n_mpp), 1:n_mpp),
                                 paste0(rep('IMR', n_imrs), 1:n_imrs),
                                 paste0(rep('WSH', n_wsh), 1:n_wsh),
                                 paste0(rep('Non_cut', n_non_cut), 1:n_non_cut),
                                 paste0(rep('Always_methylated', n_meth), 1:n_meth),
                                 paste0(rep('Always_unmethylated', n_unmeth), 1:n_unmeth)),
                          Type=c(rep('HSC_high', n_hscs),
                                 rep('MPP1_high', n_mpp1),
                                 rep('MPP2_high', n_mpp2),
                                 rep('MPP_high', n_mpp),
                                 rep('IMR', n_imrs),
                                 rep('WSH', n_wsh),
                                 rep('Non_cut', n_non_cut),
                                 rep('Always_methylated', n_meth),
                                 rep('Always_unmethylated', n_unmeth)))
write.csv(final_frame, 
          '/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/MB_input/CpGs_only.csv',
          row.names=FALSE,
          quote=FALSE)
saveRDS(all_frames, '/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/MB_input/combined_upload.RDS')
final_frame[!grepl('Non_cut', final_frame$Type), 'Start'] <- final_frame[!grepl('Non_cut', final_frame$Type), 'Start']-10
final_frame[!grepl('Non_cut', final_frame$Type), 'End'] <- final_frame[!grepl('Non_cut', final_frame$Type), 'End']+10
final_frame[grepl('Non_cut', final_frame$Type), 'Start'] <- final_frame[grepl('Non_cut', final_frame$Type), 'Start']+175
final_frame[grepl('Non_cut', final_frame$Type), 'End'] <- final_frame[grepl('Non_cut', final_frame$Type), 'End']-175
final_frame_gr <- makeGRangesFromDataFrame(final_frame)
op <- findOverlaps(final_frame_gr, final_frame_gr)
queryHits(op)[queryHits(op)!=subjectHits(op)]
subjectHits(op)[queryHits(op)!=subjectHits(op)]
to_write <- data.frame(rep('Region', nrow(final_frame)),
                       paste0(final_frame$Chromosome, ':', final_frame$Start, '-', final_frame$End),
                       final_frame$Name)
write.csv(to_write, 
            '/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/MB_input/upload_all_samples.csv',
          row.names=FALSE,
          quote=FALSE)