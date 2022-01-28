# Investigate dropout in relationship to GC content
library(BSgenome.Mmusculus.UCSC.mm10)
library(RnBeads)
genome <- BSgenome.Mmusculus.UCSC.mm10
ampli.info <- read.csv("/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/MB_output/4030-design-summary.csv")
ampli.gr <- makeGRangesFromDataFrame(ampli.info[,c('chr','amplicon_start','amplicon_end')],start.field='amplicon_start',end.field='amplicon_end')
gc.content <- apply(ampli.info,1,function(x){
  seq <- genome[[x[2]]]
  sel.seq <- seq[x[3]:x[6]]
  freq <- alphabetFrequency(sel.seq)
  (freq['C']+freq['G'])/length(sel.seq)
})
ampli.info$GC_content <- gc.content
cpgs <- makeGRangesFromDataFrame(rnb.annotation2data.frame(rnb.get.annotation("CpG", assembly = 'mm10')))
op <- findOverlaps(ampli.gr,cpgs)
cpg.counts <- plyr::count(queryHits(op))$freq/2
ampli.info$CpG_count <- cpg.counts
write.csv(ampli.info,"/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/MB_output/amplicons_annotated.csv")
