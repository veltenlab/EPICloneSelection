# Investigate dropout in relationship to GC content
library(BSgenome.Hsapiens.UCSC.hg19)
library(RnBeads)
genome <- BSgenome.Hsapiens.UCSC.hg19
ampli.info <- read.table("/users/lvelten/project/Methylome/infos/BCells/Blood.Bone.Marrow.Amplicons.design.dropout.added.tsv")
ampli.gr <- makeGRangesFromDataFrame(ampli.info[,c('chr','amplicon_start','amplicon_end')],start.field='amplicon_start',end.field='amplicon_end')
gc.content <- apply(ampli.info,1,function(x){
  seq <- genome[[x[1]]]
  sel.seq <- seq[x[2]:x[5]]
  freq <- alphabetFrequency(sel.seq)
  (freq['C']+freq['G'])/length(sel.seq)
})
ampli.info$GC_content <- gc.content
cpgs <- makeGRangesFromDataFrame(rnb.annotation2data.frame(rnb.get.annotation("CpG", assembly = 'hg19')))
op <- findOverlaps(ampli.gr,cpgs)
cpg.counts <- plyr::count(queryHits(op))$freq/2
ampli.info$CpG_count <- cpg.counts
write.csv(ampli.info,"/users/lvelten/project/Methylome/analysis/selection_pipeline/MB_output/amplicons_annotated.csv")
