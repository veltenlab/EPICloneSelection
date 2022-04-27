library(Biostrings)
library(ggplot2)
tag_seq <- reverse(DNAString('GAAATTCCGGCCAGGATCGTTCTGATGACGCTCA')) # the new antibody-B tag
mb_id <- '4318'
ampli_file <- read.csv(paste0('/users/lvelten/project/Methylome/analysis/selection_pipeline/MB_output/Designer/', mb_id, '-design-summary_annotated.csv'))
res <- c()
for(i in 1:nrow(ampli_file)){
  id <- ampli_file[i, 'AmpID']
  rev_seq <- as.character(DNAString(DNAString(ampli_file[i, 'rev_seq'])))
  fwd_seq <- as.character(tag_seq)
  stri <- paste('PRIMER_TASK=check_primers',
                'PRIMER_MAX_SIZE=36',
                paste0('SEQUENCE_ID=', id),
                paste0('SEQUENCE_PRIMER=', fwd_seq),
                paste0('SEQUENCE_PRIMER_REVCOMP=', rev_seq),
                'PRIMER_EXPLAIN_FLAG=1',
                'PRIMER_MAX_TM=100', 
                'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1',  
                'PRIMER_SECONDARY_STRUCTURE_ALIGNMENT=1', 
                '=',
                sep='\n')
  tmp <- tempdir() 
  fi <- paste0(tmp, '/primer3')
  writeLines(stri, fi)
  out <- paste0(tmp, '/out')
  primer_3_out <- system(paste('primer3_core',
                fi), 
          intern = TRUE)
  any_th <- primer_3_out[grepl('_COMPL_ANY_TH', primer_3_out)]
  any_th <- unlist(strsplit(any_th, '_COMPL_ANY_TH='))[2]
  end_th <- primer_3_out[grepl('_COMPL_END_TH', primer_3_out)]
  end_th <- unlist(strsplit(end_th, '_COMPL_END_TH='))[2]
  res <- c(res, max(as.numeric(c(end_th, any_th))))
}
to_plot <- data.frame(AmplID=ampli_file$AmpID, TM=res)
plot <- ggplot(to_plot, aes(x=TM, y=..count..))+geom_histogram()+theme_minimal()
ggsave('TM_histogram_reverse.pdf', plot)
write.csv(to_plot, 'TM_reverse.csv')
ampli_file$TM_rev <- res

res <- c()
for(i in 1:nrow(ampli_file)){
  id <- ampli_file[i, 'AmpID']
  rev_seq <- as.character(DNAString(DNAString(ampli_file[i, 'fwd_seq'])))
  fwd_seq <- as.character(tag_seq)
  stri <- paste('PRIMER_TASK=check_primers',
                'PRIMER_MAX_SIZE=36',
                paste0('SEQUENCE_ID=', id),
                paste0('SEQUENCE_PRIMER=', fwd_seq),
                paste0('SEQUENCE_PRIMER_REVCOMP=', rev_seq),
                'PRIMER_EXPLAIN_FLAG=1',
                'PRIMER_MAX_TM=100', 
                'PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1',  
                'PRIMER_SECONDARY_STRUCTURE_ALIGNMENT=1', 
                '=',
                sep='\n')
  tmp <- tempdir() 
  fi <- paste0(tmp, '/primer3')
  writeLines(stri, fi)
  out <- paste0(tmp, '/out')
  primer_3_out <- system(paste('primer3_core',
                               fi), 
                         intern = TRUE)
  any_th <- primer_3_out[grepl('_COMPL_ANY_TH', primer_3_out)]
  any_th <- unlist(strsplit(any_th, '_COMPL_ANY_TH='))[2]
  end_th <- primer_3_out[grepl('_COMPL_END_TH', primer_3_out)]
  end_th <- unlist(strsplit(end_th, '_COMPL_END_TH='))[2]
  res <- c(res, max(as.numeric(c(end_th, any_th))))
}
to_plot <- data.frame(AmplID=ampli_file$AmpID, TM=res)
plot <- ggplot(to_plot, aes(x=TM, y=..count..))+geom_histogram()+theme_minimal()
ggsave('TM_histogram_forward.pdf', plot)
write.csv(to_plot, 'TM_forward.csv')
ampli_file$TM_fwd <- res
write.csv(ampli_file, paste0('/users/lvelten/project/Methylome/analysis/selection_pipeline/MB_output/Designer/', mb_id, '-design-summary_annotated.csv'))
