library(RnBeads)
diff.tab <- read.csv('/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/RnBeads/aging/rnb_report_20211109_mm10/differential_methylation_data/diffMethTable_site_cmp1.csv')
write.csv(diff.tab[abs(diff.tab$mean.diff)>0.2, ], '/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/RnBeads/aging/selected_dmrs.csv')