# "/users/lvelten/project/SCG4SYN/LibraryDesign/Medoid_Motifs/Selected_TFBS/200423_Mouse_selected_TFBS_medoid_motifs.rds"
# sind das alle 38 Factors vom screen?
# "/users/lvelten/project/SCG4SYN/LibraryDesign/Medoid_Motifs/200713_update_medoids_and_resources/200713_mouse_medoid_TFBS_for_each_892_genes.rds"
# das sind alle medoids von den 892 tfs in mouse
# "/users/lvelten/project/SCG4SYN/LibraryDesign/Medoid_Motifs/TFs_in_DB.zip" enthält die raw daten von den 3 datenbanken
.libPaths(c(.libPaths(), '/users/mscherer/conda/envs/rnbeads/'))
library(TFBSTools)
dat <- readRDS('/users/mscherer/cluster/project/SCG4SYN/LibraryDesign/Medoid_Motifs/Selected_TFBS/200423_Mouse_selected_TFBS_medoid_motifs.rds')
names(dat) <- lapply(strsplit(names(dat), '-'), function(x)x[[2]])
max_motifs <- lapply(dat,function(pwm){
  mat <- pwm@profileMatrix
  paste0(apply(mat, 2, function(x){
    names(which.max(x))
  }), collapse='')
})
max_motifs <- data.frame(TF=names(max_motifs), Motif=unlist(max_motifs))
write.csv(max_motifs, '/users/mscherer/cluster/project/Methylome/analysis/selection_pipeline/Robert_max_TFBSs.csv')
