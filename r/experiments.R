# Run all
###############################################################################
source("./source.R")
cat("\n\n\n\n")
#experimentRunAll()

# motif M4 is too dense and it disconnects the graph
motif_names = c('Ms','Md','M1','M2','M3','M5','M6','M7','M8','M9','M10','M11','M12','M13')

eigenvectors_for_map = c(2:7)
n_vects_for_cluster = 7
n_clusters = 7
simplify_sweep = 1
for(weighted in c(TRUE, FALSE)){
  experimentMigration(motif_names, eigenvectors_for_map, n_vects_for_cluster, n_clusters, simplify_sweep, weighted)
}
