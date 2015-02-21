library(ape)
library(ClockstaRG)
library(rgl)

get_map_phylo <- function(log_file, trees_file){
    map_location <- which.max(log_file$clockRate)
    map_chrono <- trees_file[[map_location]]
    map_phylo <- map_chrono
    map_phylo$edge.length <- map_phylo$edge.length * log_file$clockRate[map_location]
    return(map_phylo)
}

get_brs_mat <- function(tree_list){
  br_lens_matrix <- matrix(NA, nrow = length(tree_list), ncol = length(tree_list[[1]]$edge.length))
  rownames(br_lens_matrix) <- names(tree_list)
  colnames(br_lens_matrix) <- paste0('br', 1:length(tree_list[[1]]$edge.length))
  for(i in 1:length(tree_list)){
    br_lens_matrix[i, ] <- ladderize(tree_list[[i]])$edge.length
  }
  return(br_lens_matrix)
}

# Reconstruct phylogram for true data
true_log <- read.table('negcontr_gmc/gmc_negcont_gmc.log', sep = ',', head = T)
true_trees <- read.nexus('negcontr_gmc/simulationduchene.trees')
true_map <- get_map_phylo(true_log, true_trees)

# Reconstruct phylograms for PPS
pps_files <- system('ls -d negcontr_gmc/pps/*', intern = T)
pps_files <- pps_files[-(grep('[.]', pps_files))]

pps_phylogs <- list()
for(i in 1:length(pps_files)){
    cat('reading data', pps_files[i], '\n')
    pps_log <- read.table(paste0(pps_files[i], '/pps.log'), sep = ',', head = T)
    pps_tree <- read.nexus(paste0(pps_files[i], '/simulationduchene.trees'))
    pps_phylogs[[i]] <- get_map_phylo(pps_log, pps_tree)
}

all_trees <- list()
all_trees[[1]] <- true_map
for(i in 1:length(pps_phylogs)){
    all_trees[[i + 1]] <- pps_phylogs[[i]]
}

brs_mat <- get_brs_mat(all_trees)

brs_quants <- sapply(1:nrow(brs_mat), function(x) abs(0.5 - sum(brs_mat[1, x] > brs_mat[-1, x]) / nrow(brs_mat)))

