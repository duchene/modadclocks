
library(phangorn)
library(ClockstaRG)
library(rgl)

# Load true data and obtain a phylogram
true_data <- as.DNAbin(read.nexus.data('al.nex'))
true_phylo <- optim.pml(pml(nj(dist.dna(true_data, model = 'JC')), phyDat(true_data)), optNni = T, optEdge = T)$tree


# Obtain phylograms for each PPS
pps_files <- system('ls -d negcontr_gmc/pps/*', intern = T)
pps_files <- pps_files[-(grep('[.]', pps_files))]

pps_phylogs <- list()

for(i in 1:length(pps_files)){
    cat('reading data', pps_files[i], '\n')
    f_name <- grep('nex', dir(pps_files[i]), value = T)
    pps_data <- as.DNAbin(read.nexus.data(paste0(pps_files[i], '/', f_name)))
    # Optimise branch lengths on the true topology
    pps_phylogs[[i]] <- optim.pml(pml(true_phylo, data = phyDat(pps_data)), model = 'JC', optNni = F, optEdge = T)$tree
#    pps_phylogs[[i]] <- nj(dist.dna(pps_data, model = 'JC'))
}

# Append the 'true' phylogram at the top of the list and obtain the branch length matirx


get_brs_mat <- function(tree_list){
  br_lens_matrix <- matrix(NA, nrow = length(tree_list), ncol = length(tree_list[[1]]$edge.length))
  rownames(br_lens_matrix) <- names(tree_list)
  colnames(br_lens_matrix) <- paste0('br', 1:length(tree_list[[1]]$edge.length))
  for(i in 1:length(tree_list)){
#    br_lens_matrix[i, ] <- ladderize(tree_list[[i]])$edge.length
      br_lens_matrix[i, ] <- tree_list[[i]]$edge.length
  }
  return(br_lens_matrix)
}

all_trees <- list()
all_trees[[1]] <- true_phylo
for(i in 1:length(pps_phylogs)){
    all_trees[[i + 1]] <- pps_phylogs[[i]]
}

brs_mat <- get_brs_mat(all_trees)

# get quantiles per row

row_quants <- sapply(1:ncol(brs_mat), function(x) abs(0.5 - sum(brs_mat[1, x] > brs_mat[-1, x])/nrow(brs_mat)))

