# This function produces phylograms from a list of chronograms and a list of ratograms.

getpostphylogs <- function(ratogs, chronos){
	phylogs <- list()
	for(i in 1:length(ratogs)){

	      tr <- ratogs[[i]]
	      tr$edge.length <- ratogs[[i]]$edge.length * chronos[[i]]$edge.length
	      phylogs[[i]] <- tr
	}

	return(phylogs)

}