# This function produces phylograms from a list of chronograms and a list of ratograms.

getpostphylogs <- function(ratogs, chronos, gmc = F){
	phylogs <- list()
	
	if(gmc == F){
	for(i in 1:length(ratogs)){

	      tr <- ratogs[[i]]
	      tr$edge.length <- ratogs[[i]]$edge.length * chronos[[i]]$edge.length
	      phylogs[[i]] <- tr
	}
	} else {
	for(i in 1:length(chronos)){
	      tr <- chronos[[i]]
	      tr$edge.length <- chronos[[i]]$edge.length * ratogs[i]
	      phylogs[[i]] <- tr
	}
	}
	class(phylogs) <- "multiPhylo"
	return(phylogs)

}