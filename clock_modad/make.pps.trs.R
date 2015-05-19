# This function processes simulated datasets created with the function make.pps.als. It takes the posterior phylograms and simulated alignments. The function creates the posterior predictive simulated distribution of phylograms. this function also requires the empirical dataset and the assumed topology.

make.pps.tr <- function(sims, empdat, topo){
	 topo$edge.length <- rep(1, length(topo$edge.length))
	 emphy <- optim.pml(pml(topo, empdat))
	 empmlik <- multlik(empdat)
	 empbl <- emphy$tree$edge.length
	 simbl <- emphy$tree$edge.length
	 simultlik <- vector()
	 for(i in 1:length(sims)){
	       ppsphy <- optim.pml(pml(topo, sims[[i]][[3]]))
	       simbl <- rbind(simbl, ppsphy$tree$edge.length)
	       simultlik[i] <- multlik(sims[[i]][[3]])
	       print(paste("Simulation", i, "processed"))
	 }
	 simbl <- simbl[2:nrow(simbl),]
	 return(list(empbl, simbl, empmlik, simultlik))

}