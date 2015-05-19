# This function processes simulated datasets created with the function make.pps.als. It takes the posterior phylograms and simulated alignments. The function creates the posterior predictive simulated distribution of phylograms. This function also requires the empirical dataset and, if desired, the true topology.

make.pps.tr <- function(sims, empdat, truetr){
	 require(MASS)
	 empmlik <- multlik(empdat)
	 if(!missing(truetr)) truelen <- sum(truetr$edge.length)
	 simultlik <- vector()
	 simtrlen <- vector()
	 topodists <- vector()
	 for(i in 1:length(sims)){
	       simultlik[i] <- multlik(sims[[i]][[3]])
	       simtrlen[i] <- sum(sims[[i]][[1]]$edge.length)
	       if(!missing(truetr)) topodists[i] <- dist.topo(truetr, sims[[i]][[1]])
	       print(paste("Simulation", i, "processed"))
	 }
	 #multlikp <- length(simultlik[which(simultlik < empmlik)]) / length(sims)
	 multlikdistr <- fitdistr(simultlik, "normal")
	 multlikp <- 1 - pnorm(empmlik, multlikdistr[[1]][1], multlikdistr[[1]][2])
	 
	 if(!missing(truetr)){
	       trlenp <- length(simtrlen[which(simtrlen < truelen)]) / length(sims)
	       res <- list(empirical.multinomial.likelihood = empmlik, pps.multinomial.likelihood = simultlik, multinomial.likelihood.p.value = multlikp, true.tree.length = truelen, pps.tree.lengths = simtrlen, tree.length.p.value = trlenp, topological.distances = topodists)
	 } else {
	       res <- list(empirical.multinomial.likelihood = empmlik, pps.multinomial.likelihood = simultlik, multinomial.likelihood.p.value = multlikp)
	 }

	 return(res)
}