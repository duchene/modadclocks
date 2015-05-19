# This is a 'wrapper' function. It takes all the information that the other functions require to perform substitution model adequacy in one step. The information required includes: The trees log file path, parameters log file path, empirical sequence data file path, sequence length, assumed tree topology, number of PP simulations to take.


adeq <- function(trees.file, log.file, empdat.file, true.topo, Nsim = 100, savetrees = F){
     empdat <- as.phyDat(as.DNAbin(read.nexus.data(empdat.file)))
     seqlen <- ncol(as.matrix(as.DNAbin(empdat)))
     sims <- make.pps.als(trees.file, log.file, Nsim, seqlen, savetrees)
     if(missing(true.topo)){
	sims <- make.pps.tr(sims, empdat)
     } else {
     	sims <- make.pps.tr(sims, empdat, true.topo)
     }     
     return(sims)
}