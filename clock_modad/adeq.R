# This is a 'wrapper' function. It takes all the information that the other functions require to perform clock and substitution model adequacy in one step. The information required includes: The trees log file path, parameters log file path, empirical sequence data file path, sequence length, assumed tree topology, number of PP simulations to take.


adeq <- function(trees.file, log.file, empdat.file, Nsim = 100){
     empdat <- as.phyDat(as.DNAbin(read.nexus.data(empdat.file)))
     seqlen <- ncol(as.matrix(as.DNAbin(empdat)))
     tree.topo <- read.nexus(trees.file)[[1]]
     sims <- make.pps.als(trees.file, log.file, Nsim, seqlen)
     sims <- make.pps.tr(sims, empdat, tree.topo)
     bls <- compile.results(sims)
     return(bls)
}