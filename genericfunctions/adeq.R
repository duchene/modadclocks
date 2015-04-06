# This is a 'wrapper' function. It takes all the information that the other functions require to perform clock and substitution model adequacy in one step. The information required includes: The trees log file path, parameters log file path, empirical sequence data file path, sequence length, assumed tree topology, number of PP simulations to take.


adeq <- function(trees.file, log.file, empdat.file, seqlen = 1000, tree.topo, Nsim = 100){
     sims <- make.pps.als(trees.file, log.file, Nsim, seqlen)
     sims <- make.pps.tr(sims, empdat.file, tree.topo)
     bls <- compile.results(sims)
     return(bls)
}