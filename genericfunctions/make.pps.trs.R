

# This function processes simulated datasets created with the function makeppsims.R by taking posterior phylograms and simulated alignments, and creating posterior predictive simulated phylograms for each. this function also requires the empirical dataset and the assumed topology.

make.pps.tr <- function(sims, empdat.file, topo, writedat = F){
	 empdat <- as.phyDat(as.DNAbin(read.nexus.data(empdat.file)))
	 topo$edge.length <- rep(1, length(topo$edge.length))
	 emphy <- optim.pml(pml(topo, empdat))
	 empmlik <- multlik(empdat)
	 empbl <- emphy$tree$edge.length
	 simbl <- emphy$tree$edge.length
	 simultlik <- vector()
	 for(i in 1:length(sims)){
	       ppsphy <- optim.pml(pml(topo, sims[[i]][[3]]))	       
	       if(writedat == T){
	       write.tree(sims[[i]][[1]], paste0("postphylog", i, ".tre"))
	       write.csv(sims[[i]][[2]], paste0("ppsphylog", i, ".tre"))
	       write.nexus.data(as.list(as.DNAbin(sims[[i]][[3]])), paste0("ppsdat", i, ".nex"))
	       dat <- readLines(paste0("ppsdat", i, ".nex"))
	       dat <- gsub("DATATYPED", "DATATYPE=D", dat)
	       writeLines(dat, paste0("ppsdat", i, ".nex"))
	       }
	       simbl <- rbind(simbl, ppsphy$tree$edge.length)
	       simultlik[i] <- multlik(sims[[i]][[3]])
	       print(paste("finished sim", i))
	      
	 }
	 simbl <- simbl[2:nrow(simbl),]
	 return(list(empbl, simbl, empmlik, simultlik))

}