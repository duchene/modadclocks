# This function takes a chronogram and a phylogram, re-roots the phylogram appropriately, and returns both ladderized.

rootnlad <- function(chron, phyl){

	 rootnode <- length(chron$tip.label) + 1
	 outgr <- chron$tip.label[Descendants(chron, Descendants(chron, rootnode, "children")[1], "tips")[[1]]]
	 chron <- ladderize(root(chron, outgroup = outgr))
	 phyl <- ladderize(root(phyl, outgroup = outgr))
	 
	 return(list(chron, phyl))
}

# This function processes simulated datasets created with the function simPhylogsBL.R by taking posterior phylograms and simulated alignmetns, and creating posterior predictive simulated phylograms for each.

prepsimsBL <- function(sims, empdat, chron, writedat = F){
	 empdat <- as.phyDat(as.DNAbin(read.nexus.data(empdat)))
	 topo <- read.tree(chron)
	 topo$edge.length <- rep(1, length(topo$edge.length))
	 emphy <- optim.pml(pml(topo, empdat))
	 empmlik <- multlik(empdat)
	 empbl <- emphy$tree$edge.length
	 simbl <- emphy$tree$edge.length
	 simultlik <- vector()
	 for(i in 1:length(sims)){
	       #print(class(sims[[i]][[3]]))	       
	       ppsphy <- optim.pml(pml(topo, sims[[i]][[3]]))
	       #trs <- rootnlad(sims[[i]][[1]], ppsphy$tree)
	       #sims[[i]][[1]] <- trs[[1]]
	       #sims[[i]][[2]] <- trs[[2]]	       
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