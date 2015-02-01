# This function processes simulated datasets created with the function simPhylogs.R by creating folders and xml files for each. The correct XML maker function must be loaded.

prepsims <- function(sims, cals = NULL){

	 for(i in 1:length(sims)){
	       system(paste0("mkdir ppsim", i))
	       setwd(paste0("ppsim", i))
	       write.tree(sims[[i]][[1]], paste0("ppsphylog", i, ".tre"))
	       write.csv(sims[[i]][[2]], paste0("ppstrdatmat", i, ".csv"))
	       write.nexus.data(as.list(as.DNAbin(sims[[i]][[3]])), paste0("ppsdat", i, ".nex"))
	       dat <- readLines(paste0("ppsdat", i, ".nex"))
	       dat <- gsub("DATATYPED", "DATATYPE=D", dat)
	       writeLines(dat, paste0("ppsdat", i, ".nex"))
	       get.xml2(cals, list(sims[[i]][[1]], as.matrix(as.DNAbin(sims[[i]][[3]]))), "pps")
	       setwd("..")
	 }

}