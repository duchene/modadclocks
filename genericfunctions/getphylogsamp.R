# This function uses a posterior of trees from beast2 to produce a random sample of phylograms from the posterior.

getphylogsamp <- function(postres, N = 100, postN = 10000, write = F, out.file = "phylogsamp.tre", gmc = F, logfile = "sim.log"){

	samp <- sample(1:postN, N)
	if(gmc == F){
	       getRatogB(postres)
	       rats <- read.nexus("ratogs.nex")
	       rats <- rats[samp]
	} else {
	       rats <- as.numeric(read.table(logfile, header = T, comment = "#")[samp, "clockRate"])      
	}

	chrons <- read.nexus(postres)
	chrons <- chrons[samp]
	phylogs <- getpostphylogs(rats, chrons, gmc = gmc)
	
	if(write == T){ write.tree(phylogs, out.file) 
	} else { return(phylogs) }

}