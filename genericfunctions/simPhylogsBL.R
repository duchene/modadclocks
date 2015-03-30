# This function takes posterior log and trees files and simulates rates thorugh each of the posterior trees using the model parameters from the corresponding line in the log file.

simPhylogsBL <- function(treesf, logdat, N = 100, subsmod = "JC", l = 1000, sample = T, ratogindir = T){
	
	trees <- read.nexus(treesf)
	logdat <- read.table(logdat, header = T, comment = "#", sep = ",")
	samp <- sample(1:length(trees), N)
	trees <- trees[samp]
	logdat <- logdat[samp,]
	sim <- list()
	if("ucldMean" %in% colnames(logdat) | "meanClockRate" %in% colnames(logdat)){
	if(ratogindir == T){
	      ratogs <- read.nexus("ratogs.tree")[samp]
	} else {
	      getRatogB(treesf)
	      ratogs <- read.nexus("ratogs.tree")[samp]
	}
	}
	for(i in 1:nrow(logdat)){
	      tr <- trees[[i]]
	      if("clockRate" %in% colnames(logdat)){
	      	     tr$edge.length <- tr$edge.length * logdat[i, "clockRate"]
		     sim[[i]] <- list(phylogram = tr)
	      } else if("ucldMean" %in% colnames(logdat) | "meanClockRate" %in% colnames(logdat)){
	      	     trr <- ratogs[[i]]
		     trr$edge.length <- trr$edge.length * tr$edge.length
	      	     sim[[i]] <- list(phylogram = trr)
	      }
	      
	      if(subsmod == "JC"){
	      	     sim[[i]][[3]] <- simSeq(sim[[i]][[1]], l = l)
	      } else if(subsmod == "GTR"){
	      	     
	      	     sim[[i]][[3]] <- simSeq(sim[[i]][[1]], Q = , bf = , l = l)
	      }

	}
	
	return(sim)

}