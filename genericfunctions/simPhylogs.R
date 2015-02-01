# This function takes posterior log and trees files and simulates rates thorugh each of the posterior trees using the model parameters from the corresponding line in the log file.

simPhylogs <- function(trees, logdat, N = 100, clockmod = "gmc", subsmod = "JC", l = 1000, logsize = 10000){
	samp <- sample(1:logsize, N)
	trees <- read.nexus(trees)[samp]
	logdat <- read.table(logdat, header = T, comment = "#")[samp,]
	sim <- list()
	for(i in 1:nrow(logdat)){
	      tr <- trees[[i]]
	      if(clockmod == "gmc"){
	      	     sim[[i]] <- simulate.clock(tr, params = list(rate = logdat[i, "clockRate"], noise = 0.000001))
	      } else if(clockmod == "ucln"){
	      	     sim[[i]] <- simulate.uncor.lnorm(tr, params = list(mean.log = log(logdat[i, "ucldMean"]), sd.log = logdat[i, "ucldStdev"]))
	      }
	      
	      if(subsmod == "JC"){
	      	     sim[[i]][[3]] <- simSeq(sim[[i]][[1]], l = l)
	      } else if(subsmod == "GTR"){
	      	     
	      	     sim[[i]][[3]] <- simSeq(sim[[i]][[1]], Q = , bf = , l = l)
	      }

	}
	
	return(sim)

}