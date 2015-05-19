# This function takes posterior log and trees files. It simulates alignments using the model parameters from the corresponding line in the log file. This function also requires the length of the alignment.

make.pps.als <- function(trees.file, log.file, N = 100, l = 1500, savetrees = F){
	
	trees <- read.nexus(trees.file)
	logdat <- read.table(log.file, header = T, comment = "[")
	endburnin <- round((length(trees) * 0.1), 0)
	samp <- sample(endburnin:length(trees), N)
	trees <- trees[samp]
	logdat <- logdat[samp,]
	if(savetrees){
	      write.tree(trees, file = paste0(trees.file, N, ".tre"))
	      write.csv(logdat, file = paste0(trees.file, N, ".csv"))
	}
	sim <- list()
	for(i in 1:nrow(logdat)){
	      sim[[i]] <- list(phylogram = trees[[i]])
	      
	      if(all(c("r.A...C.", "r.A...G.", "r.A...T.", "r.C...G.", "r.C...T.", "r.G...T.") %in% colnames(logdat))){
	      	     # GENERAL TIME REVERSIBLE (GTR)
		     #print("The substitution model is GTR")
		     basef <- c(logdat$pi.A.[i], logdat$pi.C.[i], logdat$pi.G.[i], logdat$pi.T.[i])
		     qmat <- c(logdat$r.A...C[i], logdat$r.A...G.[i], logdat$r.A...T.[i], logdat$r.C...G.[i], logdat$r.C...T.[i], logdat$r.G...T.[i])
		     #print(basef)
		     #print(qmat)

		     if("alpha" %in% colnames(logdat)){
		     	   rates = phangorn:::discrete.gamma(logdat$alpha[i], k = 4)
        		   sim_dat_all<- lapply(rates, function(r) simSeq(sim[[i]][[1]], l = round(l/4, 0), Q = qmat, bf = basef, rate = r))
        		   sim[[i]][[3]] <- c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]])
		     } else {		    
	      	       	   sim[[i]][[3]] <- simSeq(sim[[i]][[1]], Q = qmat, bf = basef, l = l)
		     }		     
		     
	      } else if("kappa" %in% colnames(logdat)){
		     # HASEGAWA-KISHINO-YANO (HKY)
		     #print("The substitution model is HKY")
		     basef <- c(logdat$pi.A.[i], logdat$pi.C.[i], logdat$pi.G.[i], logdat$pi.T.[i])
		     qmat <- c(1, 2*logdat$kappa[i], 1, 1, 2*logdat$kappa[i], 1)

		     if("alpha" %in% colnames(logdat)){
		     	   rates = phangorn:::discrete.gamma(logdat$alpha[i], k = 4)
       			   sim_dat_all<- lapply(rates, function(r) simSeq(sim[[i]][[1]], l = round(l/4, 0), Q = qmat, bf = basef, rate = r))
                           sim[[i]][[3]] <- c(sim_dat_all[[1]], sim_dat_all[[2]], sim_dat_all[[3]], sim_dat_all[[4]])
		     } else {
                           sim[[i]][[3]] <- simSeq(sim[[i]][[1]], Q = qmat, bf = basef, l = l)
	       	     }
		     
	      } else { 
	      	     # JUKES-CANTOR (JC)
		     #print("The substitution model is assumed to be JC")
	      	     sim[[i]][[3]] <- simSeq(sim[[i]][[1]], l = l)
	      }

	}
	
	return(sim)

}