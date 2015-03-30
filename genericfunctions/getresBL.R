# This function takes the results of any given run and returns the pps branch lengths, the pvalue for each brnach length and prints the global pvalue and the mean range in banch lengths.

require(phangorn)

getresBL <- function(postres, postlog, alf, chronf, l = 2000){
	 
	sims <- simPhylogsBL(postres, postlog, l = l)
	bls <- prepsimsBL(sims, alf, chronf)
	pvals <- vector()
	for(i in 1:length(bls[[1]])){
	      pvals[i] <- length(bls[[2]][,i][which(bls[[2]][,i] < bls[[1]][i])]) / nrow(bls[[2]])
	}
	
	fullpval <- 1 - (length(pvals[which(pvals < 0.05 | pvals > 0.95)]) / length(pvals))
	print(fullpval)
	meandifran1 <- mean(sapply(1:ncol(bls[[2]]), function(x) diff(range(bls[[2]][,x]))))
	meandifran2 <- mean(sapply(1:ncol(bls[[2]]), function(x) diff(range(bls[[2]][,x])) / mean(bls[[2]][,x])))
	meandifran3 <- mean(sapply(1:ncol(bls[[2]]), function(x) diff(range(bls[[2]][,x])) / bls[[1]][x]))
	meandifran4 <- mean(sapply(1:ncol(bls[[2]]), function(x) diff(quantile(bls[[2]][,x], prob = c(0.05, 0.95))) / mean(bls[[2]][,x])))
	print(c(meandifran1, meandifran2, meandifran3))
	bls[[5]] <- pvals
	bls[[6]] <- length(bls[[4]][which(bls[[4]] < bls[[3]])]) / length(bls[[4]])
	bls[[7]] <- meandifran1
	bls[[8]] <- meandifran2
	bls[[9]] <- meandifran3
	bls[[10]] <- meandifran4
	## The following is to compare the "simulated emprirical" rates to the posterior rates.
	emprat <- read.table("trdatmat.csv", header = T, sep = ",")[,6]
	if("ratogs.tree" %in% dir()){
		rats <- read.nexus("ratogs.tree")
		ratmat <- matrix(NA, ncol = length(rats[[1]]$edge.length), nrow = length(rats))
		for(i in 1:length(rats)){
                        ratmat[i,] <- rats[[i]]$edge.length
                }
	} else {
	        rats <- read.nexus(postres)
		logf <- read.table(postlog, header = T, sep = ",")
		ratmat <- matrix(NA, ncol = length(rats[[1]]$edge.length), nrow = nrow(logf))
		for(i in 1:nrow(logf)){
			ratmat[i,] <- rep(logf[i,"clockRate"], length(rats[[1]]$edge.length))
		}
	}
	bls[[11]] <- vector()
	for(i in 1:length(emprat)){
              bls[[11]][i] <- length(ratmat[,i][which(ratmat[,i] < emprat[i])]) / nrow(ratmat)
        }
	bls[[12]] <- 1 - (length(bls[[11]][which(bls[[11]] < 0.05 | bls[[11]] > 0.95)]) / length(bls[[11]]))
	names(bls) <- c("empirical_branch_lengths", "pps_branch_lengths", "empirical_multlik", "pps_multliks", "branch_length_pvalues", "mlik_pvalues", "ppsbl_range_size", "ppsbl_range_over_mean", "ppsbl_range_over_emp", "ppsbl_range_CI", "emp_postrat_comp", "emp_postrat_p")
	return(bls)

}