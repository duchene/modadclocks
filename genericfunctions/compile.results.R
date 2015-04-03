# This function takes the results of any given run and returns the pps branch lengths, the pvalue for each brnach length and prints the global pvalue and the mean range in banch lengths.

require(phangorn)

compile.results <- function(branch.results){
	
	bls <- branch.results
	pvals <- vector()
	for(i in 1:length(bls[[1]])){
	      pvals[i] <- length(bls[[2]][,i][which(bls[[2]][,i] < bls[[1]][i])]) / nrow(bls[[2]])
	}
	
	fullpval <- 1 - (length(pvals[which(pvals < 0.05 | pvals > 0.95)]) / length(pvals))
	print(fullpval)
	meandifran1 <- mean(sapply(1:ncol(bls[[2]]), function(x) diff(range(bls[[2]][,x])) / mean(bls[[2]][,x])))
	meandifran2 <- mean(sapply(1:ncol(bls[[2]]), function(x) diff(quantile(bls[[2]][,x], prob = c(0.05, 0.95))) / mean(bls[[2]][,x])))
	bls[[5]] <- pvals
	bls[[6]] <- length(bls[[4]][which(bls[[4]] < bls[[3]])]) / length(bls[[4]])
	bls[[7]] <- meandifran1
	bls[[8]] <- meandifran2
	
	names(bls) <- c("empirical_branch_lengths", "pps_branch_lengths", "empirical_multlik", "pps_multliks", "branch_length_pvalues", "mlik_pvalues", "ppsbl_range_over_mean", "ppsbl_range_CI")
	return(bls)

}