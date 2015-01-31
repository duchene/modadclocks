# This function computes the multinomial likelihood test statistic as proposed by Bollback (2002).

require(phangorn)

multlik <- function(al){

	if(class(al) != "DNAbin"){ al <- as.list(as.DNAbin(al)) } else { al <- as.list(al) }
	mat <- as.character(as.list(as.matrix(al))[[1]])
	for(i in 2:length(as.list(al))){
	      mat <- rbind(mat, as.character(as.list(as.matrix(al))[[i]]))
	}
	al <- mat
	nsites <- ncol(al)
	usites <- unique(al, MARGIN = 2)
	liks <- 0
	for(i in 1:ncol(usites)){
	      insts <- sum(apply(al, 2, identical, usites[, i]))
	      
	      liks <- liks + log((insts / nsites)^insts)
	      #print(liks)
	}

	return(liks)
}
