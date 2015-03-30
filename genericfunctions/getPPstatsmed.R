# This fuunction produces test statsitics for each analysis by searching for a posterior distribution of empirical data, and then for a pps folder that should have the pp simulatons and runs.

require(phangorn)
require(NELSI)

getPPstatsmed <- function(N = 200){
	if(length(grep("[.]log", dir(), value = T)) > 1 || length(grep("[.]trees", dir(), value = T)) > 1 || length(grep("[.]nex", dir(), value = T)) > 1){ 
		print("Theres more than one empirical dataset or posterior"); stop
	}
	empost <- read.table(grep("[.]log", dir(), value = T), header = T, comment = "#", sep = ",")
	samp <- sample(1:nrow(empost), N)
	empost <- empost[samp,]
	emptre <- read.nexus(grep("[.]trees", dir(), value = T))
	emptre <- emptre[samp]
	empal <- read.nexus.data(grep("[.]nex", dir(), value = T))
	#### Test stats for empirical data: empmlik, emplogstats
	####

	empmlik <- multlik(empal)
	emplogstats <- rbind(apply(empost, 2, median), apply(empost, 2, var))
	colnames(emplogstats) <- colnames(empost)
	if("clockRate" %in% colnames(empost)){ 
	       phylogs <- getpostphylogs(empost[,"clockRate"], emptre, gmc = T)
	       print(head(phylogs[[2]]$edge.length))
	       emptrlens <- sapply(phylogs, function(x) sum(x$edge.length))
	} else { 
	       #getRatogB(grep("[.]trees", dir(), value = T))
	       ratos <- read.nexus("ratogs.tree")[samp]
	       print(head(ratos[[2]]$edge.length))
	       phylogs <- getpostphylogs(ratos, emptre)
	       print(head(phylogs[[2]]$edge.length))
	       emptrlens <- sapply(phylogs, function(x) sum(x$edge.length))
	}
	emplogstats <- cbind(emplogstats, treeLength = c(median(emptrlens), var(emptrlens)))
	print(empmlik); print(emplogstats)    
	####
	####
	setwd("pps/")
	ppsfiles <- grep("pps", dir(), value = T)
	ppsmlik <- vector()
	ppslogstats <- list()
	for(i in 1:length(ppsfiles)){
	      #print(getwd())
	      #print(ppsfiles[i])
	      setwd(ppsfiles[i])
	      ppspost <- try(read.table(grep("[.]log", dir(), value = T), header = T, comment = "#", sep = ","))
	      if(class(ppspost) == "try-error"){ setwd(".."); next }
	      samp <- sample(1:nrow(ppspost), N)
	      ppspost <- ppspost[samp,]
	      ppstre <- try(read.nexus(grep("[.]trees", dir(), value = T)))
	      if(class(ppstre) == "try-error"){ setwd(".."); next }
	      ppstre <- ppstre[samp]
	      ppsal <- read.nexus.data(grep("[.]nex", dir(), value = T))
	      ##### Test stats for PPS data: ppsmlik, ppslogstats
	      #####

	      ppsmlik[i] <- multlik(ppsal)
              ppslogstats[[i]] <- rbind(apply(ppspost, 2, median), apply(ppspost, 2, var))
              colnames(ppslogstats[[i]]) <- colnames(ppspost)
              if("clockRate" %in% colnames(ppspost)){
              		phylogs <- getpostphylogs(ppspost[,"clockRate"], ppstre, gmc = T)
              		ppstrlens <- sapply(phylogs, function(x) sum(x$edge.length))
               } else {
			#getRatogB(grep("[.]trees", dir(), value = T))
               		ratos <- read.nexus("ratogs.tree")[samp]
               		phylogs <- getpostphylogs(ratos, ppstre)
			#print(mean(ratos[[2]]$edge.length))
			#print(mean(phylogs[[2]]$edge.length))
               		ppstrlens <- sapply(phylogs, function(x) sum(x$edge.length))
              }
	      
              ppslogstats[[i]] <- cbind(ppslogstats[[i]], treeLength = c(median(ppstrlens), var(ppstrlens)))
	      #print(ppsmlik[i]); print(ppslogstats[[i]])
	      ####
	      ####
	      setwd("..")
	      print(paste("PPS", i, "processed."))
	}
	setwd("..")
	allstats <- list(empirical = list(mult.lik = empmlik, post.stats = emplogstats), pps = list(mult.lik = ppsmlik, post.stats = ppslogstats))
	return(allstats)
	
}