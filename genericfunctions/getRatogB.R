# This function is to get ratograms from the output of beast2. It takes an input file, and returns an R object or trees.

require(phangorn)

getRatogB <- function(trees, out.file = "ratogs.tree"){

	  trs <- readLines(trees)

	  trs <- gsub("[[]&([a-z])+[=]", ":", trs)
	  trs <- gsub("[]]:([0-9]|[.])+", "", trs)

	  writeLines(trs, out.file)
}