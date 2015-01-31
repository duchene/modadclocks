
get.xml2 <- function(cal.data, sim.data, file.name = "sim"){

# BLOCK1 - Header

b1 <- "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?><beast beautitemplate=\'Standard\' beautistatus=\'\' namespace=\"beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood\" version=\"2.0\">
"

# BLOCK2 - Alignment

b2 <- c("<data", paste0("id=\"simulationduchene\""), "name=\"alignment\">")

for(i in 1:nrow(sim.data[[2]])){
    sec <- paste0("<sequence id=\"seq_", rownames(sim.data[[2]])[i], "\" taxon=\"", rownames(sim.data[[2]])[i], "\" totalcount=\"4\" value=\"", paste(toupper(sim.data[[2]][i, ]), collapse = ""), "\"/>")
    b2 <- c(b2, sec)
}
b2 <- c(b2, "</data>")

# BLOCK3 - Mathematical distributions and tree

b3 <- "<map name=\"Beta\">beast.math.distributions.Beta</map>
<map name=\"Exponential\">beast.math.distributions.Exponential</map>
<map name=\"InverseGamma\">beast.math.distributions.InverseGamma</map>
<map name=\"LogNormal\">beast.math.distributions.LogNormalDistributionModel</map>
<map name=\"Gamma\">beast.math.distributions.Gamma</map>
<map name=\"Uniform\">beast.math.distributions.Uniform</map>
<map name=\"prior\">beast.math.distributions.Prior</map>
<map name=\"LaplaceDistribution\">beast.math.distributions.LaplaceDistribution</map>
<map name=\"OneOnX\">beast.math.distributions.OneOnX</map>
<map name=\"Normal\">beast.math.distributions.Normal</map>
"

tr.topo <- sim.data[[1]]
tr.topo$edge.length <- abs(rnorm(n = length(tr.topo$edge.length)))

b3 <- c(b3, paste0("<tree id=\"Tree.t:simulationduchene\" spec=\'beast.util.TreeParser\' newick=\"", write.tree(tr.topo), "\""), "taxa=\"@simulationduchene\"/>")

# BLOCK4 - State and population model

b4 <- "<run chainLength=\"1000000\" id=\"mcmc\" spec=\"MCMC\">
    <state id=\"state\" storeEvery=\"100\">
        <input name=\'stateNode\' idref=\'Tree.t:simulationduchene\'/>
        <parameter id=\"birthRate.t:simulationduchene\" name=\"stateNode\">1.0</parameter>
        <parameter id=\"clockRate.c:simulationduchene\" name=\"stateNode\">1.0</parameter>
    </state>
"

# BLOCK5 - Prior 1

b5 <- "<distribution id=\"posterior\" spec=\"util.CompoundDistribution\">
        <distribution id=\"prior\" spec=\"util.CompoundDistribution\">
            <distribution birthDiffRate=\"@birthRate.t:simulationduchene\" id=\"YuleModel.t:simulationduchene\" spec=\"beast.evolution.speciation.YuleModel\" tree=\"@Tree.t:simulationduchene\"/>
            <prior id=\"YuleBirthRatePrior.t:simulationduchene\" name=\"distribution\" x=\"@birthRate.t:simulationduchene\">
                <Uniform id=\"Uniform.0\" name=\"distr\" upper=\"Infinity\"/>
            </prior>
            <prior id=\"ClockPrior.c:simulationduchene\" name=\"distribution\" x=\"@clockRate.c:simulationduchene\">
	    	   <Uniform id=\"Uniform.0.1\" name=\"distr\" upper=\"Infinity\"/>
	    </prior>
"

# BLOCK6 - Prior 2 - Calibrations

b6 <- vector()

for(i in 1:length(cal.data)){
    cal.name <- paste0("cal", i)
    cal.block1 <- c(paste0("<distribution id=\"", cal.name, ".prior\" spec=\"beast.math.distributions.MRCAPrior\" tree=\"@Tree.t:simulationduchene\">"), paste0("<taxonset id=\"", cal.name, "\" spec=\"TaxonSet\">")) 
    prevcals <- vector()
    if(i > 1){
	    for(j in 1:(i-1)){
    		prevcals <- c(prevcals, cal.data[[j]][[2]])
    	}
    }
    cal.block2 <- vector()
    for(j in 1:length(cal.data[[i]][[2]])){
    	if(cal.data[[i]][[2]][j] %in% prevcals){
    		cal.block2 <- c(cal.block2, paste0("<taxon idref=\"", cal.data[[i]][[2]][j], "\"/>"))
    	} else {
    		cal.block2 <- c(cal.block2, paste0("<taxon id=\"", cal.data[[i]][[2]][j], "\" spec=\"Taxon\"/>"))
    	}
    }
    cal.block3 <- c("</taxonset>", paste0("<Normal id=\"Normal.0", i, "\" name=\"distr\">"), paste0("<parameter estimate=\"true\" id=\"RealParameter.a", i, "\" name=\"mean\">", round(cal.data[[i]][[1]], 3), "</parameter>"), paste0("<parameter estimate=\"true\" id=\"RealParameter.b", i, "\" name=\"sigma\">", round((cal.data[[i]][[1]] / 10), 5), "</parameter>"), "</Normal>", "</distribution>")
    b6 <- c(b6, cal.block1, cal.block2, cal.block3)
}

b6 <- c(b6, "</distribution>")

# BLOCK7

b7 <- "<distribution id=\"likelihood\" spec=\"util.CompoundDistribution\">
            <distribution data=\"@simulationduchene\" id=\"treeLikelihood.simulationduchene\" spec=\"TreeLikelihood\" tree=\"@Tree.t:simulationduchene\">
                <siteModel id=\"SiteModel.s:simulationduchene\" spec=\"SiteModel\">
                    <parameter estimate=\"true\" id=\"mutationRate.s:simulationduchene\" name=\"mutationRate\">1.0</parameter>
                    <substModel id=\"JC69.s:simulationduchene\" spec=\"JukesCantor\"/>
                </siteModel>
                <branchRateModel clock.rate=\"@clockRate.c:simulationduchene\" spec=\"beast.evolution.branchratemodel.StrictClockModel\">
            </distribution>
        </distribution>
    </distribution>
"

# BLOCK8 - Operators

b8 <- "<operator id=\"YuleBirthRateScaler.t:simulationduchene\" parameter=\"@birthRate.t:simulationduchene\" scaleFactor=\"0.75\" spec=\"ScaleOperator\" weight=\"3.0\"/>

    <operator id=\"treeScaler.t:simulationduchene\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" tree=\"@Tree.t:simulationduchene\" weight=\"3.0\"/>

    <operator id=\"treeRootScaler.t:simulationduchene\" rootOnly=\"true\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" tree=\"@Tree.t:simulationduchene\" weight=\"3.0\"/>

    <operator id=\"UniformOperator.t:simulationduchene\" spec=\"Uniform\" tree=\"@Tree.t:simulationduchene\" weight=\"30.0\"/>
    
    <operator id=\"strictClockRateScaler.c:simulationduchene\" parameter=\"@clockRate.c:simulationduchene\" scaleFactor=\"0.7\" spec=\"ScaleOperator\" weight=\"3.0\"/>

    <operator id=\"strictClockUpDownOperator.c:simulationduchene\" scaleFactor=\"0.75\" spec=\"UpDownOperator\" weight=\"3.0\">
    	      <parameter idref=\"clockRate.c:simulationduchene\" name=\"up\"/>
	      <tree idref=\"Tree.t:simulationduchene\" name=\"down\"/>
    </operator>
"

# BLOCK 9 - Loggers

b9 <- paste0("<logger fileName=\"", file.name, ".log\" id=\"tracelog\" logEvery=\"100\" model=\"@posterior\" sanitiseHeaders=\"true\" sort=\"smart\">")

b9 <- c(b9, "<log idref=\"posterior\"/>
        <log idref=\"likelihood\"/>
        <log idref=\"prior\"/>
        <log idref=\"treeLikelihood.simulationduchene\"/>
        <log id=\"TreeHeight.t:simulationduchene\" spec=\"beast.evolution.tree.TreeHeightLogger\" tree=\"@Tree.t:simulationduchene\"/>
        <log idref=\"YuleModel.t:simulationduchene\"/>
        <parameter idref=\"birthRate.t:simulationduchene\" name=\"log\"/>
        <parameter idref=\"clockRate.c:simulationduchene" name=\"log\"/>
")

cals <- vector()
for(i in 1:length(cal.data)){
	cals <- c(cals, paste0("<log idref=\"cal", i, ".prior\"/>"))
}

b9 <- c(b9, cals, "</logger>")

b10 <- "<logger id=\"screenlog\" logEvery=\"100\">
        <log idref=\"posterior\"/>
        <log arg=\"@posterior\" id=\"ESS.0\" spec=\"util.ESS\"/>
        <log idref=\"likelihood\"/>
        <log idref=\"prior\"/>
    </logger>

    <logger fileName=\"$(tree).trees\" id=\"treelog.t:simulationduchene\" logEvery=\"100\" mode=\"tree\">
        <log branchratemodel=\"@clockRate.c:simulationduchene\" id=\"TreeWithMetaDataLogger.t:simulationduchene\" spec=\"beast.evolution.tree.TreeWithMetaDataLogger\" tree=\"@Tree.t:simulationduchene\"/>
    </logger>

</run>

</beast>
"

complete.xml <- c(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10)

writeLines(complete.xml, con = paste0(file.name, ".xml"))

}
