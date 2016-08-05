source('R/sciClone.R')
source('R/clustering.R')
source('R/object.R')

binomCluster <- function(){
  
  outputDir = "data/Cov_1000_Samples_4_Mut_100_Clone_10_PCR_Removed/sim_0/"
  numSamples = 4
    
  names = vector(mode="character", numSamples)
  vafsList = vector("list", numSamples)
  
  for (i in 1:numSamples){
    names[i] <- paste("Sample", i, sep=" ")
    # Add psuedo-counts
    vafsTable = read.table(paste(outputDir, "VAF_file_", i, ".vaf", sep=""))
    vafsTable[, 4] <- vafsTable[, 4] + 1e-12
    vafsList[[i]] <- vafsTable
  }
  
  sc = sciClone(vafs = vafsList,
                sampleNames = names,
                clusterMethod = "binomial.bmm",
                minimumDepth = 0,
                maximumClusters=50,
                doClusteringAlongMargins = FALSE)
  
  writeClusterTable(sc, outputFile = paste(outputDir, "binom_clusters", sep=""))
  
  write.table(sc@clust$a, file = paste(outputDir, "binom_params_as", sep=""), sep="\t")
  write.table(sc@clust$b, file = paste(outputDir, "binom_params_bs", sep=""), sep="\t")
}