##### Set Up the Environment #####
library(data.table)

# Have the user select the relevant files.
MSIDataPath = choose.files(multi = FALSE, caption = "Select MSI trinuc reads.",
                           filters = c(c("Bed File (*.bed)","Any files"),c("*.bed","*.*")), index = 1)
parentMutationDataPath = choose.files(multi = FALSE, caption = "Select parent mutation data.",
                                      filters = c(c("Bed File (*.bed)","Any files"),c("*.bed","*.*")), index = 1)

# The mutation types to be examined.
mutationTypes = c("C>A","C>G","C>T","T>A","T>C","T>G")

# Read in the data, and extract the mutation type (e.g. C>T)
MSIData = fread(file = MSIDataPath)
MSIMutations = MSIData[,5]
MSIMutationFrequencies = as.data.table(table(MSIMutations))
setkey(MSIMutationFrequencies,"MSIMutations")

parentMutationData = fread(file = parentMutationDataPath)
parentMutations = parentMutationData[,5]

# This function will generate the random sample of parent mutations equal equal in size to the number of
# MSIMutations.  Then, the differences between the mutation counts for the two data sets will be returned.
createComparisonRow = function(i,MSIMutationFrequencies, parentMutations) {
  
  print(paste0("Creating row ",i))
  thisSample = parentMutations[sample(.N,sum(MSIMutationFrequencies[,2]))]
  sampleFrequencies = as.data.table(table(thisSample))
  return(as.list(MSIMutationFrequencies[,N]-sampleFrequencies[,N]))
  
}

# Generate the comparison table
comparisonTable = data.table(t(sapply(1:1000, function(i) createComparisonRow(i,MSIMutationFrequencies,parentMutations))))
colnames(comparisonTable) = gsub(">","_To_",mutationTypes)

# Generates statistics for a given mutation in the comparison table to make a row in the results table.
createResultRow = function(mutation, comparisonTable) {
  mutationSubset = unlist(comparisonTable[,gsub(">","_To_",mutation), with = FALSE])
  return(list(mutation,mean(mutationSubset),median(mutationSubset),sqrt(var(mutationSubset))))
}

# Generate the results table
resultsTable = data.table(t(sapply(mutationTypes, function(x) createResultRow(x,comparisonTable))))
colnames(resultsTable) = c("Mutation","Mean","Median","Std_Dev")
