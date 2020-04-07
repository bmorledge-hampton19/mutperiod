# Given Dr. Wyrick's nucleosome mutation data and background mutation data, 
# returns a data.table of base and normalized mutation data with respect to dyad position.
library(data.table)

parseWyrickNucleosomeMutationData = function(mutationDataFilePath, backgroundDataFilePath) {
  
  # Parse data on total mutation counts by dyad position.
  rawData = readLines(mutationDataFilePath)
  relevantRows = strsplit(rawData[c(3,10,13,17)],"\t")
  nucleosomeMutationData = data.table(sapply(relevantRows, unlist))
  setnames(nucleosomeMutationData,colnames(nucleosomeMutationData),
           gsub(" ", "_", unlist(nucleosomeMutationData[1],use.names = FALSE)))
  nucleosomeMutationData = nucleosomeMutationData[-1]
  nucleosomeMutationData = nucleosomeMutationData[,lapply(.SD, as.numeric)]
  setnames(nucleosomeMutationData,1,"Dyad_Position")
  
  # Parse data on the mutational background to be used in normalization
  rawData = readLines(backgroundDataFilePath)
  relevantRows = strsplit(rawData[c(2,4,5)], "\t")
  expectedMutations = data.table(sapply(relevantRows, unlist))
  setnames(expectedMutations,colnames(expectedMutations),
           gsub(" ", "_", unlist(expectedMutations[1],use.names = FALSE)))
  expectedMutations = expectedMutations[-1]
  expectedMutations = expectedMutations[,lapply(.SD, as.numeric)]
  expectedMutations[,Total := Plus_Strand + Minus_Strand]
  
  # Normalize the data
  normalizedData = copy(nucleosomeMutationData)
  normalizedData[,Normalized_Plus_Strand := CPDs_Plus_Strand/expectedMutations$Plus_Strand]
  normalizedData[,Normalized_Minus_Strand := CPDs_Minus_Strand/expectedMutations$Minus_Strand]
  normalizedData[,Normalized_Both_Strands := CPDs_both_strands/expectedMutations$Total]
  
  # Add a column for aligned, normalized strands
  normalizedData[,Normalized_Aligned_Strands := 
                   mapply(function(plus,minus) mean(c(plus,minus)), 
                          Normalized_Plus_Strand, rev(Normalized_Minus_Strand))]
  
  return(normalizedData)
  
}


# Given my (Ben's) nucleosome mutation data and background mutation data, 
# returns a data.table of base and normalized mutation data with respect to dyad position.
parseBMHNucleosomeMutationData = function(mutationDataFilePath, backgroundDataFilePath){
 
  # Read in the data
  nucleosomeMutationData = fread(file = mutationDataFilePath)
  expectedMutations = fread(file = backgroundDataFilePath)
  
  # Create a table of normalized values from the given data
  normalizedData = data.table()
  normalizedData[,Dyad_Position := expectedMutations$Dyad_Position]
  normalizedData[,Normalized_Minus_Strand := nucleosomeMutationData$Minus_Strand_Counts/
                   expectedMutations$Expected_Mutations_Minus_Strand]
  normalizedData[,Normalized_Plus_Strand := nucleosomeMutationData$Plus_Strand_Counts/
                   expectedMutations$Expected_Mutations_Plus_Strand]
  normalizedData[,Normalized_Both_Strands := nucleosomeMutationData$Both_Strands_Counts/
                   expectedMutations$Expected_Mutations_Both_Strands]
  
  # Add a column for aligned, normalized strands
  normalizedData[,Normalized_Aligned_Strands := 
                   mapply(function(plus,minus) mean(c(plus,minus)), 
                          Normalized_Plus_Strand, rev(Normalized_Minus_Strand))]
  
  return(normalizedData)
   
}