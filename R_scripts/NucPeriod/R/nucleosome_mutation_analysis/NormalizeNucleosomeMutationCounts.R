library(data.table)

# Given a path to a file of raw mutation counts with a corresponding file of background counts,
# Read in the data and produce normalized counts, writing them to a new file if specified (by default true).
normalizeNucleosomeMutationCounts = function(rawCountsFilePath, writeNormalizedData = TRUE) {
  
  # Check for the background file corresponding to each raw counts file
  if (!endsWith(rawCountsFilePath,"nucleosome_mutation_counts.tsv")) {
    stop(paste0("Given file does not follow the naming convention for a raw counts file.  ", rawCountsFilePath))
  }
  backgroundCountsFilePath = paste0(unlist(strsplit(rawCountsFilePath,"counts.tsv"))[1],"background.tsv")
  if (!file.exists(backgroundCountsFilePath)) {
    stop(paste0("Background file not found at the expected location: ",backgroundCountsFilePath))
  }
  
  # Read in the data
  rawCounts = fread(file = rawCountsFilePath)
  backgroundCounts = fread(file = backgroundCountsFilePath)
  
  # Create a table of normalized values from the given data
  normalizedData = data.table()
  normalizedData[,Dyad_Position := backgroundCounts$Dyad_Position]
  normalizedData[,Normalized_Minus_Strand := rawCounts$Minus_Strand_Counts/
                   backgroundCounts$Expected_Mutations_Minus_Strand]
  normalizedData[,Normalized_Plus_Strand := rawCounts$Plus_Strand_Counts/
                   backgroundCounts$Expected_Mutations_Plus_Strand]
  normalizedData[,Normalized_Both_Strands := rawCounts$Both_Strands_Counts/
                   backgroundCounts$Expected_Mutations_Both_Strands]
  
  # Add a column for aligned, normalized strands
  normalizedData[,Normalized_Aligned_Strands := 
                   mapply(function(plus,minus) mean(c(plus,minus)), 
                          Normalized_Plus_Strand, rev(Normalized_Minus_Strand))]
  
  # Write the normalized data to a new file. (If desired)
  if (writeNormalizedData) {
    normalizedDataFilePath = paste0(unlist(strsplit(rawCountsFilePath,"mutation_counts.tsv"))[1],
                                    "normalized_mutation_counts.tsv")
    fwrite(normalizedData, sep = '\t', file = normalizedDataFilePath)
  }
  
  return(normalizedData)
  
}