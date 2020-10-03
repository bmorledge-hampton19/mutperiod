# Given a path to a file of raw mutation counts with a corresponding file of background counts,
# Read in the data and produce normalized counts, writing them to a new file if specified (by default true).
#' @export
normalizeNucleosomeMutationCounts = function(rawCountsFilePath, backgroundCountsFilePath,
                                             normalizedDataFilePath = NA, writeNormalizedData = TRUE) {

  # Read in the data
  rawCounts = data.table::fread(file = rawCountsFilePath)
  backgroundCounts = data.table::fread(file = backgroundCountsFilePath)

  # Make sure that the raw and background counts both use the same number of dyad positions.

  if (dim(rawCounts)[1] != dim(backgroundCounts)[1]) {
    stop("Unequal dyad positions in raw vs. background counts data.")
  }

  # Create a table of normalized values from the given data
  normalizedData = data.table::data.table()
  normalizedData[,Dyad_Position := backgroundCounts$Dyad_Position]
  normalizedData[,Normalized_Minus_Strand := mapply(normalize,rawCounts$Minus_Strand_Counts,
                                                    backgroundCounts$Expected_Mutations_Minus_Strand)]
  normalizedData[,Normalized_Plus_Strand := mapply(normalize,rawCounts$Plus_Strand_Counts,
                                                   backgroundCounts$Expected_Mutations_Plus_Strand)]
  normalizedData[,Normalized_Both_Strands := mapply(normalize,rawCounts$Both_Strands_Counts,
                                                    backgroundCounts$Expected_Mutations_Both_Strands)]

  # Add a column for aligned, normalized strands
  normalizedData[,Normalized_Aligned_Strands :=
                   mapply(function(plus,minus) mean(c(plus,minus)),
                          Normalized_Plus_Strand, rev(Normalized_Minus_Strand))]

  # Write the normalized data to a new file. (If desired)
  if (writeNormalizedData) {
    data.table::fwrite(normalizedData, sep = '\t', file = normalizedDataFilePath)
  }

  return(normalizedData)

}

# Normalizes data by dividing raw by expected, except in the case where expected is 0.
# In this case, raw will also be 0, or at least some mutations would be expected, so
# the returned value can safely be set to 0 instead of dividing by 0.
normalize = function(raw, expected) {
  if (expected == 0) return(0)
  else return(raw/expected)
}
