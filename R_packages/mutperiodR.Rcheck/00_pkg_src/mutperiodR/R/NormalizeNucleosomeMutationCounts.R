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

  # If the background file is "custom input" it will have the same column names as a raw counts file.
  # Check for the raw counts headers, and replace them with the expected headers if necessary.
  if (colnames(backgroundCounts)[2] == "Plus_Strand_Counts") {
    colnames(backgroundCounts) = c("Dyad_Position", "Expected_Mutations_Plus_Strand",
                                   "Expected_Mutations_Minus_Strand", "Expected_Mutations_Both_Strands",
                                   "Expected_Mutations_Aligned_Strands")
  }

  # Compute a factor to adjust normalized values based on the ratio of total
  # background:raw nucleosome mutation counts.
  totalCountsAdjust = sum(backgroundCounts$Expected_Mutations_Both_Strands) /
                      sum(rawCounts$Both_Strands_Counts)

  # Create a table of normalized values from the given data
  normalizedData = data.table::data.table()
  normalizedData[,Dyad_Position := backgroundCounts$Dyad_Position]
  normalizedData[,Normalized_Minus_Strand := mapply(normalize,rawCounts$Minus_Strand_Counts,
                                                    backgroundCounts$Expected_Mutations_Minus_Strand,
                                                    MoreArgs = list(totalCountsAdjust = totalCountsAdjust))]
  normalizedData[,Normalized_Plus_Strand := mapply(normalize,rawCounts$Plus_Strand_Counts,
                                                   backgroundCounts$Expected_Mutations_Plus_Strand,
                                                   MoreArgs = list(totalCountsAdjust = totalCountsAdjust))]
  normalizedData[,Normalized_Both_Strands := mapply(normalize,rawCounts$Both_Strands_Counts,
                                                    backgroundCounts$Expected_Mutations_Both_Strands,
                                                    MoreArgs = list(totalCountsAdjust = totalCountsAdjust))]
  normalizedData[,Normalized_Aligned_Strands := mapply(normalize,rawCounts$Aligned_Strands_Counts,
                                                       backgroundCounts$Expected_Mutations_Aligned_Strands,
                                                       MoreArgs = list(totalCountsAdjust = totalCountsAdjust))]

  # Write the normalized data to a new file. (If desired)
  if (writeNormalizedData) {
    data.table::fwrite(normalizedData, sep = '\t', file = normalizedDataFilePath)
  }

  return(normalizedData)

}

# Normalizes data by dividing raw by expected, except in the case where expected is 0.
# (To avoid dividing by zero)
normalize = function(raw, expected, totalCountsAdjust) {
  if (expected == 0) return(0)
  else return(raw/expected*totalCountsAdjust)
}
