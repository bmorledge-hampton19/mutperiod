# Given a path to a file of raw mutation counts with a corresponding file of background counts,
# Read in the data and produce normalized counts, writing them to a new file if specified (by default true).
#' @export
normalizeNucleosomeMutationCounts = function(rawCountsFilePath, backgroundCountsFilePath,
                                             normalizedDataFilePath = NULL, alternativeScalingFactor = NULL) {

  # Read in the data
  rawCounts = data.table::fread(file = rawCountsFilePath)
  backgroundCounts = data.table::fread(file = backgroundCountsFilePath)

  # If the background file is "custom input" it will have the same column names as a raw counts file.
  # Check for the raw counts headers, and replace them with the expected headers if necessary.
  if (colnames(backgroundCounts)[2] == "Plus_Strand_Counts") {
    colnames(backgroundCounts) = c("Dyad_Position", "Expected_Mutations_Plus_Strand",
                                   "Expected_Mutations_Minus_Strand", "Expected_Mutations_Both_Strands",
                                   "Expected_Mutations_Aligned_Strands")
  }

  # Check to see if one of the tables has half-base positions while the other does not.
  # If this turns out to be the case, lower the resolution of the half-base positions by
  # splitting them up into the two adjacent integer positions.  Also warn the user.
  if (as.integer(backgroundCounts[[1]][1]) != as.numeric(backgroundCounts[[1]][1]) &&
      as.integer(rawCounts[[1]][1]) == as.numeric(rawCounts[[1]][1])) {

    warning("Background counts use half-base positions but raw counts do not.  ",
            "Splitting background counts into adjacent integer positions and omitting a base on each end.")
    rawCounts = rawCounts[c(-1,-dim(rawCounts)[1])]

    backgroundCounts = data.table::rbindlist(lapply(1:(dim(backgroundCounts)[1] - 1),
                                                    function(x) (backgroundCounts[x] + backgroundCounts[x+1]) / 2))

  } else if (as.integer(backgroundCounts[[1]][1]) == as.numeric(backgroundCounts[[1]][1]) &&
             as.integer(rawCounts[[1]][1]) != as.numeric(rawCounts[[1]][1])) {

    warning("Raw counts use half-base positions but background counts do not.  ",
            "Splitting raw counts into adjacent integer positions and omitting a base on each end.")
    backgroundCounts = backgroundCounts[c(-1,-dim(backgroundCounts)[1])]

    rawCounts = data.table::rbindlist(lapply(1:(dim(rawCounts)[1] - 1),
                                             function(x) (rawCounts[x] + rawCounts[x+1]) / 2))

  }

  # Make sure that the raw and background counts both use the same number of dyad positions.
  if (dim(rawCounts)[1] != dim(backgroundCounts)[1]) {
    stop("Unequal dyad positions in raw vs. background counts data.")
  }

  # Compute a basic scaling factor from the ratio of total
  # background:raw nucleosome mutation counts.  (Centers final normalized counts on 1)
  scalingFactor = sum(backgroundCounts$Expected_Mutations_Both_Strands) /
                      sum(rawCounts$Both_Strands_Counts)

  # Create a table of normalized values from the given data
  normalizedData = data.table::data.table()
  normalizedData[,Dyad_Position := backgroundCounts$Dyad_Position]
  normalizedData[,Normalized_Minus_Strand := mapply(normalize,rawCounts$Minus_Strand_Counts,
                                                    backgroundCounts$Expected_Mutations_Minus_Strand,
                                                    MoreArgs = list(scalingFactor = scalingFactor))]
  normalizedData[,Normalized_Plus_Strand := mapply(normalize,rawCounts$Plus_Strand_Counts,
                                                   backgroundCounts$Expected_Mutations_Plus_Strand,
                                                   MoreArgs = list(scalingFactor = scalingFactor))]
  normalizedData[,Normalized_Both_Strands := mapply(normalize,rawCounts$Both_Strands_Counts,
                                                    backgroundCounts$Expected_Mutations_Both_Strands,
                                                    MoreArgs = list(scalingFactor = scalingFactor))]
  normalizedData[,Normalized_Aligned_Strands := mapply(normalize,rawCounts$Aligned_Strands_Counts,
                                                       backgroundCounts$Expected_Mutations_Aligned_Strands,
                                                       MoreArgs = list(scalingFactor = scalingFactor))]

  # If an alternative scaling factor was given, compute additional normalized values with this factor
  if (!is.null(alternativeScalingFactor)) {
    normalizedData[,Alternative_Normalized_Minus_Strand :=
                     Normalized_Minus_Strand/scalingFactor*alternativeScalingFactor]
    normalizedData[,Alternative_Normalized_Plus_Strand :=
                     Normalized_Plus_Strand/scalingFactor*alternativeScalingFactor]
    normalizedData[,Alternative_Normalized_Both_Strands :=
                     Normalized_Both_Strands/scalingFactor*alternativeScalingFactor]
    normalizedData[,Alternative_Normalized_Aligned_Strands :=
                     Normalized_Aligned_Strands/scalingFactor*alternativeScalingFactor]
  }

  # Write the normalized data to a new file. (If desired)
  if (!is.null(normalizedDataFilePath)) {
    data.table::fwrite(normalizedData, sep = '\t', file = normalizedDataFilePath)
  }

  return(normalizedData)

}

# Normalizes data by dividing raw by expected, except in the case where expected is 0.
# (To avoid dividing by zero)
normalize = function(raw, expected, scalingFactor) {
  if (expected == 0) return(0)
  else return(raw/expected*scalingFactor)
}
