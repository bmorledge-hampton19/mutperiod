# This function uses the deconstructSigs package to assign mutation signatures to given sets of mutation data.
#' @export
getMutSigs = function(inputFilePath) {

  # Read in the input data and format it for deconstructSigs
  inputData = fread(input = inputFilePath)
  print("Formatting input...")
  colnames(inputData) = c("Sample","chr","pos","ref","alt")
  sampleIDs = unique(inputData$Sample)
  trinucCounts = deconstructSigs::mut.to.sigs.input(inputData)

  # Deconstruct the signatures!
  print("Deconstructing Signatures...")

  # Function which returns the signatures for a given sample ID that are at least at the given threshold.
  # If multiple signatures meet the threshold, they are returned as a comma separated string.
  # If no signatures meet the threshold, "None" is returned instead.
  findDominantSigs = function(sampleID, threshold = .5) {

    baseMutSigs = deconstructSigs::whichSignatures(trinucCounts,
                                                   signatures.ref = deconstructSigs::signatures.nature2013,
                                                   contexts.needed = T, sample.id = sampleID)
    validMutSigs = names(baseMutSigs$weights)[baseMutSigs$weights >= threshold]
    validMutSigs = sapply(validMutSigs, function(x) strsplit(x, '.', fixed = TRUE)[[1]][2])
    validMutSigs = paste(validMutSigs, collapse = ',')

    if (validMutSigs == '') validMutSigs = "None"

    return(validMutSigs)

  }

  # Create and return a data table of the sample IDs and their dominant signatures.
  return(data.table(Cohort_ID = sampleIDs, Signatures = sapply(sampleIDs, findDominantSigs)))

}
