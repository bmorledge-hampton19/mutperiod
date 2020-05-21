# The class to hold the data resulting from the periodicity and asymmetry analyses.
#' @export
NucPeriodData = setClass("NucPeriodData", slots = list(periodicityResults = className("data.table","data.table"),
                                                       asymmetryResults = className("data.table","data.table"),
                                                       MSIInputs = "character", MSSInputs = "character",
                                                       wilcoxinResult = "list"))

#' @export
generateNucPeriodData = function(mutationCountsFilePaths, outputFilePath,
                                 MSIFilePaths = '', MSSFilePaths = '',
                                 enforceInputNamingConventions = FALSE) {

  # Generate a list of prefixes for the data files as identifiers.
  filePrefixes = sapply(strsplit(basename(mutationCountsFilePaths),"_nucleosome"), function(x) x[1])

  # Prep for MSI vs MSS analysis if the respective file path lists were given.
  compareMS = (MSIFilePaths != '' && MSSFilePaths != '')
  if (compareMS) {
    MSIFilePrefixes = sapply(strsplit(basename(MSIFilePaths),"_nucleosome"), function(x) x[1])
    MSSFilePrefixes = sapply(strsplit(basename(MSSFilePaths),"_nucleosome"), function(x) x[1])
  }

  peakPeriodicities = numeric(length(filePrefixes))
  periodicityPValues = numeric(length(filePrefixes))
  periodicitySNRs = numeric(length(filePrefixes))

  generalAsymmetryTValue = numeric(length(filePrefixes))
  generalAsymmetryPValue = numeric(length(filePrefixes))

  peakAsymmetryTValue = numeric(length(filePrefixes))
  peakAsymmetryPValue = numeric(length(filePrefixes))

  valleyAsymmetryTValue = numeric(length(filePrefixes))
  valleyAsymmetryPValue = numeric(length(filePrefixes))

  for (i in 1:length(filePrefixes)) {

    # Read in the data, and normalize it if necessary.

    if (enforceInputNamingConventions) {
      if (endsWith(mutationCountsFilePaths[i],"counts.tsv")) {
        normalizedData = normalizeNucleosomeMutationCounts(mutationCountsFilePaths[i])
      } else if (endsWith(mutationCountsFilePaths[i],"normalized.tsv")) {
        normalizedData = data.table::fread(file = mutationCountsFilePaths[i])
      } else {
        stop(paste("The file,", mutationCountsFilePaths,
                   "does not follow naming conventions for raw or normalized nucleosome mutation counts."))
      }
    } else {
      normalizedData = data.table::fread(file = mutationCountsFilePaths[i])
    }


    ##### Periodicity Analysis #####

    # Calculate the periodicity of the data using a Lomb-Scargle periodiagram.
    lombResult = lomb::lsp(normalizedData[,.(Dyad_Position,Normalized_Both_Strands)],
                     type = "period", from = 2, to = 50, ofac = 100, plot = F)
    plot(normalizedData[,.(Dyad_Position,Normalized_Both_Strands)],type = 'b', main = filePrefixes[i])
    # Store the relevant results!
    peakPeriodicities[i] = lombResult$peak.at[1]
    periodicityPValues[i] = lombResult$p.value

    # Calculate the SNR
    noiseBooleanVector = (lombResult$scanned < lombResult$peak.at[1] - 0.5
                          | lombResult$scanned > lombResult$peak.at[1] + 0.5)
    periodicitySNRs[i] = lombResult$peak / median(lombResult$power[noiseBooleanVector])


    ##### Asymmetry Analysis #####

    # Run asymmetry analysis on the aligned strand counts.
    generalAsymmetryResult =
      t.test(normalizedData$Normalized_Aligned_Strands[6:73],
             rev(normalizedData$Normalized_Aligned_Strands[75:142]), paired = T)
    generalAsymmetryTValue[i] = generalAsymmetryResult$statistic
    generalAsymmetryPValue[i] = generalAsymmetryResult$p.value

    # Run asymmetry analysis on the peaks and valleys of the aligned strands.
    peakAsymmetryResult =
      runExtremeAnalysisSuite(normalizedData$Normalized_Aligned_Strands, dyadPosCutoff = 68)
    peakAsymmetryTValue[i] = peakAsymmetryResult$statistic
    peakAsymmetryPValue[i] = peakAsymmetryResult$p.value


    valleyAsymmetryResult =
      runExtremeAnalysisSuite(normalizedData$Normalized_Aligned_Strands, maxes = FALSE, dyadPosCutoff = 68)
    valleyAsymmetryTValue[i] = valleyAsymmetryResult$statistic
    valleyAsymmetryPValue[i] = valleyAsymmetryResult$p.value

    # # Negative controls for asymmetry analysis on peaks and valleys
    # runExtremeAnalysisSuite(normalizedData$Normalized_Both_Strands, dyadPosCutof = 68)
    # runExtremeAnalysisSuite(normalizedData$Normalized_Both_Strands, maxes = FALSE, dyadPosCutoff = 68)

  }

  # Create data.tables for all the results.
  periodicityResults = data.table::data.table(Data_Set=filePrefixes,Peak_Periodicity=peakPeriodicities,
                                              PValue=periodicityPValues,SNR=periodicitySNRs)

  asymmetryResults = data.table::data.table(Data_Set=filePrefixes,
                                            General_Asymmetry_TValue = generalAsymmetryTValue,
                                            General_Asymmetry_PValue = generalAsymmetryPValue,
                                            Peak_Asymmetry_TValue = peakAsymmetryTValue,
                                            Peak_Asymmetry_PValue = peakAsymmetryPValue,
                                            Valley_Asymmetry_TValue = valleyAsymmetryTValue,
                                            Valley_Asymmetry_PValue = valleyAsymmetryPValue)

  # Run the SNR wilcoxin's test if necessary.
  if (compareMS) {

    MSS_SNR = periodicityResults[Data_Set %in% MSSFilePrefixes, SNR]
    MSI_SNR = periodicityResults[Data_Set %in% MSIFilePrefixes, SNR]

    wilcoxinResult = wilcox.test(MSS_SNR, MSI_SNR)

  } else wilcoxinResult = list()

  # Create the data object to return
  nucPeriodData = NucPeriodData(periodicityResults = periodicityResults, asymmetryResults = asymmetryResults,
                                MSIInputs = MSIFilePaths, MSSInputs = MSSFilePaths, wilcoxinResult = wilcoxinResult)
  save(nucPeriodData, file = outputFilePath)
  return(nucPeriodData)

}
