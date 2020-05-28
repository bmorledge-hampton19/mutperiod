#' @export
generateNucPeriodData = function(mutationCountsFilePaths, outputFilePath,
                                 MSIFilePaths = '', MSSFilePaths = '',
                                 enforceInputNamingConventions = FALSE,
                                 outputGraphs = FALSE) {

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

    print(paste("Working with", filePrefixes[i]))

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

    print("Running periodicity analysis...")

    # Calculate the periodicity of the data using a Lomb-Scargle periodiagram.
    lombResult = lomb::lsp(normalizedData[,.(Dyad_Position,Normalized_Both_Strands)],
                     type = "period", from = 2, to = 50, ofac = 100, plot = outputGraphs)
    if (outputGraphs) {
      plot(normalizedData[,.(Dyad_Position,Normalized_Both_Strands)],type = 'b', main = filePrefixes[i])
    }
    # Store the relevant results!
    peakPeriodicities[i] = lombResult$peak.at[1]
    periodicityPValues[i] = lombResult$p.value

    # Calculate the SNR
    noiseBooleanVector = (lombResult$scanned < lombResult$peak.at[1] - 0.5
                          | lombResult$scanned > lombResult$peak.at[1] + 0.5)
    periodicitySNRs[i] = lombResult$peak / median(lombResult$power[noiseBooleanVector])


    ##### Asymmetry Analysis #####

    print("Running asymmetry analysis")

    # Run asymmetry analysis on the aligned strand counts.
    generalAsymmetryResult =
      t.test(normalizedData$Normalized_Aligned_Strands[6:73],
             rev(normalizedData$Normalized_Aligned_Strands[75:142]), paired = T)
    generalAsymmetryTValue[i] = generalAsymmetryResult$statistic
    generalAsymmetryPValue[i] = generalAsymmetryResult$p.value

    # Run asymmetry analysis on the peaks and valleys of the aligned strands.
    # If there is not enough data to derive peaks and valleys (e.g. for individual donors)
    # the result will be NA.
    peakAsymmetryResult =
      runExtremeAnalysisSuite(normalizedData$Normalized_Aligned_Strands, dyadPosCutoff = 68)
    if (any(!is.na(peakAsymmetryResult))) {
      peakAsymmetryTValue[i] = peakAsymmetryResult$statistic
      peakAsymmetryPValue[i] = peakAsymmetryResult$p.value
    } else {
      peakAsymmetryTValue[i] = NA
      peakAsymmetryPValue[i] = NA
    }

    valleyAsymmetryResult =
      runExtremeAnalysisSuite(normalizedData$Normalized_Aligned_Strands, maxes = FALSE, dyadPosCutoff = 68)
    if (any(!is.na(valleyAsymmetryResult))) {
      valleyAsymmetryTValue[i] = valleyAsymmetryResult$statistic
      valleyAsymmetryPValue[i] = valleyAsymmetryResult$p.value
    } else {
      valleyAsymmetryTValue[i] = NA
      valleyAsymmetryPValue[i] = NA
    }

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

    print("Comparing periodicity results based on microsatellite stability.")

    MSS_SNR = periodicityResults[Data_Set %in% MSSFilePrefixes, SNR]
    MSI_SNR = periodicityResults[Data_Set %in% MSIFilePrefixes, SNR]

    if ( length(MSS_SNR) + length(MSI_SNR) > nrow(periodicityResults)) {
      warning("The sum of the MSS and MSI Data is greater than the total number of SNR values")
    } else if ( length(MSS_SNR) + length(MSI_SNR) < nrow(periodicityResults)) {
      warning("The sum of the MSS and MSI Data is less than the total number of SNR values")
    }

    wilcoxinResult = wilcox.test(MSS_SNR, MSI_SNR)
    print(paste("wilcoxin test p-value:", wilcoxinResult$p.value))

  }

  # Create the data object to return
  if (compareMS) {
    nucPeriodData = list(periodicityResults = periodicityResults, asymmetryResults = asymmetryResults,
                         MSIInputs = MSIFilePrefixes, MSSInputs = MSSFilePrefixes, wilcoxinResult = wilcoxinResult)
  } else {
    nucPeriodData = list(periodicityResults = periodicityResults, asymmetryResults = asymmetryResults)
  }
  save(nucPeriodData, file = outputFilePath)

}
