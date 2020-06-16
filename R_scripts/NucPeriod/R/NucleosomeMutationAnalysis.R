#' @export
generateNucPeriodData = function(mutationCountsFilePaths, outputFilePath,
                                 MSIFilePaths = '', MSSFilePaths = '',
                                 dyadPosCutoff = 60,
                                 nucleosomeMutationCutoff = 5000,
                                 enforceInputNamingConventions = FALSE,
                                 outputGraphs = FALSE) {

  # Get the raw mutation counts for each given mutation counts file.
  # If naming conventions aren't enforced, it is assumed that the given mutation counts file will suffice.
  if (enforceInputNamingConventions) {
    rawCountsFilePaths = sapply(strsplit(mutationCountsFilePaths,"_nucleosome_mutation_counts"),
                                function(x) paste0(x[1],"_nucleosome_mutation_counts.tsv"))
  } else rawCountsFilePaths = mutationCountsFilePaths

  dataSetNames = getDataSetNames(rawCountsFilePaths, enforceInputNamingConventions)
  rawNucleosomeMutationCounts = sapply(rawCountsFilePaths, getRawNucleosomeMutationCounts)

  mutationCountsTable = data.table::data.table(File_Path = mutationCountsFilePaths,
                                               Data_Set = dataSetNames,
                                               Raw_Nucleosome_Mutation_Counts = rawNucleosomeMutationCounts)

  # If a cutoff greater than 0 was given for nucleosome mutation counts, enforce it.
  if (nucleosomeMutationCutoff > 0) {

    filteredCountsTable = mutationCountsTable[Raw_Nucleosome_Mutation_Counts >= nucleosomeMutationCutoff]
    validFilePaths = filteredCountsTable$File_Path
    validDataSetNames = filteredCountsTable$Data_Set

  } else {
    validFilePaths = mutationCountsPaths
    validDataSetNames = dataSetNames
  }

  # Prep for MSI vs MSS analysis if the respective file path lists were given.
  compareMS = (MSIFilePaths != '' && MSSFilePaths != '')
  if (compareMS) {
    MSIDataSetNames = sapply(strsplit(basename(MSIFilePaths),"_nucleosome"), function(x) x[1])
    MSSDataSetNames = sapply(strsplit(basename(MSSFilePaths),"_nucleosome"), function(x) x[1])
  }

  peakPeriodicities = numeric(length(validFilePaths))
  periodicityPValues = numeric(length(validFilePaths))
  periodicitySNRs = numeric(length(validFilePaths))

  generalAsymmetryTValue = numeric(length(validFilePaths))
  generalAsymmetryPValue = numeric(length(validFilePaths))

  peakAsymmetryTValue = numeric(length(validFilePaths))
  peakAsymmetryPValue = numeric(length(validFilePaths))

  valleyAsymmetryTValue = numeric(length(validFilePaths))
  valleyAsymmetryPValue = numeric(length(validFilePaths))

  for (i in 1:length(validFilePaths)) {

    print(paste("Working with", validDataSetNames[i]))

    # Read in the data, and normalize it if necessary.

    if (enforceInputNamingConventions) {
      if (endsWith(validFilePaths[i],"counts.tsv")) {
        normalizedData = normalizeNucleosomeMutationCounts(validFilePaths[i])
      } else if (endsWith(validFilePaths[i],"normalized.tsv")) {
        normalizedData = data.table::fread(file = validFilePaths[i])
      } else {
        stop(paste("The file,", validFilePaths[i],
                   "does not follow naming conventions for raw or normalized nucleosome mutation counts."))
      }
    } else {
      normalizedData = data.table::fread(file = validFilePaths[i])
    }


    ##### Periodicity Analysis #####

    print("Running periodicity analysis...")

    # Calculate the periodicity of the data using a Lomb-Scargle periodiagram.
    lombResult = lomb::lsp(normalizedData[Dyad_Position >= -dyadPosCutoff & Dyad_Position <= dyadPosCutoff,
                                          .(Dyad_Position,Normalized_Both_Strands)],
                     type = "period", from = 2, to = 50, ofac = 100, plot = outputGraphs)
    if (outputGraphs) {
      plot(normalizedData[Dyad_Position >= -dyadPosCutoff & Dyad_Position <= dyadPosCutoff,
                          .(Dyad_Position,Normalized_Both_Strands)],
           type = 'b', main = validDataSetNames[i])
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
      t.test(normalizedData$Normalized_Aligned_Strands[(74-dyadPosCutoff):73],
             rev(normalizedData$Normalized_Aligned_Strands[75:(74+dyadPosCutoff)]), paired = T)
    generalAsymmetryTValue[i] = generalAsymmetryResult$statistic
    generalAsymmetryPValue[i] = generalAsymmetryResult$p.value

    # Run asymmetry analysis on the peaks and valleys of the aligned strands.
    # If there is not enough data to derive peaks and valleys (e.g. for individual donors)
    # the result will be NA.
    peakAsymmetryResult = runExtremeAnalysisSuite(normalizedData$Normalized_Aligned_Strands,
                                                  dyadPosCutoff = dyadPosCutoff)
    if (any(!is.na(peakAsymmetryResult))) {
      peakAsymmetryTValue[i] = peakAsymmetryResult$statistic
      peakAsymmetryPValue[i] = peakAsymmetryResult$p.value
    } else {
      peakAsymmetryTValue[i] = NA
      peakAsymmetryPValue[i] = NA
    }

    valleyAsymmetryResult = runExtremeAnalysisSuite(normalizedData$Normalized_Aligned_Strands,
                                                    maxes = FALSE, dyadPosCutoff = dyadPosCutoff)
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
  periodicityResults = data.table::data.table(Data_Set=validDataSetNames,Peak_Periodicity=peakPeriodicities,
                                              PValue=periodicityPValues,SNR=periodicitySNRs)

  asymmetryResults = data.table::data.table(Data_Set=validDataSetNames,
                                            General_Asymmetry_TValue = generalAsymmetryTValue,
                                            General_Asymmetry_PValue = generalAsymmetryPValue,
                                            Peak_Asymmetry_TValue = peakAsymmetryTValue,
                                            Peak_Asymmetry_PValue = peakAsymmetryPValue,
                                            Valley_Asymmetry_TValue = valleyAsymmetryTValue,
                                            Valley_Asymmetry_PValue = valleyAsymmetryPValue)

  # Run the SNR wilcoxin's test if necessary.
  if (compareMS) {

    print("Comparing periodicity results based on microsatellite stability.")

    MSS_SNR = periodicityResults[Data_Set %in% MSSDataSetNames, SNR]
    MSI_SNR = periodicityResults[Data_Set %in% MSIDataSetNames, SNR]

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
    nucPeriodData = list(countsData = mutationCountsTable,
                         periodicityResults = periodicityResults, asymmetryResults = asymmetryResults,
                         MSIInputs = MSIDataSetNames, MSSInputs = MSSDataSetNames, wilcoxinResult = wilcoxinResult)
  } else {
    nucPeriodData = list(countsData = mutationCountsTable,
                         periodicityResults = periodicityResults, asymmetryResults = asymmetryResults)
  }
  save(nucPeriodData, file = outputFilePath)

}

getDataSetNames = function(mutationCountsFilePaths, enforceInputNamingConventions) {

  if (enforceInputNamingConventions) {
    return(sapply(strsplit(basename(mutationCountsFilePaths),"_nucleosome"), function(x) x[1]))
  } else {
    return(sapply(strsplit(basename(mutationCountsFilePaths),'.', fixed = TRUE), function(x) x[1]))
  }

}

getRawNucleosomeMutationCounts = function(countsFilePath) {

  countsTable = data.table::fread(file = countsFilePath)
  return(sum(countsTable$Both_Strands_Counts))

}
