#' @export
generateMutperiodData = function(mutationCountsFilePaths, expectedPeakPeriodicities = NA,
                                 filePathGroup1 = '', filePathGroup2 = '',
                                 nucleosomeDyadPosCutoff = 60,
                                 nucleosomeMutationCutoff = 5000,
                                 enforceInputNamingConventions = FALSE,
                                 alignStrands = FALSE,
                                 outputGraphs = FALSE,
                                 overridePeakPeriodicityWithExpectedPeak = FALSE) {

  # Get the raw mutation counts for each given mutation counts file.
  # If naming conventions aren't enforced, it is assumed that the given mutation counts file will suffice.
  # Otherwise, if conventions are enforced and a normalized file is given,
  # the path to the file with raw counts is generated.
  if (enforceInputNamingConventions) {
    rawCountsFilePaths = sapply(mutationCountsFilePaths, getRawCountsFilePath)
  } else rawCountsFilePaths = mutationCountsFilePaths
  rawInputFilePaths = rawCountsFilePaths[rawCountsFilePaths == mutationCountsFilePaths]

  dataSetNames = getDataSetNames(mutationCountsFilePaths, enforceInputNamingConventions)
  rawNucleosomeMutationCounts = sapply(rawCountsFilePaths, getRawNucleosomeMutationCounts,
                                       nucleosomeDyadPosCutoff)

  rawCountsTable = data.table::data.table(Associated_File_Path = mutationCountsFilePaths,
                                          Associated_Data_Set = dataSetNames,
                                          Raw_Nucleosome_Mutation_Counts = rawNucleosomeMutationCounts,
                                          Expected_Peak_Periodicity = expectedPeakPeriodicities)
  # Sort the newly created table by data set names.
  data.table::setorder(rawCountsTable,Associated_Data_Set)

  # If a cutoff greater than 0 was given for nucleosome mutation counts, enforce it.
  if (nucleosomeMutationCutoff > 0) {

    filteredCountsTable = rawCountsTable[mapply(filterCounts, Raw_Nucleosome_Mutation_Counts,
                                                nucleosomeMutationCutoff, Associated_Data_Set)]
    validFilePaths = filteredCountsTable$Associated_File_Path
    expectedPeakPeriodicities = filteredCountsTable$Expected_Peak_Periodicity
    validDataSetNames = filteredCountsTable$Associated_Data_Set

  } else {
    validFilePaths = rawCountsTable$Associated_File_Path
    expectedPeakPeriodicities = rawCountsTable$Expected_Peak_Periodicity
    validDataSetNames = rawCountsTable$Associated_Data_Set
  }

  # Generate a list of raw counts tables from the set of valid file paths.
  validRawFilePaths = unique(sapply(validFilePaths, getRawCountsFilePath))
  rawNucleosomeCountsTables = lapply(validRawFilePaths, function(x) data.table::fread(file = x))
  names(rawNucleosomeCountsTables) = getDataSetNames(validRawFilePaths, enforceInputNamingConventions)

  # Prep for group comparison if the respective file path lists were given.
  compareGroups = (filePathGroup1 != '' && filePathGroup2 != '')
  if (compareGroups) {
    group1DataSetNames = sort(sapply(strsplit(basename(filePathGroup1),"_nucleosome"), function(x) x[1]))
    group2DataSetNames = sort(sapply(strsplit(basename(filePathGroup2),"_nucleosome"), function(x) x[1]))
  }

  # Determine which of the valid file paths are normalized and initialize the associated lists.
  normalizedValidFilePaths = validFilePaths[!(validFilePaths %in% rawInputFilePaths)]
  normalizedNucleosomeCountsTables = lapply(normalizedValidFilePaths,
                                            function(x) data.table::fread(file = x))
  names(normalizedNucleosomeCountsTables) = validDataSetNames[!(validFilePaths %in% rawInputFilePaths)]

  peakPeriodicities = numeric(length(validFilePaths))
  relevantPeriodicityPowers = numeric(length(validFilePaths))
  periodicitySNRs = numeric(length(validFilePaths))

  for (i in 1:length(validFilePaths)) {

    # Make sure that if no files validated, no analysis attempts to run.
    if (length(validFilePaths) == 0) {
      warning("No valid files given.  Returning blank analysis.")
      break()
    }

    print(paste("Working with", validDataSetNames[i]))

    # Read in the data.
    nucleosomeCountsData = data.table::fread(file = validFilePaths[i])

    # Adjust the dyadPosCutoff if we have translational periodicity data.
    if (grepl("nuc-group", basename(validFilePaths[i]), fixed = TRUE)) {
      dyadPosCutoff = 1000
      lombFrom = 50
      lombTo = 250
    } else {
      dyadPosCutoff = nucleosomeDyadPosCutoff
      lombFrom = 5
      lombTo = 25
    }

    # Determine whether or not the current data is normalized or not, and
    # derive the necessary counts vector from it.
    if (validFilePaths[[i]] %in% rawInputFilePaths) {
      if (alignStrands) {
        counts = nucleosomeCountsData[Dyad_Position >= -dyadPosCutoff & Dyad_Position <= dyadPosCutoff,
                                      .(Dyad_Position,Aligned_Strands_Counts)]
      } else {
        counts = nucleosomeCountsData[Dyad_Position >= -dyadPosCutoff & Dyad_Position <= dyadPosCutoff,
                                      .(Dyad_Position,Both_Strands_Counts)]
      }
    } else {
      if (alignStrands) {
        counts = nucleosomeCountsData[Dyad_Position >= -dyadPosCutoff & Dyad_Position <= dyadPosCutoff,
                                      .(Dyad_Position,Normalized_Aligned_Strands)]
      } else {
        counts = nucleosomeCountsData[Dyad_Position >= -dyadPosCutoff & Dyad_Position <= dyadPosCutoff,
                                      .(Dyad_Position,Normalized_Both_Strands)]
      }
    }

    ##### Periodicity Analysis #####

    print("Running periodicity analysis...")

    # Calculate the periodicity of the data using a Lomb-Scargle periodiagram.
    lombResult = lomb::lsp(counts, type = "period", from = lombFrom, to = lombTo,
                           ofac = 100, plot = FALSE)

    # Get the power associated with the periodicity being used (either peak or expected)
    if (overridePeakPeriodicityWithExpectedPeak) {
      relevantPeriodicity = expectedPeakPeriodicities[i]
      closestPeriodicity = lombResult$scanned[which.min(abs(lombResult$scanned - relevantPeriodicity))]
      if (abs(closestPeriodicity - relevantPeriodicity) > 0.1) {
        stop(paste0("No scanned periodicities exist within 0.1 units of the given expected periodicity, ",
                    relevantPeriodicity,".  Closest scanned periodicity is at ",closestPeriodicity,"."))
      }
      relevantPower = lombResult$power[lombResult$scanned == closestPeriodicity]
    } else {
      relevantPeriodicity = lombResult$peak.at[1]
      relevantPower = lombResult$peak
    }



    # Store the relevant periodicity as well as its power and p-value
    peakPeriodicities[i] = lombResult$peak.at[1]
    relevantPeriodicityPowers[i] = relevantPower

    # Calculate the SNR, then store it.
    noiseBooleanVector = (lombResult$scanned < relevantPeriodicity - 0.5
                          | lombResult$scanned > relevantPeriodicity + 0.5)
    periodicitySNRs[i] = relevantPower / median(lombResult$power[noiseBooleanVector])

  }

  # Create data.tables for all the results.
  periodicityResults = data.table::data.table(Data_Set = validDataSetNames,
                                              Peak_Periodicity = peakPeriodicities,
                                              Expected_Peak_Periodicity = expectedPeakPeriodicities,
                                              Power = relevantPeriodicityPowers,
                                              SNR = periodicitySNRs)
  data.table::setorder(periodicityResults, Data_Set)

  # Run the SNR wilcoxon's test if necessary.
  if (compareGroups) {

    print("Comparing periodicity results between given groups.")

    group1SNR = periodicityResults[Data_Set %in% group1DataSetNames, SNR]
    group2SNR = periodicityResults[Data_Set %in% group2DataSetNames, SNR]

    if ( length(group1SNR) + length(group2SNR) > nrow(periodicityResults)) {
      warning("The number of the group1 and group2 combined SNR values is greater than the total number of SNR values")
    } else if ( length(group1SNR) + length(group2SNR) < nrow(periodicityResults)) {
      warning("The number of the group1 and group2 combined SNR values is less than the total number of SNR values")
    }

    wilcoxonResult = wilcox.test(group1SNR, group2SNR)
    print(paste("wilcoxon test p-value:", wilcoxonResult$p.value))

  }

  # Create a string to represent the type of relevant periodicity
  if (overridePeakPeriodicityWithExpectedPeak) {
    periodicityUsedForPowerAnalysis = "Expected Peak Periodicity"
  } else {
    periodicityUsedForPowerAnalysis = "Peak Periodicity"
  }

  # Create a string to represent the strand data used.
  if (alignStrands) {
    strandPositioning = "Aligned"
  } else {
    strandPositioning = "Antiparallel"
  }

  # Create the data object to return
  mutperiodData = list(normalizedNucleosomeCountsTables = normalizedNucleosomeCountsTables,
                       rawNucleosomeCountsTables = rawNucleosomeCountsTables,
                       periodicityUsedForPowerAnalysis = periodicityUsedForPowerAnalysis,
                       strandPositioning = strandPositioning,
                       periodicityResults = periodicityResults)

  # Add the group comparison results if requested.
  if (compareGroups) {
    mutperiodData = append(mutperiodData, list(group1Inputs = group1DataSetNames,
                                               group2Inputs = group2DataSetNames,
                                               wilcoxonResult = wilcoxonResult))
  }

  return(mutperiodData)

}

getRawCountsFilePath = function(rawOrNormalizedCountsFilePath) {

  fileName = basename(rawOrNormalizedCountsFilePath)

  if (grepl("raw_nucleosome",fileName)) {
    return(rawOrNormalizedCountsFilePath)
  } else {
    dataSetName = strsplit(fileName,"singlenuc|dinuc|trinuc|quadrunuc|pentanuc|hexanuc|custom_context")[[1]][1]
    if (grepl("nuc-group",fileName,fixed = TRUE)) {
      dyadRadius = "nuc-group_"
    }  else dyadRadius = ''
    if (grepl("linker+",fileName,fixed = TRUE)) {
      linkerOffset = grep("linker+",strsplit(fileName,'_')[[1]],value = TRUE)
      linkerOffset = paste0(linkerOffset,'_')
    } else {
      linkerOffset = ''
    }
    return(file.path(dirname(rawOrNormalizedCountsFilePath),
                     paste0(dataSetName,linkerOffset,dyadRadius,"raw_nucleosome_mutation_counts.tsv")))
  }

}

getDataSetNames = function(mutationCountsFilePaths, enforceInputNamingConventions) {

  if (enforceInputNamingConventions) {
    return(sapply(strsplit(basename(mutationCountsFilePaths),"_nucleosome_mutation_counts"), function(x) x[1]))
  } else {
    return(sapply(strsplit(basename(mutationCountsFilePaths),'.', fixed = TRUE), function(x) x[1]))
  }

}

getRawNucleosomeMutationCounts = function(countsFilePath, dyadPosCutoff) {

  countsTable = data.table::fread(file = countsFilePath)
  return(sum(countsTable[Dyad_Position >= -dyadPosCutoff & Dyad_Position <= dyadPosCutoff,
                         Both_Strands_Counts]))

}

filterCounts = function(counts, cutoff, dataGroup) {

  if (counts >= cutoff) return(TRUE)
  else {
    print(paste(dataGroup,"is not valid with only",counts,"total counts and will be filtered out."))
    return(FALSE)
  }

}

# A helper function for when I want to quickly get lomb results with SNR.
# If the lomb results were computed elsewhere, they can be provided with the precomputedLombResult argument
# (In the above case, all other required values are irrelevant and should be set to NA)
#' @export
getPeakPeriodicityAndSNR = function(counts, lombFrom, lombTo, plot = FALSE, precomputedLombResult = NA) {

  # Calculate the periodicity of the data using a Lomb-Scargle periodiagram. (Or retrieve the precomputed result)
  if (!is.na(precomputedLombResult)) {
    lombResult = precomputedLombResult
  } else {
    lombResult = lomb::lsp(counts, type = "period", from = lombFrom, to = lombTo,
                           ofac = 100, plot = plot)
  }

  # Get the peak periodicity and its associated SNR
  peakPeriodicity = lombResult$peak.at[1]
  noiseBooleanVector = (lombResult$scanned < peakPeriodicity - 0.5
                        | lombResult$scanned > peakPeriodicity + 0.5)
  SNR = lombResult$peak / median(lombResult$power[noiseBooleanVector])

  return(c(peakPeriodicity,SNR))

}
