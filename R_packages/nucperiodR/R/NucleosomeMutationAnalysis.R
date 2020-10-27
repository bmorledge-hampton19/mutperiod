#' @export
generateNucPeriodData = function(mutationCountsFilePaths, outputFilePath,
                                 filePathGroup1 = '', filePathGroup2 = '',
                                 nucleosomeDyadPosCutoff = 60,
                                 nucleosomeMutationCutoff = 5000,
                                 enforceInputNamingConventions = FALSE,
                                 outputGraphs = FALSE) {

  # Get the raw mutation counts for each given mutation counts file.
  # If naming conventions aren't enforced, it is assumed that the given mutation counts file will suffice.
  # Otherwise, if conventions are enforced and a normalized file is given,
  # the path to the file with raw counts is generated.
  if (enforceInputNamingConventions) {
    rawCountsFilePaths = sapply(mutationCountsFilePaths, getRawCountsFilePath)
  } else rawCountsFilePaths = mutationCountsFilePaths

  dataSetNames = getDataSetNames(mutationCountsFilePaths, enforceInputNamingConventions)
  rawNucleosomeMutationCounts = sapply(rawCountsFilePaths, getRawNucleosomeMutationCounts,
                                       nucleosomeDyadPosCutoff)

  rawCountsTable = data.table::data.table(Associated_File_Path = mutationCountsFilePaths,
                                          Associated_Data_Set = dataSetNames,
                                          Raw_Nucleosome_Mutation_Counts = rawNucleosomeMutationCounts)
  # Sort the newly created table by data set names.
  data.table::setorder(rawCountsTable,Associated_Data_Set)

  # If a cutoff greater than 0 was given for nucleosome mutation counts, enforce it.
  if (nucleosomeMutationCutoff > 0) {

    filteredCountsTable = rawCountsTable[mapply(filterCounts, Raw_Nucleosome_Mutation_Counts,
                                                nucleosomeMutationCutoff, Associated_Data_Set)]
    validFilePaths = filteredCountsTable$Associated_File_Path
    validDataSetNames = filteredCountsTable$Associated_Data_Set

  } else {
    validFilePaths = rawCountsTable$Associated_File_Path
    validDataSetNames = rawCountsTable$Associated_Data_Set
  }

  # Generate a list of raw counts tables from the set of valid file paths.
  validRawFilePaths = unique(sapply(validFilePaths, getRawCountsFilePath))
  rawNucleosomeCountsTables = lapply(validRawFilePaths, function(x) data.table::fread(file = x))
  names(rawNucleosomeCountsTables) = getDataSetNames(validRawFilePaths, enforceInputNamingConventions)

  # Prep for group comparison if the respective file path lists were given.
  compareGroups = (filePathGroup1 != '' && filePathGroup2 != '')
  if (compareGroups) {
    group1DataSetNames = sapply(strsplit(basename(filePathGroup1),"_nucleosome"), function(x) x[1])
    group2DataSetNames = sapply(strsplit(basename(filePathGroup2),"_nucleosome"), function(x) x[1])
  }

  normalizedNucleosomeCountsTables = vector("list",length(validFilePaths))
  names(normalizedNucleosomeCountsTables) = validDataSetNames

  peakPeriodicities = numeric(length(validFilePaths))
  periodicityPValues = numeric(length(validFilePaths))
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
    normalizedNucleosomeCountsTables[[i]] = nucleosomeCountsData

    # Adjust the dyadPosCutoff if we have translational periodicity data.
    if (grepl("nuc-group", basename(validFilePaths[i]), fixed = TRUE)) {
      dyadPosCutoff = 1000
      lombFrom = 50
      lombTo = 250
    } else {
      dyadPosCutoff = nucleosomeDyadPosCutoff
      lombFrom = 7
      lombTo = 20
    }

    ##### Periodicity Analysis #####

    print("Running periodicity analysis...")

    # Calculate the periodicity of the data using a Lomb-Scargle periodiagram.
    lombResult = lomb::lsp(nucleosomeCountsData[Dyad_Position >= -dyadPosCutoff & Dyad_Position <= dyadPosCutoff,
                                                  .(Dyad_Position,Normalized_Both_Strands)],
                     type = "period", from = lombFrom, to = lombTo, ofac = 100, plot = outputGraphs)
    if (outputGraphs) {
      plot(nucleosomeCountsData[Dyad_Position >= -dyadPosCutoff & Dyad_Position <= dyadPosCutoff,
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

  }

  # Create data.tables for all the results.
  periodicityResults = data.table::data.table(Data_Set=validDataSetNames,Peak_Periodicity=peakPeriodicities,
                                              PValue=periodicityPValues,SNR=periodicitySNRs)

  # Run the SNR wilcoxin's test if necessary.
  if (compareGroups) {

    print("Comparing periodicity results between given groups.")

    group1SNR = periodicityResults[Data_Set %in% group1DataSetNames, SNR]
    group2SNR = periodicityResults[Data_Set %in% group2DataSetNames, SNR]

    if ( length(group1SNR) + length(group2SNR) > nrow(periodicityResults)) {
      warning("The number of the group1 and group2 combined SNR values is greater than the total number of SNR values")
    } else if ( length(group1SNR) + length(group2SNR) < nrow(periodicityResults)) {
      warning("The number of the group1 and group2 combined SNR values is less than the total number of SNR values")
    }

    wilcoxinResult = wilcox.test(group1SNR, group2SNR)
    print(paste("wilcoxin test p-value:", wilcoxinResult$p.value))

  }

  # Create the data object to return
  nucPeriodData = list(normalizedNucleosomeCountsTables = normalizedNucleosomeCountsTables,
                       rawNucleosomeCountsTables = rawNucleosomeCountsTables,
                       periodicityResults = periodicityResults)

  # Add the group comparison results if requested.
  if (compareGroups) {
    nucPeriodData = append(nucPeriodData, list(group1Inputs = group1DataSetNames,
                                               group2Inputs = group2DataSetNames,
                                               wilcoxinResult = wilcoxinResult))
  }

  save(nucPeriodData, file = outputFilePath)

}

getRawCountsFilePath = function(rawOrNormalizedCountsFilePath) {

  fileName = basename(rawOrNormalizedCountsFilePath)

  if (grepl("raw_nucleosome",fileName)) {
    return(rawOrNormalizedCountsFilePath)
  } else {
    dataSetName = strsplit(fileName,"singlenuc|trinuc|pentanuc")[[1]][1]
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
