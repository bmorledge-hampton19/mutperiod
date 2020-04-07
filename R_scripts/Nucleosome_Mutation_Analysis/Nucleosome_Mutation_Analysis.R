#Ensure that BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("lomb")
BiocManager::install("quantmod")
BiocManager::install("bspec")

library(lomb)
library(data.table)
library(bspec)

source("Parse_Nucleosome_Mutation_Data.R")
source("Assymetry_Analysis.R")

filePaths = choose.files(caption = "Select Raw Nucleosome Mutation Data.",
                         filters = c(c("Tab Separated Values (*.tsv)","Any files"),c("*.tsv","*.*")), index = 1)

# Generate a list of prefixes for the data files so that iterating through them is easier later.
filePrefixes = sapply(strsplit(basename(filePaths),"_nucleosome"), function(x) x[1])

peakPeriodicities = numeric(length(filePrefixes))
periodicityPValues = numeric(length(filePrefixes))
periodicitySNRs = numeric(length(filePrefixes))

peakAssymetryTValue = numeric(length(filePrefixes))
peakAssymetryPValue = numeric(length(filePrefixes))

valleyAssymetryTValue = numeric(length(filePrefixes))
valleyAssymetryPValue = numeric(length(filePrefixes))

for (i in 1:length(filePrefixes)) {
  
  ##### Parse Data #####
  
  # Parse the data into a normalized format
  normalizedData = parseBMHNucleosomeMutationData(
    paste0("Data/Raw Counts/",filePrefixes[i],"_nucleosome_mutation_counts.tsv"),
    paste0("Data/Background Data/",filePrefixes[i],"_nucleosome_mutation_background.tsv"))
  
  # Write the normalized data to a new file.
  fwrite(normalizedData, sep = '\t',
              file = paste0("Data/Normalized Counts/",filePrefixes[i],
                            "ESAD-UK_nucleosome_mutation_counts_normalized.txt"))
  
  ##### Periodicity Analysis #####
  
  # Calculate the periodicity of the data using a Lomb-Scargle periodiagram.
  lombResult = lsp(normalizedData[,.(Dyad_Position,Normalized_Both_Strands)], type = "period", plot = FALSE)
  plot(normalizedData[,.(Dyad_Position,Normalized_Both_Strands)],type = 'b', main = filePrefixes[i])
  # Store the relevant results!
  peakPeriodicities[i] = lombResult$peak.at[1]
  periodicityPValues[i] = lombResult$p.value
  
  # Calculate the SNR
  #timeSeries = as.ts(normalizedData[,.(Dyad_Position,Normalized_Both_Strands)])
  #test = bspec(timeSeries)
  #PSDEstimate = welchPSD(timeSeries, seglength = 10)
  #SNRResult = snr(timeSeries, PSDEstimate$power)
  
  
  ##### Assymetry Analysis #####
  
  # BEWARE THE TEXAS SHARPSHOOTER!
  
  # # Run assymetry analysis on the peaks and valleys of each strand.
  # runExtremeAnalysisSuite(normalizedData$Normalized_Plus_Strand, dyadPosCutoff = 68)
  # runExtremeAnalysisSuite(normalizedData$Normalized_Minus_Strand, dyadPosCutoff = 68)
  # runExtremeAnalysisSuite(normalizedData$Normalized_Plus_Strand, maxes = FALSE, dyadPosCutoff = 68)
  # runExtremeAnalysisSuite(normalizedData$Normalized_Minus_Strand, maxes = FALSE, dyadPosCutoff = 68)
  
  # Run assymetry analysis on the peaks and valleys of the aligned strands.
  peakAssymetryResult = 
    runExtremeAnalysisSuite(normalizedData$Normalized_Aligned_Strands, dyadPosCutoff = 68)
  peakAssymetryTValue[i] = peakAssymetryResult$statistic
  peakAssymetryPValue[i] = peakAssymetryResult$p.value
  
  
  valleyAssymetryResult = 
    runExtremeAnalysisSuite(normalizedData$Normalized_Aligned_Strands, maxes = FALSE, dyadPosCutoff = 68)
  valleyAssymetryTValue[i] = valleyAssymetryResult$statistic
  valleyAssymetryPValue[i] = valleyAssymetryResult$p.value
  
  # # Negative controls for assymetry analysis on peaks and valleys
  # runExtremeAnalysisSuite(normalizedData$Normalized_Both_Strands, dyadPosCutof = 68)
  # runExtremeAnalysisSuite(normalizedData$Normalized_Both_Strands, maxes = FALSE, dyadPosCutoff = 68)
  
}

# Create data.tables for all the results.
periodicityResults = data.table(Data_Set=filePrefixes,Peak_Periodicity=peakPeriodicities,PValue=periodicityPValues)

extremeAssymetryResults = data.table(Data_Set=filePrefixes, 
                                     Peak_Assymetry_TValue = peakAssymetryTValue, 
                                     Peak_Assymetry_PValue = peakAssymetryPValue,
                                     Valley_Assymetry_TValue = valleyAssymetryTValue, 
                                     Valley_Assymetry_PValue = valleyAssymetryPValue)

