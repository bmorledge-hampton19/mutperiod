#Ensure that BiocManager is installed
if (!requireNamespace("BiocManager", quietly = T))
  install.packages("BiocManager")
if (!requireNamespace("lomb", quietly = T))
  BiocManager::install("lomb")
if (!requireNamespace("quantmod", quietly = T))
  BiocManager::install("quantmod")
if (!requireNamespace("bspec", quietly = T))
  BiocManager::install("bspec")

library(lomb)
library(data.table)
library(bspec)

source("Normalize_Nucleosome_Mutation_Counts_Module.R")
source("Asymmetry_Analysis.R")

rawCountsFilePaths = choose.files(caption = "Select Raw Nucleosome Mutation Data.",
                                  filters = c(c("Tab Separated Values (*.tsv)","Any files"),
                                              c("*.tsv","*.*")), index = 1)
dataDir = dirname(dirname(rawCountsFilePaths[1]))

# Generate a list of prefixes for the data files so that iterating through them is easier later.
filePrefixes = sapply(strsplit(basename(rawCountsFilePaths),"_nucleosome"), function(x) x[1])

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
  
  ##### Parse and normalize data #####
  
  normalizeNucleosomeMutationCounts(rawCountsFilePaths[i])
  
  ##### Periodicity Analysis #####
  
  # Calculate the periodicity of the data using a Lomb-Scargle periodiagram.
  lombResult = lsp(normalizedData[,.(Dyad_Position,Normalized_Both_Strands)], 
                   type = "period", from = 2, to = 50, ofac = 100, plot = F)
  plot(normalizedData[,.(Dyad_Position,Normalized_Both_Strands)],type = 'b', main = filePrefixes[i])
  # Store the relevant results!
  peakPeriodicities[i] = lombResult$peak.at[1]
  periodicityPValues[i] = lombResult$p.value
  
  # Calculate the SNR
  noiseBooleanVector = (lombResult$scanned < lombResult$peak.at[1] - 0.5
                        | lombResult$scanned > lombResult$peak.at[1] + 0.5)
  periodicitySNRs[i] = lombResult$peak / median(lombResult$power[noiseBooleanVector])
  
  # Attempted to use bspec package...
  #timeSeries = as.ts(normalizedData[,.(Normalized_Both_Strands)],start = -73, end = 73)
  #PSDEstimate = welchPSD(timeSeries, seglength = 10)
  #SNRResult = snr(timeSeries, PSDEstimate$power)
  
  
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
periodicityResults = data.table(Data_Set=filePrefixes,Peak_Periodicity=peakPeriodicities,
                                PValue=periodicityPValues,SNR=periodicitySNRs)

asymmetryResults = data.table(Data_Set=filePrefixes, 
                                     General_Asymmetry_TValue = generalAsymmetryTValue,
                                     General_Asymmetry_PValue = generalAsymmetryPValue,
                                     Peak_Asymmetry_TValue = peakAsymmetryTValue, 
                                     Peak_Asymmetry_PValue = peakAsymmetryPValue,
                                     Valley_Asymmetry_TValue = valleyAsymmetryTValue, 
                                     Valley_Asymmetry_PValue = valleyAsymmetryPValue)
