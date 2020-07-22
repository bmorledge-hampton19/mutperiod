##### Set Up the Environment #####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("deconstructSigs")

library(deconstructSigs)
library(data.table)

# Get the input data from the user.  The MSI donor list is optional.  If it is not selected, data will not be
# subsetted by microsatellite instability.
inputDataPath = choose.files(multi = FALSE, caption = "Select MutSig Input Data",
                             filters = c(c("Tab Separated Files (*.tsv)","Any files"),c("*.tsv","*.*")), index = 1)
MSIDonorListPath = choose.files(multi = FALSE, caption = "Select List of MSI Donors (Or cancel if no subsetting desired)",
                            filters = c(c("Text Files (*.txt)","Any files"),c("*.txt","*.*")), index = 1)

# Read in the input data and format it for deconstructSigs
inputData = fread(input = inputDataPath)
colnames(inputData) = c("Sample","chr","pos","ref","alt")

# If an MSI donor list was provided, also subset the data by microsatellite stability.
subsetByMSI = F
if (length(MSIDonorListPath) != 0) {
  subsetByMSI = T
  
  MSIDonorList = fread(input = MSIDonorListPath, header = F)
  colnames(MSIDonorList) = "MSI_Donors"
  
  MSIInputData = inputData[Sample %in% MSIDonorList$MSI_Donors]
  MSSInputData = inputData[!Sample %in% MSIDonorList$MSI_Donors]
     
}

# For now, lump all donors together into one group for each data set.
inputData$Sample = 1
if (subsetByMSI) {
  MSIInputData$Sample = 1
  MSSInputData$Sample = 1
}

# Deconstruct the signatures!
baseMutSigs = whichSignatures(mut.to.sigs.input(inputData), 1, contexts.needed = T)
makePie(baseMutSigs, sub = "No Omissions")
plotSignatures(baseMutSigs, sub = "No Omissions")

if (subsetByMSI) {
  MSIMutSigs = whichSignatures(mut.to.sigs.input(MSIInputData), 1, contexts.needed = T)
  makePie(MSIMutSigs, sub = "MSI")
  plotSignatures(MSIMutSigs, sub = "MSI")
  
  MSSMutSigs = whichSignatures(mut.to.sigs.input(MSSInputData), 1, contexts.needed = T)
  makePie(MSSMutSigs, sub = "MSS")
  plotSignatures(MSSMutSigs, sub = "MSS")
}