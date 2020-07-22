# Set Up the Environment
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("MSIseq", quietly = TRUE))
  BiocManager::install("MSIseq")

suppressPackageStartupMessages(library(MSIseq))
suppressPackageStartupMessages(library(data.table))

selectInputAndRun = function(){
  
  # Select a file to use for MSIseq
  inputDataPath = choose.files(multi = FALSE, caption = "Select MSIseq Input Data",
                               filters = c(c("Tab Separated Files (*.tsv)","Any files"),
                                           c("*.tsv","*.*")), index = 1)
  
  findMSIDonors(inputDataPath)
  
}

findMSIDonors = function(inputDataPath) {
  
  # The number of base pairs (in megabases) in the hg19 genome (excluding mitochondria).
  captureLength = 3096
  
  ## download the Hg19repeats annotation file and load it
  print("Downloading genome information...")
  url <- "http://steverozen.net/data/Hg19repeats.rda" # NOTE: DON'T CHANGE THIS LINE WITHOUT CHANGING CAPTURE LENGTH
  file <- basename(url)
  if (!file.exists(file)) download.file(url, file)
  load(file)
  
  # Get the input data
  inputData = fread(inputDataPath)
  setnames(inputData,colnames(inputData),c("Chrom","Start_Position","End_Position",
                                           "Variant_Type","Tumor_Sample_Barcode"))
  
  # Remove mitochondrial mutations.
  inputData = inputData[Chrom != "chrM"]
  
  # Get mutation counts for input data.
  print("Running analysis...")
  mutationCounts = Compute.input.variables(inputData,Hg19repeats,captureLength)
  result = as.data.table(MSIseq.classify(mutationCounts))
  
  # Check for polE deficiency.
  PolEDeficientResults = result[Likely_POLE_deficiency == "Yes"]
  if (nrow(PolEDeficientResults) > 0) {
    print(paste(nrow(PolEDeficientResults),"PolE deficient result(s)."))
  }
  
  # Export the MSI-H results to a text file.
  MSIResults = result[MSI_status == "MSI-H",Tumor_Sample_Barcode]
  write(as.character(MSIResults),sep = '\n', file = paste0(dirname(inputDataPath),'/',
                                                           unlist(strsplit(basename(inputDataPath),'_'))[1],
                                                           "_MSI_donors.txt"))
}

# Get potential arguments from command line calls and use them to determine what action to take.
args = commandArgs(trailingOnly = T)
if (length(args) == 0) {
  selectInputAndRun()
} else if (length(args) == 1) {
  findMSIDonors(args[1])
} else {
  stop("Invalid number of arguments passed.  Expected 1 argument for MSIseq input data or no arguments
       to select input manually")
}