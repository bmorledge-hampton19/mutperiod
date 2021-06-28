library(mutperiodR)

writeResults = function(mutperiodData, outputFilePath) {
  
  fileExtension = substr(outputFilePath, nchar(outputFilePath)-3, nchar(outputFilePath))
  if (fileExtension == ".rda") {
    save(mutperiodData, file = outputFilePath)
  } else if (fileExtension == ".tsv") {
    data.table::fwrite(mutperiodData$periodicityResults, file = outputFilePath, sep = '\t')
  } else {
    stop(paste0("Invalid file extension for output file.  Should be .tsv or .rda but got: ",fileExtension))
  }
  
}

selectInputAndRun = function(){
  
  # Select a file to use for MSIseq
  inputDataPaths = choose.files(multi = TRUE, caption = "Select counts files",
                                filters = c(c("Tab Separated Files (*.tsv)","Any files"),
                                            c("*.tsv","*.*")), index = 1)
  
  mutperiodData = generateMutperiodData(inputDataPaths, outputFilePath = "nucMutationAnalysis.RData",
                                        enforceInputNamingConventions = TRUE)
  
  save(mutperiodData, file = "nucMutationAnalysis.rda")
  
  print("Results outputted to nucMutationAnalysis.rda")
  
}


# Gets a file path from command line which contains inputs for the generateMutperiodData function.
# If there are no command line arguments, the above function is run instead.
args = commandArgs(trailingOnly = T)

if (length(args) == 1) {
  
  # Read in inputs from the given file path.
  inputFile = file(args[1],'r')
  inputs = readLines(inputFile)
  close(inputFile)
  
  # Call the function to generate the mutperiodData object based on the inputs in the given file.
  if (length(inputs) == 5) {
    
    # Five inputs runs the analysis without grouped comparison.
    mutperiodData = generateMutperiodData(unlist(strsplit(inputs[1],'$',fixed = TRUE)),
                                          expectedPeakPeriodicities = 
                                            as.numeric(unlist(strsplit(inputs[4],'$',fixed = TRUE))), 
                                          enforceInputNamingConventions = TRUE, 
                                          overridePeakPeriodicityWithExpectedPeak = as.logical(inputs[3]),
                                          alignStrands = as.logical(inputs[5]))
    
  } else if (length(inputs) == 7) {
    
    # Seven inputs runs the analysis with grouped comparison.
    mutperiodData = generateMutperiodData(unlist(strsplit(inputs[1],'$',fixed = TRUE)),
                                          expectedPeakPeriodicities = 
                                            as.numeric(unlist(strsplit(inputs[6],'$',fixed = TRUE))), 
                                          filePathGroup1 = unlist(strsplit(inputs[3],'$',fixed = TRUE)),
                                          filePathGroup2 = (strsplit(inputs[4],'$',fixed = TRUE)),
                                          enforceInputNamingConventions = TRUE,
                                          overridePeakPeriodicityWithExpectedPeak = as.logical(inputs[5]),
                                          alignStrands = as.logical(inputs[7]))
    
  } else {
    stop("Invalid number of arguments in input file.  Expected 5 argument for analysis without grouped comparison, 
         or 7 for analysis with grouped comparison.")
  }
  
  writeResults(mutperiodData, inputs[2])
  
} else if (length(args) == 0) {
  selectInputAndRun()
} else {
  stop("Invalid number of command line arguments passed.  Expected 1 argument for input data file path, 
        or no arguments to select input manually")
}


