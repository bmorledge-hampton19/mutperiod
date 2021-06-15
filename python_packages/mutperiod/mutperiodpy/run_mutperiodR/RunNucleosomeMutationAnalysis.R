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
  if (length(inputs) == 4) {
    
    # Four inputs should mean that the input file contains a string of mutation counts file paths separated by '$'
    # followed by the path to the output file, whether or not to override the peak periodicity, and a string of
    # expected counts seperated by '$'.
    mutperiodData = generateMutperiodData(unlist(strsplit(inputs[1],'$',fixed = TRUE)),
                                          as.numeric(unlist(strsplit(inputs[4],'$',fixed = TRUE))), 
                                          enforceInputNamingConventions = TRUE, 
                                          relevantPeriodicityDefault = as.logical(inputs[3]))
    
  } else if (length(inputs) == 6) {
    
    # Five inputs should mean that the input file contains the counts file paths and the output file path
    # along with 2 file path groups for comparison and the default relevant periodicity (in that order).
    mutperiodData = generateMutperiodData(unlist(strsplit(inputs[1],'$',fixed = TRUE)),
                                          as.numeric(unlist(strsplit(inputs[6],'$',fixed = TRUE))),
                                          unlist(strsplit(inputs[3],'$',fixed = TRUE)),
                                          unlist(strsplit(inputs[4],'$',fixed = TRUE)),
                                          enforceInputNamingConventions = TRUE,
                                          relevantPeriodicityDefault = as.logical(inputs[5]))
    
  } else {
    stop("Invalid number of arguments in input file.  Expected 4 argument for mutation counts, output path, 
         and information on expected counts, or 6 arguments for mutation counts, output path,
         2 file path groups for comparison, and information on expected counts.")
  }
  
  writeResults(mutperiodData, inputs[2])
  
} else if (length(args) == 0) {
  selectInputAndRun()
} else {
  stop("Invalid number of command line arguments passed.  Expected 1 argument for input data file path, 
        or no arguments to select input manually")
}


