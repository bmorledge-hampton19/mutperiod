library(nucperiodR)

selectInputAndRun = function(){
  
  # Select a file to use for MSIseq
  inputDataPaths = choose.files(multi = TRUE, caption = "Select counts files",
                                filters = c(c("Tab Separated Files (*.tsv)","Any files"),
                                            c("*.tsv","*.*")), index = 1)
  
  generateNucPeriodData(inputDataPaths, outputFilePath = "nucMutationAnalysis.RData",
                        enforceInputNamingConventions = TRUE)
  
  print("Results outputted to nucMutationAnalysis.RData")
  
}


# Gets a file path from command line which contains inputs for the generateNucPeriodData function.
# If there are no command line arguments, the above function is run instead.
args = commandArgs(trailingOnly = T)

if (length(args) == 1) {
  
  # Read in inputs from the given file path.
  inputFile = file(args[1],'r')
  inputs = readLines(inputFile)
  close(inputFile)
  
  # Call the function to generate the NucPeriodData object based on the inputs in the given file.
  if (length(inputs) == 2) {
    
    # Two inputs should mean that the input file contains a string of mutation counts file paths separated by '$'
    # followed by the path to the output file.
    generateNucPeriodData(unlist(strsplit(inputs[1],'$',fixed = TRUE)), inputs[2], 
                          enforceInputNamingConventions = TRUE)
    
  } else if (length(inputs) == 4) {
    
    # Two inputs should mean that the input file contains the counts file paths and the output file path
    # along with MSI and MSS designations (in that order).
    generateNucPeriodData(unlist(strsplit(inputs[1],'$',fixed = TRUE)),inputs[2],
                          unlist(strsplit(inputs[3],'$',fixed = TRUE)),
                          unlist(strsplit(inputs[4],'$',fixed = TRUE)),
                          enforceInputNamingConventions = TRUE)
    
  } else {
    stop("Invalid number of arguments in input file.  Expected 2 argument for mutation counts and output path, or
          4 arguments for mutation counts, output path, and MSI/MSS designations.")
  }
  
} else if (length(args) == 0) {
  selectInputAndRun()
} else {
  stop("Invalid number of command line arguments passed.  Expected 1 argument for input data file path, 
        or no arguments to select input manually")
}

  