library(NucPeriod)

selectInputAndRun = function(){
  
  # Select a file to use for MSIseq
  inputDataPaths = choose.files(multi = TRUE, caption = "Select counts files",
                                filters = c(c("Tab Separated Files (*.tsv)","Any files"),
                                            c("*.tsv","*.*")), index = 1)
  
  generateNucPeriodData(inputDataPaths, outputFilePath = "nucMutationAnalysis.RData",
                        enforceInputNamingConventions = TRUE)
  
  print("Results outputted to nucMutationAnalysis.RData")
  
}


# Get potential arguments from command line calls and use them to determine what action to take.
args = commandArgs(trailingOnly = T)

if (length(args) == 2) {
  
  # Two arguments should mean that the script was passed a string of input file paths separated by '$'
  # followed by the path to the output file.
  generateNucPeriodData(unlist(strsplit(args[1],'$')), args[2], enforceInputNamingConventions = TRUE)
  
} else if (length(args) == 4) {
  
  # Four arguments should mean that script was passed the input file paths and the output file path
  # along with MSI and MSS designations (in that order).
  generateNucPeriodData(unlist(strsplit(args[1],'$')),args[2],
                        unlist(strsplit(args[3],'$')),unlist(strsplit(args[4],'$')))
  
} else if (length(args) == 0) {
  
  selectInputAndRun()

} else {
  stop("Invalid number of arguments passed.  Expected 2 argument for mutation counts and output path, 
        4 arguments for mutation counts, output path, and MSI/MSS designations,
        or no arguments to select input manually")
}