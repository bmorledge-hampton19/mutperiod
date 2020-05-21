library(NucPeriod)

selectInputAndRun = function(){
  
  # Select a file to use for MSIseq
  inputDataPaths = choose.files(multi = TRUE, caption = "Select raw counts files",
                               filters = c(c("Tab Separated Files (*.tsv)","Any files"),
                                           c("*.tsv","*.*")), index = 1)
  
  sapply(inputDataPaths, normalizeNucleosomeMutationCounts)
  
}


# Get potential arguments from command line calls and use them to determine what action to take.
args = commandArgs(trailingOnly = T)
if (length(args) == 1) {
  normalizeNucleosomeMutationCounts(args[1])
} else if (length(args) == 0) {
  selectInputAndRun()
} else {
  stop("Invalid number of arguments passed.  Expected 1 argument for raw mutation counts or no arguments
       to select input manually")
}