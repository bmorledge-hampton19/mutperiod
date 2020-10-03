library(NucPeriod)

selectInputAndRun = function(){
  
  # Select the raw counts file
  rawCountsFilePath = choose.files(multi = FALSE, caption = "Select raw counts files",
                                   filters = c(c("Tab Separated Files (*.tsv)","Any files"),
                                               c("*.tsv","*.*")), index = 1)
  
  # Select the background file
  backgroundCountsFilePath = choose.files(multi = FALSE, caption = "Select background counts files",
                                          filters = c(c("Tab Separated Files (*.tsv)","Any files"),
                                                      c("*.tsv","*.*")), index = 1)
  
  normalizedCounts = normalizeNucleosomeMutationCounts(rawCountsFilePath, 
                                                       backgroundCountsFilePath,
                                                       writeNormalizedData = FALSE)
  
}


# Get potential arguments from command line calls and use them to determine what action to take.
args = commandArgs(trailingOnly = T)
if (length(args) == 3) {
  normalizeNucleosomeMutationCounts(args[1],args[2],args[3])
} else if (length(args) == 0) {
  selectInputAndRun()
} else {
  stop(paste("Invalid number of arguments passed.  Expected 3 argument for raw mutation counts file,",
             "background counts file, and normalized counts file, or no arguments to select input manually"))
}