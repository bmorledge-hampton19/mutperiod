library(mutperiodR)

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
# Three arguments represents the raw counts file path, background counts file path, and output file path.
# Four represents the same as three but with an alternative scaling factor added in.
if (length(args) == 3) {
  normalizeNucleosomeMutationCounts(args[1],args[2],args[3])
} else if (length(args) == 4) {
  normalizeNucleosomeMutationCounts(args[1],args[2],args[3], as.numeric(args[4]))
} else if (length(args) == 0) {
  selectInputAndRun()
} else {
  stop(paste("Invalid number of arguments passed.  Expected 3 argument for raw mutation counts file,",
             "background counts file, and normalized counts file, 4 arguments if an alternative scaling",
             "factor is included, or no arguments to select input manually"))
}