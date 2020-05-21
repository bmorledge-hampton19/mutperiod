library(NucPeriod)

selectInputAndRun = function(){
  
  # Select a file to use for MSIseq
  inputDataPath = choose.files(multi = FALSE, caption = "Select MSIseq Input Data",
                               filters = c(c("Tab Separated Files (*.tsv)","Any files"),
                                           c("*.tsv","*.*")), index = 1)
  
  findMSIDonors(inputDataPath)
  
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