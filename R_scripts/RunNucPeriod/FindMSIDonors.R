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
} else if (length(args) == 2) {
  findMSIDonors(args[1], args[2])
} else {
  stop("Invalid number of arguments passed.  Expected 2 arguments for MSIseq input data 
       and mutation group name or no arguments to select input manually")
}