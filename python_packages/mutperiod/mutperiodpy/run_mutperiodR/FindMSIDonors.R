library(mutperiodR)

selectInputAndRun = function(){
  
  # Select a file to use for MSIseq
  inputDataPath = choose.files(multi = FALSE, caption = "Select MSIseq Input Data",
                               filters = c(c("Tab Separated Files (*.tsv)","Any files"),
                                           c("*.tsv","*.*")), index = 1)
  
  write(as.character(findMSIDonors(inputDataPath)), "output.txt", sep = '\n')
  
}

# Get potential arguments from command line calls and use them to determine what action to take.
args = commandArgs(trailingOnly = T)
if (length(args) == 0) {
  selectInputAndRun()
} else if (length(args) == 2) {
  write(as.character(findMSIDonors(args[1])), args[2], sep = '\n')
} else {
  stop("Invalid number of arguments passed.  Expected 2 arguments for MSIseq input data path
       and output data path or no arguments to select input manually")
}