library(mutperiodR)
library(data.table)

selectInputAndRun = function(){
  
  # Get the input data from the user.
  inputDataPath = choose.files(multi = FALSE, caption = "Select MutSig Input Data", index = 1,
                               filters = c(c("Tab Separated Files (*.tsv)","Any files"),c("*.tsv","*.*")))
  
  fwrite(getMutSigs(inputDataPath), "output.tsv", sep = '\t')
  
}

# Get potential arguments from command line calls and use them to determine what action to take.
args = commandArgs(trailingOnly = T)
if (length(args) == 0) {
  selectInputAndRun()
} else if (length(args) == 2) {
  fwrite(getMutSigs(args[1]), args[2], sep = '\t')
} else {
  stop("Invalid number of arguments passed.  Expected 2 arguments for deconstructSigs input data path
       and output data path or no arguments to select input manually")
}




