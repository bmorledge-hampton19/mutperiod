library(NucPeriod)

selectInputAndRun = function(){
  
  # Select a file to use for MSIseq
  inputDataPaths = choose.files(multi = TRUE, caption = "Select counts files",
                                filters = c(c("Tab Separated Files (*.tsv)","Any files"),
                                            c("*.tsv","*.*")), index = 1)
  
  generateNucPeriodData(inputDataPaths, enforceInputNamingConventions = TRUE)
  
}


# Get potential arguments from command line calls and use them to determine what action to take.
args = commandArgs(trailingOnly = T)

if (length(args) == 1) {
  
  # One argument should mean that the script was passed a string of file paths separated by '$'
  generateNucPeriodData(unlist(strsplit(args[1],'$')), enforceInputNamingConventions = TRUE)
  
} else if (length(args) == 3) {
  
  # 3 arguments should mean that script was passed the file paths along with MSI and MSS designations (in that order).
  generateNucPeriodData(unlist(strsplit(args[1],'$')),unlist(strsplit(args[2],'$')),unlist(strsplit(args[3],'$')))
  
} else if (length(args) == 0) {
  
  selectInputAndRun()

} else {
  stop("Invalid number of arguments passed.  Expected 1 argument for mutation counts, 3 arguments for
        mutation counts and MSI/MSS designations or no arguments to select input manually")
}