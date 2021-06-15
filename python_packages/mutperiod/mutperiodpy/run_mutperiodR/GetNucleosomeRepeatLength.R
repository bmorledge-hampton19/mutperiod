# This fairly simple script takes input from the command line to run a nucleosome map's counts against itself
# through an lsp to determine the peak periodicity (which should be the nucleosome repeat length) and record it
# at the specified output file path.
library(mutperiodR)

args = commandArgs(trailingOnly = T)

if (length(args) == 2) {
  
  mutperiodData = generateMutperiodData(args[1], enforceInputNamingConventions = TRUE)
  fwrite(mutperiodData$periodicityResults$Relevant_Periodicity[1], args[2])

} else {
  stop("Invalid arguments.  Expected 2 arguments: The path to the nucleosome map self-counts, 
       and the output file path.")
}