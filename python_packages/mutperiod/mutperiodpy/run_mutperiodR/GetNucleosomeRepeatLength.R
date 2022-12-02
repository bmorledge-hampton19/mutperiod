# This fairly simple script takes input from the command line to run a nucleosome map's counts against itself
# through an lsp to determine the peak periodicity (which should be the nucleosome repeat length) and record it
# at the specified output file path.
library(mutperiodR)
library(data.table)

args = commandArgs(trailingOnly = T)

if (length(args) == 2) {
  
  mutperiodData = generateMutperiodData(args[1], enforceInputNamingConventions = TRUE)
  lombResult = getNRL(args[1], returnFullLombResult = TRUE)
  NRL = lombResult$peak.at[1]
  noiseBooleanVector = (lombResult$scanned < NRL - 0.5
                        | lombResult$scanned > NRL + 0.5)
  SNR = lombResult$peak / median(lombResult$power[noiseBooleanVector])
  fwrite(list(NRL, paste("SNR:",SNR)), args[2], sep = '\n')

} else {
  stop("Invalid arguments.  Expected 2 arguments: The path to the nucleosome map self-counts, 
       and the output file path.")
}