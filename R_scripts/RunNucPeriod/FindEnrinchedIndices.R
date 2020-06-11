library(data.table)
library(Biostrings)

# Gets the frequency of each base at a certain postion from the end of a set of sequences.
getBaseIndexFrequency = function(sequences, positionFromEnd) {
  
  bases = invisible(sapply(sequences, function(x) substr(x, nchar(x)-positionFromEnd+1, nchar(x)-positionFromEnd+1)))
  return(table(factor(bases, c("A","C","G","T"))))
  
}

# Prints the position of the maximum value of a given base and a vector of frequencies (ordered).
getMaxFrequencyPosition = function(base, frequencies) {
  
  maxValue = max(frequencies)
  maxPosition = match(maxValue,frequencies)
  print(paste0("Base ",base," is most common at position ",maxPosition," from the end with a count of ",maxValue,"."))
  
}


# Read in the fasta file
fastaFile = file.choose()
print(paste0("Reading from ",basename(fastaFile),"..."))
fastaTable = fread(file = fastaFile, header = FALSE, col.names = "Sequences")
fastaTable = fastaTable[!startsWith(Sequences,">")]

# Count the occurences of each base at a specific position from the end of each sequence.
baseFrequencies = data.table(bases = c("A","C","G","T"))
print("Calculating base frequencies at positions up to 12 from the end of each sequence...")
sapply(c(1:12), function(x) baseFrequencies[,paste0(x,"_From_End") := getBaseIndexFrequency(fastaTable$Sequences,x)])

# Print the maxima positions for each base.
invisible(sapply(baseFrequencies$bases, function(x) getMaxFrequencyPosition(x, baseFrequencies[bases == x,-1])))
