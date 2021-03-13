library(data.table)
library(tools)
library(Biostrings)

#choice = menu(c("Single Nucleotides", "Dipys"), title = "Search for...")
#What is happening?
choice = 2

if (choice == 1) {
  dipys = FALSE
} else if (choice == 2) {
  dipys = TRUE
} else {
  stop("User chose to exit.")
}

startTime = Sys.time()

# Gets the frequency of each base at a certain postion from the start or end of a set of sequences.
# Exactly one of the "position" arguments should be provided.  Otherwise, an error is raised.
getBaseIndexFrequency = function(sequences, positionFromEnd = NA, positionFromStart = NA, dipys = FALSE) {
  
  if ( (is.na(positionFromEnd) && is.na(positionFromStart)) || 
       (!is.na(positionFromEnd) && !is.na(positionFromStart)) ) {
    stop("Error: Exactly one position argument should be given, not both or neither.")
  }
  
  if (!is.na(positionFromEnd)) {
    bases = sapply(sequences, function(sequence) {
      
      sequenceLength = nchar(sequence)
      if (sequenceLength < positionFromEnd) {
        stop(paste0("Desired position (",positionFromEnd," from end) ",
                    "exceeds sequence length (",sequenceLength,")"))
      }
      
      if (!dipys) {       
        return(substr(sequence, sequenceLength-positionFromEnd+1, sequenceLength-positionFromEnd+1))
      } else {
        return(substr(sequence, sequenceLength-positionFromEnd, sequenceLength-positionFromEnd+1))
      }
      
    })
    
  } else if (!is.na(positionFromStart)) {
    bases = sapply(sequences, function(sequence) {
      
      sequenceLength = nchar(sequence)
      if (sequenceLength < positionFromStart) {
        stop(paste0("Desired position (",positionFromStart," from start) ",
                    "exceeds sequence length (",sequenceLength,")"))
      }
      
      if (!dipys) { 
        return(substr(sequence, positionFromStart, positionFromStart))
      } else {
        return(substr(sequence, positionFromStart, positionFromStart+1))
      }
      
    })
  }
  
  if (!dipys) {
    countsTable = table(factor(bases, c("A","C","G","T")))
    return(countsTable/length(bases))
  } else {
    countsTable = table(factor(bases, c("CC","CT","TT","TC")))
    return(sum(countsTable)/length(bases))
  }
  
}

# Calculates the base frequencies for a given set of sequences, at given values, 
# and return the resulting baseFrequency table.
calculateBaseFrequencyTable = function(sequences, fromStartValues, fromEndValues, dipys = FALSE) {
  
  if (!dipys) {
    baseFrequencies = data.table(base = c("A","C","G","T"))
  } else {
    baseFrequencies = data.table(base = "Dipys")
  }
  
  for (i in fromStartValues) { 
    print(paste0("Calculating base frequency at position ",i," base(s) from start."))
    if (!dipys) {
      featurePos = i
    } else {
      featurePos = i + 0.5
    }
    baseFrequencies[,paste0(featurePos,"_From_Start") := getBaseIndexFrequency(sequences, positionFromStart = i, 
                                                                      dipys = dipys)]
  }
    
  
  for (i in fromEndValues) { 
    print(paste0("Calculating base frequency at position ",i," base(s) from end."))
    if (!dipys) {
      featurePos = i
    } else {
      featurePos = i + 0.5
    }
    baseFrequencies[,paste0(featurePos,"_From_End") := getBaseIndexFrequency(sequences, positionFromEnd = i, 
                                                                             dipys = dipys)]
  }
  
  return(baseFrequencies)
  
}

# Returns, as text, the position of the maximum frequency value of a given base in a base frequency table
# followed by the frequency itself.
getMaxFrequencyInfo = function(thisBase, frequencyTable) {
  
  relevantFrequencies = frequencyTable[base == thisBase, -1]
  maxValue = max(relevantFrequencies)
  maxFrequencyLocation = colnames(relevantFrequencies)[match(maxValue,relevantFrequencies)]
  
  return(list(thisBase, maxValue, maxFrequencyLocation))
  
}


sequenceLengthRange = 10:40

# Read in the fasta file
fastaFile = file.choose()
print(paste0("Reading from ",basename(fastaFile),"..."))
fastaTable = fread(file = fastaFile, header = FALSE, col.names = "Sequences")
fastaTable = fastaTable[!startsWith(Sequences,">")]

# Split the fasta sequences based on length
print("Splitting fasta sequences by length...")
fastaTableSequenceLengths = nchar(fastaTable$Sequences)

sequencesByLength = sapply(sequenceLengthRange, function(i) {
  print(paste0("Getting sequences of Length ",i,"..."))
  return(fastaTable$Sequences[fastaTableSequenceLengths == i])
})
names(sequencesByLength) = as.character(sequenceLengthRange)

# Calculate a base frequency table for each sequence length.
print("Calculating base frequency table for each sequence length...")
baseFrequencyTables = vector("list", length(sequenceLengthRange))
names(baseFrequencyTables) = as.character(sequenceLengthRange)
for (i in sequenceLengthRange) {
  print(paste0("Calculating base frequency table for sequence length of ",i,"..."))
  baseFrequencyTables[[as.character(i)]] = 
    calculateBaseFrequencyTable(sequencesByLength[[as.character(i)]], 
                                fromStartValues = numeric(), fromEndValues = 1:10, dipys = dipys)
}

# Print the position of the maximum frequency value for each of the bases in each frequency table.
maxValues = numeric(length(sequenceLengthRange))
names(maxValues) = as.character(sequenceLengthRange)
maxValueLocations = character(length(sequenceLengthRange))
names(maxValueLocations) = as.character(sequenceLengthRange)

maxFrequencyInfoTable = data.table(Sequence_Length = sequenceLengthRange)

if (!dipys) {
  features = c('A','C','T','G')
} else {
  features = c("Dipys")
}
for (feature in features) {
  
  for (i in sequenceLengthRange) {
    
    maxFrequencyInfo = getMaxFrequencyInfo(feature, baseFrequencyTables[[as.character(i)]])
    maxValues[as.character(i)] = maxFrequencyInfo[[2]]
    maxValueLocations[as.character(i)] = maxFrequencyInfo[[3]]
    
    print(paste0("For fragments of length ",i, 
                 ", maximum frequency of ", feature, " is at position ", 
                 maxValueLocations[as.character(i)], " with a frequency of ", 
                 maxValues[as.character(i)], '.'))
  }
  
  maxFrequencyInfoTable[,paste0("Max_",feature,"_Frequency") := maxValues]
  maxFrequencyInfoTable[,paste0("Max_",feature,"_Frequency_Location") := maxValueLocations]
  
}

stopTime = Sys.time()
print(stopTime - startTime)

if (!dipys) {
  outputFilePath = file = file.path(dirname(fastaFile),
                                    paste0(file_path_sans_ext(basename(fastaFile)),"_max_base_frequency.csv"))
} else {
  outputFilePath = file = file.path(dirname(fastaFile),
                                    paste0(file_path_sans_ext(basename(fastaFile)),"_max_dipy_frequency.csv"))
}

fwrite(maxFrequencyInfoTable, file = outputFilePath)
