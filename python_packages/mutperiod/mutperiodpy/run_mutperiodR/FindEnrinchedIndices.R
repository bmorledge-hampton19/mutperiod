library(data.table)
library(tools)
library(Biostrings)

#choice = menu(c("Single Nucleotides", "Dipys"), title = "Search for...")
choice = 2

if (choice == 1) {
  dipys = FALSE
} else if (choice == 2) {
  dipys = TRUE
} else {
  stop("User chose to exit.")
}


# Calculates the base frequencies for a given set of sequences, at given values, 
# and return the resulting baseFrequency table.
calculateBaseFrequencyTable = function(sequences, fromStartValues, fromEndValues, dipys = FALSE) {
  
  # Do some setup based on whether or not this function is counting individual bases or dipys.
  if (!dipys) {
    checkBases = c("A","C","G","T")
    dipyOffset = 0
  } else {
    checkBases = c("CC","CT","TT","TC")
    dipyOffset = 0.5
  }
  
  
  # Get base frequencies by position for the "fromStartValues"
  if (length(fromStartValues) > 0) {
    
    # First, get raw bases by position.
    basesByLocationFromStart = sapply(sequences, function(sequence) { 
      return(sapply(fromStartValues, function(fromStartValue) {
        
        if(!dipys) {
          return(substr(sequence, fromStartValue, fromStartValue))
        } else {
          return(substr(sequence, fromStartValue, fromStartValue+1))
        }
        
        
      }))
    })
    
    # Then, convert to frequencies.
    baseFrequenciesFromStart = apply(basesByLocationFromStart, 1, function(bases) {table(factor(bases, checkBases))})
    if (dipys) {
      baseFrequenciesFromEnd = t(as.matrix(colSums(baseFrequenciesFromStart)))
      rownames(baseFrequenciesFromEnd) = "Dipys"
    }
    baseFrequenciesFromStart = data.table(baseFrequenciesFromStart, keep.rownames = TRUE)
    
    # Add column names as necessary.
    colnames(baseFrequenciesFromStart) = c("Bases",sapply(fromStartValues, function(fromStartValue) {
      paste0(fromStartValue + dipyOffset,"_From_Start")
    }))
    
  # If there are no "fromStart" values, at least set up the first column.
  } else {
    if (!dipys) {
      baseFrequenciesFromStart = data.table(Bases = checkBases)
    } else {
      baseFrequenciesFromStart = data.table(Bases = "Dipys")
    }
    
  }
  
  
  # Do the same (mostly) for the "fromEndValues"
  if (length(fromEndValues) > 0) {
    
    basesByLocationFromEnd = sapply(sequences, function(sequence) { 
      sequenceLength = nchar(sequence)
      return(sapply(fromEndValues, function(fromEndValue) {
        
        if(!dipys) {
          return(substr(sequence, sequenceLength-fromEndValue+1, sequenceLength-fromEndValue+1))
        } else {
          return(substr(sequence, sequenceLength-fromEndValue, sequenceLength-fromEndValue+1))
        }
        
        
      }))
    })
    baseFrequenciesFromEnd = apply(basesByLocationFromEnd, 1, function(bases) {table(factor(bases, checkBases))})
    if (dipys) {
      baseFrequenciesFromEnd = t(as.matrix(colSums(baseFrequenciesFromEnd)))
    }
    baseFrequenciesFromEnd = data.table(baseFrequenciesFromEnd)
    
    colnames(baseFrequenciesFromEnd) = sapply(fromEndValues, function(fromEndValue) {
      paste0(fromEndValue + dipyOffset,"_From_End")
    })
   
    # One major difference here:  Bind the two tables together.
    baseFrequencyTable = cbind(baseFrequenciesFromStart, baseFrequenciesFromEnd)
    
  # If there are no "fromEndValues" Just make the fromStart table the main table!
  } else {
    baseFrequencyTable = baseFrequenciesFromStart
  }
  
  return(baseFrequencyTable)
  
}

# Returns, as text, the position of the maximum frequency value of a given base in a base frequency table
# followed by the frequency itself.
getMaxFrequencyInfo = function(thisBase, frequencyTable) {
  
  relevantFrequencies = frequencyTable[Bases == thisBase, -1]
  maxValue = max(relevantFrequencies)
  maxFrequencyLocation = colnames(relevantFrequencies)[match(maxValue,relevantFrequencies)]
  
  return(list(thisBase, maxValue, maxFrequencyLocation))
  
}


sequenceLengthRange = 10:40

# Read in the fasta file
fastaFile = file.choose()
startTime = Sys.time()
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
  if (length(sequencesByLength[[as.character(i)]]) > 0) {
    baseFrequencyTables[[as.character(i)]] = 
      calculateBaseFrequencyTable(sequencesByLength[[as.character(i)]], 
                                  fromStartValues = numeric(), fromEndValues = 1:10, dipys = dipys)
  }
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
    if (length(baseFrequencyTables[[as.character(i)]]) > 0) {
      maxFrequencyInfo = getMaxFrequencyInfo(feature, baseFrequencyTables[[as.character(i)]])
      maxValues[as.character(i)] = maxFrequencyInfo[[2]]
      maxValueLocations[as.character(i)] = maxFrequencyInfo[[3]]
      
      print(paste0("For fragments of length ",i, 
                   ", maximum frequency of ", feature, " is at position ", 
                   maxValueLocations[as.character(i)], " with a frequency of ", 
                   maxValues[as.character(i)], '.'))
    }
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
