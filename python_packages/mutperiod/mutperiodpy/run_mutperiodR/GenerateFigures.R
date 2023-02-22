library(data.table)
library(ggplot2)

ROTATIONAL = 1
TRANSLATIONAL = 2

# Default text scaling
defaultTextScaling = theme(plot.title = element_text(size = 20, hjust = 0.5),
                           axis.title = element_text(size = 22), axis.text = element_text(size = 18),
                           legend.title = element_text(size = 22), legend.text = element_text(size = 18),
                           strip.text = element_text(size = 22))

# Blank background theme
blankBackground = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Actual domain colors
domainColors = c("BLACK" = "black", "black" = "black", "BLUE" = "blue", "blue" = "blue",
                 "GREEN" = "forestgreen", "green" = "forestgreen",
                 "RED" = "red", "red" = "red", "YELLOW" = "gold2", "yellow" = "gold2",
                 "GRAY" = "gray", "gray" = "gray")

# Creates figures for the given data and exports them to one or many pdf files.
# Files can be given in tsv form, or as an rda file containing a mutperiodData object.
generateFigures = function(tsvFilePaths = list(), tsvExpectedPeriods = list(), rdaFilePaths = list(), 
                           exportDir = getwd(), exportFileName = "plots.pdf", oneFile = TRUE, 
                           omitOutliers = TRUE, smoothNucGroup = TRUE, useAlignedStrands = FALSE) {
  
  ### Set up the positions for coloring periodic positions
  # Minor in vs. out positions from Cui and Zhurkin, with 1 added to each side and internal half-base positions included.
  minorInPositions = list(4:9-74, 14:19-74, 25:30-74, 36:41-74, 46:51-74, 56:61-74, 66:71-74,
                          5:9-74.5, 15:19-74.5, 26:30-74.5, 37:41-74.5, 47:51-74.5, 57:61-74.5, 67:71-74.5)
  minorOutPositions = list(9:14-74, 20:24-74, 31:35-74, 41:46-74 ,51:56-74, 61:66-74,
                           10:14-74.5, 21:24-74.5, 32:35-74.5, 42:46-74.5, 52:56-74.5, 62:66-74.5)
  
  # Retrieve and name a list of data tables from any given tsv files.
  if (length(tsvFilePaths) > 0) {
    
    tsvDerivedTables = lapply(tsvFilePaths, fread)
    names(tsvDerivedTables) = lapply(tsvFilePaths, function(x) strsplit(basename(x), 
                                                                        "_nucleosome_mutation_counts")[[1]][1])
    names(tsvExpectedPeriods) = names(tsvDerivedTables)
    
  } else tsvDerivedTables = list()
  
  # Retrieve the named list of nucleosome counts tables from the rda file.
  if (length(rdaFilePaths) > 0) {
    
    # Make a function for loading in the rda data, and returning the nucleosome counts tables
    retrieveRDATables = function(rdaFilePath) {
      
      if (load(rdaFilePath) != "mutperiodData") {
        stop(paste("rda file", basename(rdaFilePath), "does not contain a mutperiodData object."))
      }

      normTables = mutperiodData$normalizedNucleosomeCountsTables

      rawTablesInAnalysis = names(mutperiodData$rawNucleosomeCountsTables) %in% 
                              mutperiodData$periodicityResults$Data_Set
      rawTables = mutperiodData$rawNucleosomeCountsTables[rawTablesInAnalysis]
      
      return(c(normTables,rawTables))
      
    }
    
    # Apply the above function to the given file paths and combine them into one list.
    rdaDerivedTables = do.call(c, lapply(rdaFilePaths,retrieveRDATables))
    
    # Make another function for retrieving the expected periods from the given files...
    retrieveRDAExpectedPeriods = function(rdaFilePath) {
      
      load(rdaFilePath)
      return(setNames(mutperiodData$periodicityResults$Expected_Peak_Periodicity,
                      mutperiodData$periodicityResults$Data_Set))
      
    }
    
    # ...and apply it!
    rdaExpectedPeriods = do.call(c, lapply(rdaFilePaths,retrieveRDAExpectedPeriods))
    
  } else {
    rdaDerivedTables = list()
    rdaExpectedPeriods = list()
  }
 
  # Combine all the retrieved tables and expected periods
  countsTables = c(tsvDerivedTables, rdaDerivedTables)
  expectedPeriods = c(tsvExpectedPeriods, rdaExpectedPeriods)
  
  # Set up the export of the files for a single file if requested.
  if (oneFile) pdf(file = file.path(exportDir,exportFileName), width = 10.8)
  
  # Smoothing function for use in plotCounts later
  # Average across 11 base pairs centered on the given position.
  smoothValues = function(middlePos, data, dataCol, averagingRadius = 5) {

    positionsToAverage = (middlePos-averagingRadius):(middlePos+averagingRadius)
    valuesToAverage = data[Dyad_Position %in% positionsToAverage][[dataCol]]
    return(mean(valuesToAverage))
    
  }
  
  # Takes a single counts table and title and plots it!
  # If requested, opens up a new pdf stream and omits outliers
  plotCounts = function(dataSetName) {
    
    # Retrieve the counts table using the data set name.
    countsTable = countsTables[[dataSetName]]
    
    # Open a new stream if the user requested multiple export files.
    if (!oneFile) pdf(file = file.path(exportDir,paste0(dataSetName,".pdf")), width = 10.8)
    
    # Set up plot margins.
    par(mar = c(5,5,4,1))
    
    # Determine whether the data is rotational, rotational+linker, or translational.
    rotational = FALSE
    rotationalPlus = FALSE
    translational = FALSE
    if (min(countsTable$Dyad_Position) >= -73) {
      rotational = TRUE
    } else if (min(countsTable$Dyad_Position) > -999) {
      rotational = TRUE
      rotationalPlus = TRUE
    } else translational = TRUE
    
    # Get the relevant data from the counts table.
    if (!useAlignedStrands) {
      
      xlab = "Position Relative to Dyad (bp)"
    
      if ("Normalized_Both_Strands" %in% colnames(countsTable)) {
        dataCol = "Normalized_Both_Strands"
        ylab = "Normalized Counts"
      } else if ("Both_Strands_Counts" %in% colnames(countsTable)) {
        dataCol = "Both_Strands_Counts"
        ylab = "Raw Counts"
      } else stop(paste("Counts table for",dataSetName,"does not have the expected data columns.",
                        "Expected a column titled \"Normalized_Both_Strands\" or \"Both_Strands_Counts\""))
      
    } else {
      
      xlab = "Position Relative to Dyad (bp, strand aligned)"
      
      if ("Normalized_Aligned_Strands" %in% colnames(countsTable)) {
        dataCol = "Normalized_Aligned_Strands"
        ylab = "Normalized Counts"
      } else if ("Aligned_Strands_Counts" %in% colnames(countsTable)) {
        dataCol = "Aligned_Strands_Counts"
        ylab = "Raw Counts"
      } else stop(paste("Counts table for",dataSetName,"does not have the expected data columns.",
                        "Expected a column titled \"Normalized_Aligned_Strands\" or \"Aligned_Strands_Counts\""))
      
    }
    
    if (omitOutliers && smoothNucGroup && translational) {
      warning(paste0("Combining outlier filtering with smoothing may have undesired consequences as the ",
                     "window for averaging is not adjusted due to omitted outliers."))
    }
    
    # Omit outliers if necessary
    if (omitOutliers) {
      omissionVector = !countsTable[[dataCol]] %in% boxplot(countsTable[[dataCol]], plot = FALSE)$out
      countsTable = countsTable[omissionVector]
    }
    
    # If smoothing is requested for nuc group tables, and we have such a table, smooth away!
    if (translational && smoothNucGroup) {
      countsTable[, (dataCol) := sapply(countsTable$Dyad_Position, smoothValues, 
                                        data = countsTable, dataCol = dataCol)]
    }
    
    if (rotational) {
      # Color rotational positioning
      countsTable[abs(Dyad_Position) %in% abs(unlist(minorInPositions)),
                 Periodic_Position_Color := "cyan3"]
      countsTable[abs(Dyad_Position) %in% abs(unlist(minorOutPositions)),
                 Periodic_Position_Color := "darkorange"]
      countsTable[is.na(Periodic_Position_Color), Periodic_Position_Color := "Black"]
      
    }
    
    if (rotationalPlus) {
      # Color linker DNA in linker+ plots.
      countsTable[Dyad_Position <= -73 | Dyad_Position >= 73, Periodic_Position_Color := "deeppink2"]
      
    }
    
    if (translational) {
      # For translational coloring, get the linker and nucleosome positions based on the expected period.
      nucRepLen = round(expectedPeriods[[dataSetName]])
      
      nucleosomePositions = sapply(1:10, function(x) return(append((-73+x*nucRepLen):(73+x*nucRepLen),
                                                                   (-72.5+x*nucRepLen):(72.5+x*nucRepLen))))
      nucleosomePositions = append(nucleosomePositions, c(0:73, 0.5:72.5))
      linkerPositions = sapply(0:8, function(x) return( append((74+x*nucRepLen):(-74+(x+1)*nucRepLen),
                                                               (73.5+x*nucRepLen):(-73.5+(x+1)*nucRepLen))))
      
      # Color translational positioning
      countsTable[Dyad_Position %in% nucleosomePositions | -Dyad_Position %in% nucleosomePositions,
                  Periodic_Position_Color := "darkblue"]
      countsTable[Dyad_Position %in% linkerPositions | -Dyad_Position %in% linkerPositions,
                  Periodic_Position_Color := "deeppink2"]
      
    }
    
    # Plot it!
    print(paste("Generating plot for",dataSetName))
    
    periodicityPlot = ggplot(countsTable, aes_string("Dyad_Position", dataCol,
                                                    color = "Periodic_Position_Color"))
    if (packageVersion("ggplot2") >= 3.4) {
      periodicityPlot = periodicityPlot + geom_path(linewidth = 1.25, lineend = "round", aes(group = 1))
    } else {
      periodicityPlot = periodicityPlot + geom_path(size = 1.25, lineend = "round", aes(group = 1))
    }
    
    if (rotationalPlus) {
      periodicityPlot = periodicityPlot +
        scale_color_identity(name = '', guide = "legend",
                             breaks = c("cyan3", "darkorange", "deeppink2"),
                             #breaks = c("#1bcc44", "#993299", "Gold"),
                             labels = c("Minor-in", "Minor-out", "Linker"))
    } else if (rotational) {
      periodicityPlot = periodicityPlot +
        scale_color_identity(name = '', guide = "legend",
                             breaks = c("cyan3", "darkorange"),
                             #breaks = c("#1bcc44", "#993299"),
                             labels = c("Minor-in", "Minor-out"))
    } else if (translational) {
      periodicityPlot = periodicityPlot +
        scale_color_identity(name = '', guide = "legend",
                             breaks = c("darkblue", "deeppink2"),
                             #breaks = c("#0571b0", "#ca0020"),
                             labels = c("Nucleosome","Linker"))
    }
  
  periodicityPlot = periodicityPlot +
    labs(title = dataSetName, x = xlab, y = ylab) +
    defaultTextScaling + blankBackground
  
  print(periodicityPlot)
    
    # Make sure to close any open streams.
    if (!oneFile) dev.off()
    
  }
  
  # Pass the counts tables and their names to the plotting function
  captureOutput = sapply(sort(names(countsTables)), plotCounts)
  
  # Close the "oneFile" pdf stream if its open.
  if (oneFile) captureOutput = dev.off()
   
}


# Get arguments from the command line which contain parameters for the generateFigures function.
# If there are no command line arguments, the above function is run with default parameters and the user is
# prompted to select rda files manually instead.
args = commandArgs(trailingOnly = T)

# If there are any inputs, there should be four total, one with the path to the file with information on input and 
# output files and three to set other parameters of the graphing process.
if (length(args) == 4) {
  
  # Read in inputs from the given file path.
  inputFile = file(args[1],'r')
  fileInputs = readLines(inputFile)
  close(inputFile)
  
  # Retrieve information for the function from the inputs file.
  if (length(fileInputs) == 5) {
    
    tsvFilePaths = strsplit(fileInputs[1],'$',fixed = TRUE)[[1]]
    tsvExpectedPeriods = as.numeric(strsplit(fileInputs[2],'$',fixed = TRUE)[[1]])
    rdaFilePaths = strsplit(fileInputs[3],'$',fixed = TRUE)[[1]]
    
    exportDir = fileInputs[4]
    
    # The presence of a file name on line 5 indicates that all graphs are to be exported to that one file.
    if (fileInputs[5] != '') {
      oneFile = TRUE
      exportFileName = fileInputs[5]
    } else {
      oneFile = FALSE
      exportFileName = NA
    }
    
  } else {
    stop(paste("Invalid number of arguments in input file.  Expected 5 argument for tsv file paths,",
          "expected tsv periods, rda file paths, output file directory, and output file name."))
  }
  
  generateFigures(tsvFilePaths, tsvExpectedPeriods, rdaFilePaths, exportDir, exportFileName, oneFile, 
                  as.logical(args[2]), as.logical(args[3]), as.logical(args[4]))
  
} else if (length(args) == 0) {
  generateFigures(rdaFilePaths = choose.files(multi = TRUE, caption = "Select mutperiod Result files",
                                              filters = c(c("R Data Files (*.rda)","Any files"),
                                                          c("*.rda","*.*")), index = 1))
} else {
  stop(paste("Invalid number of command line arguments passed.  Expected 6 arguments for input data file path and", 
             "5 other parameters of the plot function (Starting with omitOutliers)."))
}
