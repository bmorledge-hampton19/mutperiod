library(data.table)

# Creates figures for the given data and exports them to one or many pdf files.
# Files can be given in tsv form, or as an rda file containing a nucPeriodData object.
generateFigures = function(tsvFilePaths = list(), rdaFilePaths = list(), exportDir = getwd(), 
                           exportFileName = "plots.pdf", oneFile = TRUE, omitOutliers = TRUE,
                           dyadPosCutoff = 1000) {
  
  # Retrieve and name a list of data tables from any given tsv files.
  if (length(tsvFilePaths) > 0) {
    
    tsvDerivedTables = lapply(tsvFilePaths, fread)
    names(tsvDerivedTables) = lapply(tsvFilePaths, function(x) strsplit(basename(x), 
                                                                        "_nucleosome_mutation_counts")[[1]][1])
    
  } else tsvDerivedTables = list()
  
  # Retrieve the named list of nucleosome counts tables from the rda file.
  if (length(rdaFilePaths) > 0) {
    
    # Make a function for loading in the rda data, and returning the nucleosome counts tables
    retrieveRDATables = function(rdaFilePath) {
      
      if (load(rdaFilePath) != "nucPeriodData") {
        stop("rda file does not contain a nucPeriodData object.")
      }
      
      return(nucPeriodData$nucleosomeCountsTables)
      
    }
    
    # Apply the above function to the given file paths and combine them into one list.
    rdaDerivedTables = do.call(c, lapply(rdaFilePaths,retrieveRDATables))
    
  } else rdaDerivedTables = list()
 
  # Combine all the retrieved tables.
  countsTables = c(tsvDerivedTables, rdaDerivedTables)
  
  # Set up the export of the files for a single file if requested.
  if (oneFile) pdf(file = file.path(exportDir,exportFileName), width = 10.8)
  
  # Takes a single counts table and title and plots it!
  # If requested, opens up a new pdf stream and omits outliers
  plotCounts = function(countsTable, title) {
    
    # Open a new stream if the user requested multiple export files.
    if (!oneFile) pdf(file = file.path(exportDir,paste0(title,".pdf")), width = 10.8)
    
    # Get the relevant data from the counts table.
    dyadPositions = countsTable$Dyad_Position
    if ("Normalized_Both_Strands" %in% colnames(countsTable)) {
      counts = countsTable$Normalized_Both_Strands
    } else if ("Both_Strands_Counts" %in% colnames(countsTable)) {
      counts = countsTable$Both_Strands_Counts
    } else stop(paste("Counts table for",title,"does not have the expected data columns.",
                      "Expected a column titled \"Normalized_Both_Strands\" or \"Both_Strands_Counts\""))
    
    # Omit outliers if necessary
    if (omitOutliers) {
      omissionVector = !counts %in% boxplot(counts, plot = FALSE)$out
      dyadPositions = dyadPositions[omissionVector]
      counts = counts[omissionVector]
    }
      
    # Plot it!
    print(paste("Generating plot for",title))
    plot(dyadPositions, counts,
         type = 'b', main = title)
    
    # Make sure to close any open streams.
    if (!oneFile) dev.off()
    
  }
  
  # Pass the counts tables and their names to the plotting function
  captureOutput = mapply(plotCounts,countsTables,names(countsTables))
  
  # Close the "oneFile" pdf stream if its open.
  if (oneFile) dev.off()
   
}


# Get arguments from the command line which contain parameters for the generateFigures function.
# If there are no command line arguments, the above function is run with default parameters and the user is
# prompted to select rda files manually instead.
args = commandArgs(trailingOnly = T)

# If there are any inputs, there should be two, one with the path to the file with information on input and 
# output files and one to tell whether or not to omit outliers.
if (length(args) == 2) {
  
  # Read in inputs from the given file path.
  inputFile = file(args[1],'r')
  fileInputs = readLines(inputFile)
  close(inputFile)
  
  # Retrieve information for the function from the inputs file.
  if (length(fileInputs) == 4) {
    
    tsvFilePaths = strsplit(fileInputs[1],'$',fixed = TRUE)[[1]]
    rdaFilePaths = strsplit(fileInputs[2],'$',fixed = TRUE)[[1]]
    
    exportDir = fileInputs[3]
    
    # The presence of a file name on line 4 indicates that all graphs are to be exported to that one file.
    if (fileInputs[4] != '') {
      oneFile = TRUE
      exportFileName = fileInputs[4]
    } else {
      oneFile = FALSE
      exportFileName = NA
    }
    
  } else {
    stop(paste("Invalid number of arguments in input file.  Expected 4 argument for tsv file paths, rda file paths,",
          "output file directory, and output file name."))
  }
  
  generateFigures(tsvFilePaths, rdaFilePaths, exportDir, exportFileName, oneFile, args[2])
  
} else if (length(args) == 0) {
  generateFigures(rdaFilePaths = choose.files(multi = TRUE, caption = "Select NucPeriod Result files",
                                              filters = c(c("R Data Files (*.rda)","Any files"),
                                                          c("*.rda","*.*")), index = 1))
} else {
  stop(paste("Invalid number of command line arguments passed.  Expected 2 arguments for input data file path and", 
             "whether or not to omit outilers or no arguments to select input manually."))
}
