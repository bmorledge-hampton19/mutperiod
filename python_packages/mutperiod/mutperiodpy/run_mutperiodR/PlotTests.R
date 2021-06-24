# Don't forget to import the data.table library!
library(data.table)
library(ggplot2)

# Read in data.
load(choose.files(multi = FALSE))
dataSetNames = as.list(setNames(nm = mutperiodData$periodicityResults$Data_Set))

# Column names...
# ...for normalized (CHOOSE ONE)
dataCol = "Normalized_Both_Strands"

dataCol = "Normalized_Aligned_Strands"

# ...for raw (CHOOSE ONE)
dataCol = "Both_Strands_Counts"

dataCol = "Aligned_Strands_Counts"

# Where to trim rotational only plots to (default of 60 is what is used in the lsp analysis)
rotationalOnlyCutoff = 60

# Default plot values
title = "PLACEHOLDER TITLE"
yAxisLabel = "Feature Counts"
ylim = NULL

# from Pich et al.
minorInPositions = list(3:7, 13:17, 24:28, 34:38, 44:48, 55:59, 65:69)
minorOutPositions = list(-2:2, 8:12, 19:23, 29:33, 39:43, 49:54, 60:64)

# Directly from Cui and Zhurkin 2010
minorInPositions = list(5:8-74, 15:18-74, 26:29-74, 37:40-74, 47:50-74, 57:60-74, 67:70-74)
minorOutPositions = list(10:14-74, 21:23-74, 32:34-74, 42:45-74 ,52:55-74, 62:65-74)

# From Cui and Zhurkin, 2010 with 1 added to each side and half-base positions included.
minorInPositions = list(4:9-74, 14:19-74, 25:30-74, 36:41-74, 46:51-74, 56:61-74, 66:71-74,
                        5:9-74.5, 15:19-74.5, 26:30-74.5, 37:41-74.5, 47:51-74.5, 57:61-74.5, 67:71-74.5)
minorOutPositions = list(9:14-74, 20:24-74, 31:35-74, 41:46-74 ,51:56-74, 61:66-74,
                         10:14-74.5, 21:24-74.5, 32:35-74.5, 42:46-74.5, 52:56-74.5, 62:66-74.5)

# Set up plot margins.
par(mar = c(5,5,4,1))


# Average across 11 base pairs centered on the given position.
smoothValues = function(middlePos, data, dataCol, averagingRadius = 5) {
 
  positionsToAverage = (middlePos-averagingRadius):(middlePos+averagingRadius)
  valuesToAverage = data[Dyad_Position %in% positionsToAverage][[dataCol]]
  return(mean(valuesToAverage))
  
}


colorInRange = function(range, color, data, dataCol, includeNegative = TRUE) {

  lines(data[Dyad_Position %in% range, Dyad_Position], 
       data[Dyad_Position %in% range][[dataCol]], type = 'l', lwd = 3, col = color)
  
  if (includeNegative) {
    
    lines(data[Dyad_Position %in% -range, Dyad_Position], 
         data[Dyad_Position %in% -range][[dataCol]], type = 'l', lwd = 3, col = color)
    
  }
  
}


plottingSuite = function(dataSetName) {
  
  # Get the relevant counts and periodicity data for the given data set name.
  if (dataSetName %in% names(mutperiodData$normalizedNucleosomeCountsTables)) {
    countsData = mutperiodData$normalizedNucleosomeCountsTables[[dataSetName]]
  } else if (dataSetName %in% names(mutperiodData$rawNucleosomeCountsTables)) {
    countsData = mutperiodData$rawNucleosomeCountsTables[[dataSetName]]
  } else stop("Unknown data set name.")
  
  periodicityData = as.list(mutperiodData$periodicityResults[Data_Set == dataSetName])
  
  # Determine whether the data is rotational, rotational+linker, or translational.
  rotational = FALSE
  rotationalPlus = FALSE
  translational = FALSE
  if (min(countsData$Dyad_Position) >= -72) { 
    rotational = TRUE
  } else if (min(countsData$Dyad_Position) > -999) {
    rotational = TRUE
    rotationalPlus = TRUE
  } else translational = TRUE
  
  # If only rotational, trim to the cutoff value
  if (rotational && !rotationalPlus) {
    countsData = countsData[Dyad_Position >= -rotationalOnlyCutoff & Dyad_Position <= rotationalOnlyCutoff]
  }
  
  # Smooth if translational
  if (translational) {
    countsData = copy(countsData)
    countsData[, (dataCol) := sapply(countsData$Dyad_Position, smoothValues, data = countsData, dataCol = dataCol)]
  }
  
  # Basic Plotting Template (ylim = c(min,max) to set axis bounds)
  plot(countsData$Dyad_Position, countsData[[dataCol]], type = 'l', main = title,
       ylab = yAxisLabel, xlab = "Position Relative to Dyad (bp)",
       cex.lab = 2, cex.main = 1.75, cex.axis = 1.5, lwd = 3, col = "black", ylim = ylim)
  
  # Color code
  
  if (rotational) {
    # Color rotational positioning
    captureOutput = sapply(minorInPositions, colorInRange, color = "#1bcc44", data = countsData, dataCol = dataCol)
    captureOutput = sapply(minorOutPositions, colorInRange, color = "#993299", data = countsData, dataCol = dataCol)
  }
  
  if (rotationalPlus) {
    # Color linker DNA in linker+ plots.
    captureOutput = sapply(list(min(countsData$Dyad_Position):-73, 73:max(countsData$Dyad_Position)), colorInRange,
                           color = "Gold", data = countsData, dataCol = dataCol)
  }
  
  if (translational) {
    
    # Derive linker and nucleosome positions from the expected period of the data.
    nucRepLen = round(periodicityData$Expected_Peak_Periodicity)
    
    nucleosomePositions = lapply(1:10, function(x) return(append((-73+x*nucRepLen):(73+x*nucRepLen), 
                                                                 (-72.5+x*nucRepLen):(72.5+x*nucRepLen))))
    nucleosomePositions = append(nucleosomePositions, list(0:73, 0.5:72.5))
    linkerPositions = lapply(0:8, function(x) return( append((73+x*nucRepLen):(-73+(x+1)*nucRepLen), 
                                                             (72.5+x*nucRepLen):(-72.5+(x+1)*nucRepLen))))
    
    # Color translational positioning
    captureOutput = sapply(nucleosomePositions, colorInRange, color = "#0571b0", data = countsData, dataCol = dataCol)
    captureOutput = sapply(linkerPositions, colorInRange, color = "#ca0020", data = countsData, dataCol = dataCol)
  }
  
}


# Plot a bunch of figures together using facets, stratified by timepoint on one axis and domains on the other.
expectedTimepoints = c("10m", "30m", "8h", "16h", "24h")
expectedDomains = c("BLACK", "BLUE", "GREEN", "RED", "YELLOW")

addTimepointAndDomainInfo = function(dataSetName) {
  
  timepoint = names(which(sapply(expectedTimepoints, function(x) grepl(x,dataSetName))))
  domain = names(which(sapply(expectedDomains, function(x) grepl(x,dataSetName))))
  
  if (length(timepoint) == 0 || length(domain) == 0) return(data.table())
  if (length(timepoint) > 1) stop(paste("Multiple timepoints found in "),dataSetName)
  if (length(domain) > 1) stop(paste("Multiple domains found in "),dataSetName)
  
  if (dataSetName %in% names(mutperiodData$normalizedNucleosomeCountsTables)) {
    countsData = mutperiodData$normalizedNucleosomeCountsTables[[dataSetName]]
  } else if (dataSetName %in% names(mutperiodData$rawNucleosomeCountsTables)) {
    countsData = mutperiodData$rawNucleosomeCountsTables[[dataSetName]]
  } else stop("Unknown data set name.")
  
  if (grepl("nuc-group", dataSetName, fixed = TRUE)) {
    countsData = copy(countsData)
    countsData[, (dataCol) := sapply(countsData$Dyad_Position, smoothValues, data = countsData, dataCol = dataCol)]
  }
  
  return(countsData[, c("Timepoint", "Domain") := list(rep(timepoint, .N), rep(domain, .N))])
  
}

isTranslationalDataSets = grepl("nuc-group", dataSetNames, fixed = TRUE)
isRep1DataSet = grepl("rep1", dataSetNames)

# Gets all data for translational rep1 data sets with timepoint and domain columns.
stratifiedCountsData = rbindlist(lapply(dataSetNames[isTranslationalDataSets & isRep1DataSet], 
                                                     addTimepointAndDomainInfo))

# Gets all data for rotational rep1 data sets with timepoint and domain columns.
stratifiedCountsData = rbindlist(lapply(dataSetNames[!isTranslationalDataSets & isRep1DataSet], 
                                                  addTimepointAndDomainInfo))

ggplot(stratifiedCountsData,
       aes_string("Dyad_Position", dataCol, color = "Domain")) +
  scale_color_manual(values = c("BLACK" = "black", "BLUE" = "blue", "GREEN" = "forestgreen",
                                "RED" = "red", "YELLOW" = "gold2"), guide = FALSE) +
  geom_line() +
  labs(title = title, x = "Position Relative to Dyad (bp)", y = yAxisLabel) +
  facet_grid(factor(Timepoint, levels = expectedTimepoints)~Domain) +
  coord_cartesian(ylim = ylim) +
  scale_x_continuous(breaks = c(round(min(stratifiedCountsData$Dyad_Position)/2),0,
                                round(max(stratifiedCountsData$Dyad_Position)/2))) +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.title = element_text(size = 15), axis.text = element_text(size = 12),
        strip.text = element_text(size = 15))

# Plot Plus and minus strands on the same graph (normalized).
plot(data$Dyad_Position, data$Normalized_Plus_Strand, type = 'l', 
     main = "Plus (blue) vs. Minus (green)",
     ylab = "Normalized Repair Read Counts", xlab = "Position Relative to Dyad (bp)",
     cex.lab = 2, cex.main = 1.75, lwd = 3, col = "blue")
lines(data$Dyad_Position, data$Normalized_Minus_Strand, type = 'l',
      lwd = 3, col = "light green")

# Plot Plus and minus strands on the same graph (raw).
plot(data$Dyad_Position, data$Plus_Strand_Counts, type = 'l', 
     main = "Plus (blue) vs. Minus (green)",
     ylab = "Raw Repair Read Counts", xlab = "Position Relative to Dyad (bp)",
     cex.lab = 2, cex.main = 1.75, lwd = 3, col = "blue")
lines(data$Dyad_Position, data$Minus_Strand_Counts, type = 'l',
      lwd = 3, col = "light green")


# Plot nucleosome self-counts...
data[Dyad_Position == 0, Both_Strands_Counts := 0]
data[,Counts_Per_Thousand_Total := Both_Strands_Counts/sum(data$Both_Strands_Counts) * 1000]

# ... With base R plotting
plot(data$Dyad_Position, data$Counts_Per_Thousand_Total, type = 'l', main = title,
     ylab = yAxisLabel, xlab = "Position Relative to Dyad (bp)",
     cex.lab = 2, cex.main = 1.75, cex.axis = 1.5, lwd = 1.5, col = "black", ylim = ylim)
lines(data$Dyad_Position, data$Counts_Per_Thousand_Total, lwd = 1.5, col = "gold2")

#... With ggplot
ggplot(data, aes_string("Dyad_Position", dataCol, color = "Color")) +
  scale_color_manual(values = c("BLACK" = "black", "BLUE" = "blue", "GREEN" = "forestgreen",
                                "RED" = "red", "YELLOW" = "gold2"), guide = FALSE) +
  geom_line() +
  labs(title = title, x = "Position Relative to Dyad (bp)", y = yAxisLabel) +
  coord_cartesian(ylim = ylim) +
  theme(plot.title = element_text(size = 20, hjust = 0.5),
        axis.title = element_text(size = 15), axis.text = element_text(size = 12),
        strip.text = element_text(size = 15))

### Create grouped comparison figure

group1DataSetNames = sapply(strsplit(basename(mutperiodData$group1Inputs),"_nucleosome"), function(x) x[1])
group2DataSetNames = sapply(strsplit(basename(mutperiodData$group2Inputs),"_nucleosome"), function(x) x[1])

group1SNR = mutperiodData$periodicityResults[Data_Set %in% group1DataSetNames, SNR]
group2SNR = mutperiodData$periodicityResults[Data_Set %in% group2DataSetNames, SNR]

groupedSNRs = data.table(group = c(rep("MSS", length(group1SNR)),rep("MSI", length(group2SNR))),
                         SNR = c(group1SNR, group2SNR))

# ggplot dotplot (Not used or updated really...)
ggplot(groupedSNRs, aes(group, SNR)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 10, fill = NA, stroke = 2) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 2, colour = "red") +
  labs(title = "Distribution of SNR Values for Rotational Periodicities",
       x = "Microsatellite Stability") +
  scale_y_continuous(trans = "log10", breaks = trans_breaks("log10", function(x) 10^x)) +
  theme(title = element_text(size = 20))

# ggplot jittered scatter plot (Log scale)
ggplot(groupedSNRs, aes(group, SNR)) + 
  geom_jitter(width = 0.2, height = 0, shape = 1, size = 2) + 
  stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 2, colour = "red") +
  labs(title = "Translational Periodicity SNR Values",
       x = "Microsatellite Stability", y = "SNR") +
  scale_y_continuous(trans = "log10", breaks = c(1,10,100,1000)) + annotation_logticks(sides = 'l') +
  theme(plot.title = element_text(size = 20, hjust = 0.5), axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 15), axis.title.x = element_blank())

# ggplot jittered scatter plot (Linear scale)
ggplot(groupedSNRs, aes(group, SNR)) + 
  geom_jitter(width = 0.2, shape = 1, size = 2) + 
  stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 2, colour = "red") +
  labs(title = "Peak Translational Periodicity Values",
       x = "Microsatellite Stability", y = "Peak Periodicity") +
  coord_cartesian(ylim = c(1,1000)) +
  theme(plot.title = element_text(size = 20, hjust = 0.5), axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 15), axis.title.x = element_blank())

# base r boxplot
boxplot(group1SNR, group2SNR, main = "SNR distribution for rotational mutation data",
        names = c("MSS", "MSI"), ylab = "Signal to Noise Ratio (SNR)",
        cex.lab = 2, cex.axis = 1.5, cex.main = 2)
