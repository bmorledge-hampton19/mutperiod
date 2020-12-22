# Don't forget to import the data.table library!
library(data.table)
library(ggplot2)

# Read in data.
load(choose.files(multi = FALSE))

# Column names...
# ...for normalized
dataCol = "Normalized_Both_Strands"
dataCol = "Normalized_Aligned_Strands"

# ...for raw
dataCol = "Both_Strands_Counts"

# Trim to the desired region
nucPosCutoff = 60
data = data[Dyad_Position >= -nucPosCutoff & Dyad_Position <= nucPosCutoff]

# from Pich et al.
minorInPositions = list(3:7, 13:17, 24:28, 34:38, 44:48, 55:59, 65:69)
minorOutPositions = list(-2:2, 8:12, 19:23, 29:33, 39:43, 49:54, 60:64)

# Directly from Cui and Zhurkin
minorInPositions = list(5:8-74, 15:18-74, 26:29-74, 37:40-74, 47:50-74, 57:60-74, 67:70-74)
minorOutPositions = list(10:14-74, 21:23-74, 32:34-74, 42:45-74 ,52:55-74, 62:65-74)

# From Cui and Zhurkin, with 1 added to each side.
minorInPositions = list(4:9-74, 14:19-74, 25:30-74, 36:41-74, 46:51-74, 56:61-74, 66:71-74)
minorOutPositions = list(9:14-74, 20:24-74, 31:35-74, 41:46-74 ,51:56-74, 61:66-74)

# Nucleosome vs. Linker Positions:
nucleosomePositions = append(lapply(1:10, function(x) return( (-73+x*192):(73+x*192) )),list(0:73))
linkerPositions = lapply(0:8, function(x) return( (73+x*192):(119+x*192) ))

# Set up plot margins.
par(mar = c(5,5,4,1))

# Remove outliers (roughly)
data = data[Normalized_Both_Strands < 2]

# Average across 11 base pairs centered on the given position.
smoothValues = function(middlePos, dataCol, averagingRadius = 5) {
 
  positionsToAverage = (middlePos-averagingRadius):(middlePos+averagingRadius)
  valuesToAverage = data[Dyad_Position %in% positionsToAverage][[dataCol]]
  return(mean(valuesToAverage))
  
}
data[, (dataCol) := sapply(data$Dyad_Position, smoothValues, dataCol)]

# Basic Plotting Template (ylim = c(min,max) to set axis bounds)
plot(data$Dyad_Position, data[[dataCol]], type = 'l', main = "MSS Subset Translational Periodicity",
     ylab = "Normalized Mutation Counts", xlab = "Position Relative to Dyad (bp)",
     cex.lab = 2, cex.main = 1.75, lwd = 3, col = "black", ylim = c(0.9,1.1))

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

# Plotting for paper
plot(data$Dyad_Position, data[[dataCol]], type = 'l', main = NULL,
     ylab = "Normalized Mutation Counts", xlab = "Position Relative to Dyad (bp)",
     cex.lab = 1.75, cex.axis = 1.75, lwd = 3, col = "black", ylim = c(0.9,1.1))

plot(data$Dyad_Position, data[[dataCol]], type = 'l', main = NULL,
     yaxt = 'n', ylab = '', xlab = "Position Relative to Dyad (bp)",
     cex.lab = 1.75, cex.axis = 1.75, lwd = 3, col = "black", ylim = c(0.9,1.1))

colorInRange = function(range, color, dataCol, includeNegative = TRUE) {

  lines(data[Dyad_Position %in% range, Dyad_Position], 
       data[Dyad_Position %in% range][[dataCol]], type = 'l', lwd = 3, col = color)
  
  if (includeNegative) {
    
    lines(data[Dyad_Position %in% -range, Dyad_Position], 
         data[Dyad_Position %in% -range][[dataCol]], type = 'l', lwd = 3, col = color)
    
  }
  
}

# Color rotational positioning
captureOutput = sapply(minorInPositions, colorInRange, color = "#1bcc44", dataCol = dataCol)
captureOutput = sapply(minorOutPositions, colorInRange, color = "#993299", dataCol = dataCol)

# Color translational positioning
captureOutput = sapply(linkerPositions, colorInRange, color = "#ca0020", dataCol = dataCol)
captureOutput = sapply(nucleosomePositions, colorInRange, color = "#0571b0", dataCol = dataCol)


# Create grouped comparison figure

group1DataSetNames = sapply(strsplit(basename(nucPeriodData$group1Inputs),"_nucleosome"), function(x) x[1])
group2DataSetNames = sapply(strsplit(basename(nucPeriodData$group2Inputs),"_nucleosome"), function(x) x[1])

group1SNR = nucPeriodData$periodicityResults[Data_Set %in% group1DataSetNames, SNR]
group2SNR = nucPeriodData$periodicityResults[Data_Set %in% group2DataSetNames, SNR]

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

# ggplot jittered scatter plot
ggplot(groupedSNRs, aes(group, SNR)) + 
  geom_jitter(width = 0.2, shape = 1, size = 2) + 
  stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 2, colour = "red") +
  labs(title = "Translational Periodicity SNR Values",
       x = "Microsatellite Stability", y = "SNR") +
  scale_y_continuous(trans = "log10", breaks = c(10,100,1000)) + annotation_logticks(sides = 'l') +
  theme(plot.title = element_text(size = 20, hjust = 0.5), axis.title = element_text(size = 15),
        axis.text.x = element_text(size = 15), axis.title.x = element_blank())

# base r boxplot
boxplot(group1SNR, group2SNR, main = "SNR distribution for rotational mutation data",
        names = c("MSS", "MSI"), ylab = "Signal to Noise Ratio (SNR)",
        cex.lab = 2, cex.axis = 1.5, cex.main = 2)
