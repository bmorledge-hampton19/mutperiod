par(mar = c(5,5,4,1))
MSIData = MSIData[Both_Strands_Counts > 2500]
plot(MSSData$Dyad_Position, MSSData$Normalized_Aligned_Strands, type = 'b', main = "Singlenuc Normalized Aligned tXR Seq Counts",
     ylab = "Normalized Reads Counts", xlab = "Dyad Position", cex.lab = 1.75, cex.main = 1.75)
