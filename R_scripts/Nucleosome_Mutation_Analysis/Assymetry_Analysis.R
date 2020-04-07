###### NOTE ######
# A lot of these fnctions are not really optimized for R as they incrementally append to vectors.
# If these functions are ever used in any large scale computing, they should probably be rewritten.

# Runs all of the three functions below at once.
runExtremeAnalysisSuite = function(data, maxes = TRUE, dyadPosCutoff = 73) {
  
  extremeDyadLocations = getDyadLocationsOfExtremes(data, maxes = maxes, dyadPosCutoff = dyadPosCutoff)
  extremes = extremeDyadLocationsToValue(data, extremeDyadLocations)
  names(extremes) = extremeDyadLocations
  analyzeExtremeAssymetry(extremes)
  
}

# Given a vector of values, returns the dyad position of either the peaks or valleys in the data.
# Locations are smoothed in the event of closely grouped peaks.
# The cutoff value excludes dyad positions with a greater absolute value than the cutoff.
getDyadLocationsOfExtremes = function(data, maxes=TRUE, dyadPosCutoff = 73) {
  
  if (length(data) != 147) stop("Error:  Expected data points for dyad positions -73 through 73")
  
  # Look for local extremes in a sliding 11 nucleotide long area.
  extremes = numeric()
  for (i in 1:137) {
    
    if (i - 1 - 73 >= -dyadPosCutoff && i - 1 - 73 + 10 <= dyadPosCutoff) {
      currentWindow = data[i:(i+10)]
      
      if (maxes) {
        currentExtreme = match(max(currentWindow),currentWindow) + i - 2 - 73
      } else { 
        currentExtreme = match(min(currentWindow),currentWindow) + i - 2 - 73
      }
        
      if (!currentExtreme %in% extremes) {
        extremes = append(extremes,currentExtreme)
      }
    }
    
  }
  
  #Now, take any groups of extremes that are at most 2 nucleotides apart and average them.
  smoothedExtremes = numeric()
  extremePos = 1
  while(extremePos < length(extremes)) {
    
    extremesToAverage = extremes[extremePos]
    
    while(extremePos < length(extremes) && extremes[extremePos] - extremes[extremePos+1] > -3) {
      
      extremePos = extremePos + 1
      extremesToAverage = append(extremesToAverage,extremes[extremePos])
      
    }
    
    smoothedExtremes = append(smoothedExtremes,mean(extremesToAverage))
    extremePos = extremePos + 1
    
    if (extremePos == length(extremes)) {
      smoothedExtremes = append(smoothedExtremes,extremes[extremePos])
    }
    
  }
  
  return(smoothedExtremes)
  
}

# This function Converts the dyad location values returned from the above function to their corresponding
# mutation values in a given data set.
extremeDyadLocationsToValue = function(data, extremeDyadLocations) {
  
  extremeValues = numeric(length(extremeDyadLocations))
  
  #If we have a location ending in ".5", average the surrounding two values.
  for (i in 1:length(extremeDyadLocations)) {
    
    extremeDataIndex = extremeDyadLocations[i] + 1 + 73
    
    if (extremeDataIndex - .5 == as.integer(extremeDataIndex)) {
      extremeValues[i] = mean(data[extremeDataIndex+.5],data[extremeDataIndex-.5])
    } else {
      extremeValues[i] = data[as.integer(extremeDataIndex)]
    }
    
  }
  
  return(extremeValues)
  
}

# Given a named vector of extreme values, split the vector into two halves (roughly centered about the dyad)
# and perform a paired t-test on the obtained values.
analyzeExtremeAssymetry = function(extremes) {
  
  # Two vectors to store each of the halves of the data.
  # The floor and celing functions work to exclude the middle extreme from either half
  # if there are an odd number of total extremes.
  firstHalf = extremes[ 1 : (floor(length(extremes)/2)) ]
  secondHalf = rev(extremes[ (ceiling(length(extremes)/2)+1) : length(extremes) ])

  # Run the paired t-test
  return(t.test(firstHalf,secondHalf))
  
}

####### Deprecated... Probably ######
# #Maybe like this?
# plusStrandMutationRatios = numeric(73)
# minusStrandMutationRatios = numeric(73)
# 
# for (j in 1:73) {
#   plusStrandMutationRatios[j] = normalizedData[Dyad_Position == j, Normalized_Plus_Strand] / 
#     normalizedData[Dyad_Position == -j, Normalized_Plus_Strand]
#   minusStrandMutationRatios[j] = normalizedData[Dyad_Position == j, CPDs_Minus_Strand] / 
#     normalizedData[Dyad_Position == -j, CPDs_Minus_Strand]
# }
# 
# print(paste0("Mean Ratio for plus strands: (+/-) ", mean(plusStrandMutationRatios)))
# print(paste0("Mean Ratio for minus strands: (+/-) ", mean(minusStrandMutationRatios)))
# 
# 
# #Or maybe like this?
# 
# plusStrandMutationRatio = sum(normalizedData[75:147,Normalized_Plus_Strand]) /
#   sum(normalizedData[1:73,Normalized_Plus_Strand])
# 
# minusStrandMutationRatio = sum(normalizedData[75:147,Normalized_Minus_Strand]) /
#   sum(normalizedData[1:73,Normalized_Minus_Strand])
# 
# print(paste0("Summed Ratio for pluus strands: (+/-) ", plusStrandMutationRatio))
# print(paste0("Summed Ratio for minus strands: (+/-) ", minusStrandMutationRatio))
# 
# 
# # Maybe even this?
# t.test(normalizedData[75:147,Normalized_Plus_Strand],paired = TRUE,normalizedData[1:73,Normalized_Plus_Strand])
# t.test(normalizedData[75:147,Normalized_Minus_Strand],paired = TRUE,normalizedData[1:73,Normalized_Minus_Strand])
# 
# 


# # Look for peaks in each normalized data set
# findPeaks(normalizedData[,Normalized_Plus_Strand]) # Not super helpful?