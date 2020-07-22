# Find the number of tumor samples (specimens?) per donor in a given data set.

library(data.table)

# Read in the file.
dataPath = file.choose()
data = fread(input = dataPath)

# Find the problem children!
donorInfo = data[, .(sample_IDs = list(unique(icgc_sample_id)), 
                     Sample_Count = length(unique(icgc_specimen_id))), by = icgc_donor_id]
multiSampleDonors = donorInfo[Sample_Count > 1]