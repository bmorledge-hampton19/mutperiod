##### Set Up the Environment #####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("MSIseq")

library(MSIseq)

# Test Data
data(NGStestdata)
data(NGStestseqLen)

## download the Hg19repeats annotation file and load it
url <- "http://steverozen.net/data/Hg19repeats.rda"
file <- basename(url)
download.file(url, file)
load("Hg19repeats.rda")

## get mutation counts for test data
testMutationNum<-Compute.input.variables(NGStestdata,repeats=Hg19repeats, seq.len = NGStestseqLen)

# Determine the presence of microsattelite instability
result <- MSIseq.classify(mutationNum = testMutationNum)
