library(data.table)

originalMutSigsDataPath = choose.files(multi = FALSE, caption = "Select MutSig Input Data", index = 1,
                                       filters = c(c("Tab Separated Files (*.tsv)","Any files"),c("*.tsv","*.*")))
originalMutSigsData = fread(originalMutSigsDataPath)

mSICohorts = fread(choose.files(multi = FALSE, caption = "Select file with MSI cohorts", index = 1,
                                filters = c(c("Text File (*.txt)","Any files"),c("*.txt","*.*"))), header = FALSE)[[1]]
mSIOrMSS = function(donorID) {
  if (donorID %in% mSICohorts) return("MSI")
  else return("MSS")
}

originalMutSigsData[,V1 := sapply(V1, mSIOrMSS)]

fwrite(originalMutSigsData, file.path(dirname(originalMutSigsDataPath),"MS_deconstructSigs_Data.tsv"), 
       sep = '\t', col.names = FALSE)
