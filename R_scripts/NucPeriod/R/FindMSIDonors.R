# This function uses the MSIseq package to generate a list of MSI donors.
#' @export
findMSIDonors = function(inputDataPath, verbose = TRUE) {

  # The number of base pairs (in megabases) in the hg19 genome (excluding mitochondria).
  captureLength = 3096

  ## download the Hg19repeats annotation file and load it
  if (verbose) print("Downloading genome information...")
  url <- "http://steverozen.net/data/Hg19repeats.rda" # NOTE: DON'T CHANGE THIS LINE WITHOUT CHANGING CAPTURE LENGTH
  file <- basename(url)
  if (!file.exists(file)) download.file(url, file)
  load(file)

  # Get the input data
  inputData = fread(inputDataPath)
  setnames(inputData,colnames(inputData),c("Chrom","Start_Position","End_Position",
                                           "Variant_Type","Tumor_Sample_Barcode"))

  # Remove mitochondrial mutations.
  inputData = inputData[Chrom != "chrM"]

  # Get mutation counts for input data.
  if (verbose) print("Running analysis...")
  mutationCounts = Compute.input.variables(inputData,Hg19repeats,captureLength)
  result = as.data.table(MSIseq.classify(mutationCounts))

  # Check for polE deficiency.
  PolEDeficientResults = result[Likely_POLE_deficiency == "Yes"]
  if (nrow(PolEDeficientResults) > 0) {
    print(paste(nrow(PolEDeficientResults),"PolE deficient result(s)."))
  }

  # Export the MSI-H results to a text file.
  if (verbose) print("Writing Results...")
  MSIResults = result[MSI_status == "MSI-H",Tumor_Sample_Barcode]
  write(as.character(MSIResults),sep = '\n', file = paste0(dirname(inputDataPath),'/',
                                                           unlist(strsplit(basename(inputDataPath),'_'))[1],
                                                           "_MSI_donors.txt"))
}
