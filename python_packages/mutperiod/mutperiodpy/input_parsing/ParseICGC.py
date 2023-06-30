# This script reads one or more "simple somatic mutation" data file(s) from ICGC and 
# writes information on single base substitution mutations to a new bed file or files for further analysis.

import os, gzip, sys, subprocess
from typing import IO, List

from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog, Selections
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (DataTypeStr, generateFilePath, getDataDirectory, checkDirs,
                                                                  generateMetadata, getIsolatedParentDir, rScriptsDirectory,
                                                                  InputFormat, getAcceptableChromosomes)
from mutperiodpy.helper_scripts.CustomErrors import *
from benbiohelpers.FileSystemHandling.DirectoryHandling import getFilesInDirectory                                                                  
from mutperiodpy.input_parsing.ParseCustomBed import parseCustomBed
from benbiohelpers.CustomErrors import *

# This class represents the mutation data obtained from ICGC in a more precise form.
# It also contains functions to represent that data in the exact format required for bed files or MSIseq, if they are valid.
class ICGCMutation:

    DNABases = ('A','C','T','G') # A list of single DNA nucleotides to check against.

    def __init__(self, choppedUpLine):

        # The relevant data extracted from a given line in the ICGC file.
        self.donorID = choppedUpLine[1]
        self.chromosome = "chr" + choppedUpLine[8] # Add "chr" to the chromosome name to follow conventions.
        if self.chromosome == "chrMT": self.chromosome = "chrM" # Standardize the naming of mitochondrial chromosomes.
        self.startPos = str(int(choppedUpLine[9]) - 1) # Make the start pos 0 based.
        self.endPos = choppedUpLine[10]
        self.mutationType = choppedUpLine[13]
        self.mutatedFrom = choppedUpLine[15]
        self.mutatedTo = choppedUpLine[16]
        if choppedUpLine[11] == "1": self.strand = '+'
        elif choppedUpLine[11] == "0": self.strand = '-'
        else: raise UserInputError("Expected strand field in 12th column containing either \"1\" for plus strand "
                               "or \"0\" for minus strand.  Found \"" + str(choppedUpLine[11]) + "\" instead.")


# This class takes an ICGCFile object and, when iterated over, returns exactly once each mutation present in a donor that is
# the result of whole genome sequencing using GRCh37 as the reference genome.
# The mutation is returned as an ICGCMutation object.
class ICGCIterator:

    def __init__(self, ICGCFile: IO, genomeFilePath):
        
        self.ICGCFile = ICGCFile # The file containing the ICGC data
        self.acceptableChromosomes = getAcceptableChromosomes(genomeFilePath)
        self.finishedDonors = list() # A list of donors to make sure we don't encounter one more than once.
        self.currentDonor = '' # The donorID currently being read for mutation data.
        self.currentDonorMutations = dict() # A dictionary of mutations unique to the current donor, to avoid writing duplicate mutations.
        self.previousDonorMutations: List[ICGCMutation] = list() # A list that keeps the last donor's mutations in order to write data on individual donors all at once.

    # Make the class iteratable, returning unique data lines one at a time.
    def __iter__(self):
        return self
    def __next__(self):

        # Read in lines until an acceptable one is found.
        while(True):

            # Read in the next line in the file.
            line = self.ICGCFile.readline()
            # Split the line into its individual data components.  The relevant components (with indices) are:
            #   0: mutation ID
            #   1: donor ID
            #   8: chromosome
            #   9: mutation start position (inclusive)
            #   10: mutation end position (inclusive)
            #   13: mutation type
            #   15: mutated from
            #   16: mutated to
            #   11: strand
            #   12: reference genome version
            #   33: sequencing method
            choppedUpLine = str(line,"utf-8").strip().split('\t')
            if len(choppedUpLine) < 34: raise StopIteration # Check to see if we have reached the end of the file.


            if (choppedUpLine[12] == "GRCh37" and # Is the reference genome hg19?
                choppedUpLine[33] == "WGS" and # Was whole genome sequencing used to generate the data?
                ("chr" + choppedUpLine[8]) in self.acceptableChromosomes): # Is the chromosome acceptable?

                # Is this a new donor?
                if choppedUpLine[1] != self.currentDonor:
                    self.finishedDonors.append(self.currentDonor)
                    self.currentDonor = choppedUpLine[1]
                    self.previousDonorMutations = list(self.currentDonorMutations.values())
                    self.currentDonorMutations.clear()
                    if self.currentDonor in self.finishedDonors:
                        raise UserInputError("Donor " + self.currentDonor + " is present in more than one block of data!")
                    print("Reading and writing from donor",self.currentDonor)

                # Have we seen this mutation in this donor?
                if choppedUpLine[0] not in self.currentDonorMutations:
                    newMutation = ICGCMutation(choppedUpLine)
                    self.currentDonorMutations[choppedUpLine[0]] = newMutation
                    return newMutation


# Handles the basic parsing of the script.
def parseICGC(ICGCFilePaths: List[str], genomeFilePath, separateDonors, 
              stratifyByMS, stratifyByMutSig):

    # Make sure suggested dependencies are installed as necessary.
    if stratifyByMS:
        failedToLoadMSISeq = False
        print("Verifying MSIseq installation...")
        try: subprocess.run(("Rscript",os.path.join(rScriptsDirectory,"TestMSIseq.R")), check = True)
        except subprocess.CalledProcessError: failedToLoadMSISeq = True
        if failedToLoadMSISeq: raise MissingMSISeqError

    if stratifyByMutSig:
        failedToLoadDeconstructSigs = False
        print("Verifying deconstructSigs installation...")
        try: subprocess.run(("Rscript",os.path.join(rScriptsDirectory,"TestDeconstructSigs.R")), check = True)
        except subprocess.CalledProcessError: failedToLoadDeconstructSigs = True
        if failedToLoadDeconstructSigs: raise MissingDeconstructSigsError

    outputBedFilePaths = list()

    if len(ICGCFilePaths) == 0: raise UserInputError("No ICGC files were found to parse.")

    # Run the parser for each ICGC file given.
    for ICGCFilePath in ICGCFilePaths:

        print("\nWorking in:",os.path.split(ICGCFilePath)[1])

        if not ICGCFilePath.endswith(".gz"):
            raise InvalidPathError(ICGCFilePath, "Given ICGC file is not gzipped (.gz file format):")
        if not "simple_somatic_mutation" in os.path.basename(ICGCFilePath):
            raise InvalidPathError(ICGCFilePath, "Given ICGC file path does not have \"simple_somatic_mutation\" in the name:",
                                   "Note: if a directory was specified to search for ICGC input files, "
                                   "all files ending in .tsv.gz are selected.")

        # Get some important file system paths for the rest of the function and generate metadata.
        dataDirectory = os.path.dirname(ICGCFilePath)
        intermediateFilesDir = os.path.join(dataDirectory,"intermediate_files")
        checkDirs(intermediateFilesDir)

        generateMetadata(getIsolatedParentDir(ICGCFilePath), getIsolatedParentDir(genomeFilePath), 
                         os.path.basename(ICGCFilePath), InputFormat.ICGC, os.path.dirname(ICGCFilePath))

        # Generate the output file.
        outputBedFilePath = generateFilePath(directory = intermediateFilesDir, dataGroup = getIsolatedParentDir(ICGCFilePath),
                                             dataType = DataTypeStr.customInput, fileExtension = ".bed")

        # Write the relevant information from the ICGC file to the output file.
        print("Writing data to custom bed format.")
        with gzip.open(ICGCFilePath, 'r') as ICGCFile:
            with open(outputBedFilePath, 'w') as outputBedFile:
                for mutation in ICGCIterator(ICGCFile, genomeFilePath):

                    # Change the formatting if a deletion or insertion is given.              
                    if mutation.mutatedFrom == '-': 
                        mutation.mutatedFrom = '*'
                        # NOTE: We are making the assumption that the given base pos (1-based) is after the insertion, not before.
                        mutation.startPos = str(int(mutation.startPos) - 1) 

                    elif mutation.mutatedTo == '-': 
                        mutation.mutatedTo = '*'

                    outputBedFile.write('\t'.join((mutation.chromosome, mutation.startPos, mutation.endPos, mutation.mutatedFrom,
                                                   mutation.mutatedTo, mutation.strand, mutation.donorID)) + '\n')

        outputBedFilePaths.append(outputBedFilePath)  

    # Pass the parsed bed files to the custom bed parser for even more parsing! (Hooray for modularization!)
    print("\nPassing data to custom bed parser...")
    parseCustomBed(outputBedFilePaths, genomeFilePath, stratifyByMS, stratifyByMutSig, separateDonors, True)


def parseArgs(args):
    
    # If only the subcommand was given, run the UI.
    if len(sys.argv) == 2: 
        main(); return

    
    # Otherwise, check to make sure valid arguments were passed:
    if args.genome_file is None: raise UserInputError("No genome file was given.")
    checkIfPathExists(args.genome_file)
    genomeFilePath = os.path.abspath(args.genome_file)

    # Get the ICGC files from the given paths, searching directories if necessary.
    finalICGCPaths = list()
    for ICGCFilePath in args.ICGCFilePaths:
        if os.path.isdir(ICGCFilePath):
            finalICGCPaths += [os.path.abspath(filePath) for filePath in getFilesInDirectory(ICGCFilePath, ".tsv.gz")]
        elif checkIfPathExists(ICGCFilePath): finalICGCPaths.append(os.path.abspath(ICGCFilePath))

    # Run the parser.
    parseICGC(list(set(finalICGCPaths)), genomeFilePath, args.stratify_by_donors, 
              args.stratify_by_microsatellite, args.stratify_by_mut_sigs)


def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory(), title = "Parse ICGC Data")
    dialog.createMultipleFileSelector("ICGC Mutation Files:",0,".tsv.gz",("gzip files",".gz"))
    dialog.createFileSelector("Genome Fasta File:",1,("Fasta Files",".fa"))
    dialog.createCheckbox("Create individual bed files for each donor.",2, 0)
    dialog.createCheckbox("Stratify results by microsatellite stability", 3, 0)
    dialog.createCheckbox("Stratify results by mutation signature", 4, 0)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    ICGCFilePaths = list(selections.getFilePathGroups())[0] # A list of ICGC mutation file paths
    genomeFilePath = list(selections.getIndividualFilePaths())[0]
    separateDonors = list(selections.getToggleStates())[0]
    stratifyByMS  = list(selections.getToggleStates())[1]
    stratifyByMutSig = list(selections.getToggleStates())[2]

    parseICGC(ICGCFilePaths, genomeFilePath, separateDonors, 
              stratifyByMS, stratifyByMutSig)

if __name__ == "__main__": main()