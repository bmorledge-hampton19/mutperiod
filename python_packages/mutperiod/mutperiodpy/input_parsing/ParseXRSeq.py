# This script takes the repair map data from tXR-Seq, calls BPDE lesions, and writes them to a bed file.
# file of estimated lesion locations.

from typing import List
import os, subprocess
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog, Selections
from mutperiodpy.helper_scripts.UsefulBioinformaticsFunctions import bedToFasta, FastaFileIterator
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (getIsolatedParentDir, generateFilePath, getDataDirectory,
                                                                  DataTypeStr, generateMetadata, InputFormat, getAcceptableChromosomes)


# Estimates (Most likely with perfect accuracy) the minimum adjusted counts value that is then assumed to represent one count.
def getMinAdjustedCountValue(tXRSeqBedGraphReadsFilePath):

    # Keep track of how many rows have been read since the last minimum was set.
    rowsSinceLastMin = 0
    minAdjustedCountValue = None

    with open(tXRSeqBedGraphReadsFilePath, 'r') as tXRSeqReadsFile:
        for line in tXRSeqReadsFile:

            # Split each line on tab characters.
            choppedUpLine = line.strip().split("\t")

            # Check for a new min.
            if minAdjustedCountValue is None or float(choppedUpLine[3]) < minAdjustedCountValue:
                rowsSinceLastMin = 0
                minAdjustedCountValue = float(choppedUpLine[3])
            else: rowsSinceLastMin += 1

            # If we've gone 1,000 lines without finding a new min, we've probably got the true min.
            # (Experimentally, it looks like the reads with the minimum value encompass over half of the reads, so...)
            if rowsSinceLastMin >= 1000:
                break

    return minAdjustedCountValue


# Create a new reads file which filters out any reads with a length greater than 28 or less than 26.
# Combine + and - reads into one file with a +/- column
# Also, convert normalized counts to actual counts and swap out the counts column for duplicated entries with multiple counts.
# The input format should be a bedGraph file.
def trimBedGraphTXRSeqData(tXRSeqBedGraphReadsFilePathPair: List[str], trimmedReadsFilePath, minAdjustedCountValue, 
                           acceptableLengths, acceptableChromosomes, roundingErrorLeniency = 0.01):
    
   with open(trimmedReadsFilePath, 'w') as trimmedReadsFile:

       # Trim/Modify all the entries from the file path pair.
       for tXRSeqBedGraphReadsFilePath in tXRSeqBedGraphReadsFilePathPair:

           plusOrMinus = tXRSeqBedGraphReadsFilePath.rsplit('.',1)[0][-1] # '+' or '-'

           with open(tXRSeqBedGraphReadsFilePath, 'r') as tXRSeqReadsFile:
               for line in tXRSeqReadsFile:

                    choppedUpLine = line.strip().split("\t")

                    # Make sure the lesion is in a valid chromosome.  Otherwise, skip it.
                    if not choppedUpLine[0] in acceptableChromosomes: continue

                    # Make sure we have a valid read length
                    readLength = int(choppedUpLine[2]) - int(choppedUpLine[1])
                    if readLength in acceptableLengths:

                        # determine the actual counts for the current entry to enter it that many times in the output file.
                        counts = float(choppedUpLine[3]) / minAdjustedCountValue

                        # Make sure the conversion to actual counts produced a reasonably whole number.
                        if abs(round(counts) - counts) > roundingErrorLeniency:
                            raise ValueError(choppedUpLine[3] + " cannot be converted to actual counts as it is not evenly divisible by" +
                                             str(minAdjustedCountValue) + ". (Rounding error > " + str(roundingErrorLeniency) + ")")
                        
                        # Write the data
                        # readID = choppedUpLine[0] + ':' + choppedUpLine[1] + '-' + choppedUpLine[2] + '(' + plusOrMinus + ')'
                        for i in range(round(counts)):
                            trimmedReadsFile.write('\t'.join((choppedUpLine[0],choppedUpLine[1],choppedUpLine[2],'NA','NA',plusOrMinus)) + '\n')


# Create a new reads file which filters out any reads with a length greater than 28 or less than 26.
# Also, replace the columns with the read ID and score with NA.
# The input format should be a bed file.
def trimBedTXRSeqData(tXRSeqBedReadsFilePath: str, trimmedReadsFilePath, acceptableLengths, acceptableChromosomes):
    
   with open(trimmedReadsFilePath, 'w') as trimmedReadsFile:

       # Trim/Modify all the entries from the file path.
        with open(tXRSeqBedReadsFilePath, 'r') as tXRSeqReadsFile:
            for line in tXRSeqReadsFile:

                choppedUpLine = line.strip().split("\t")

                # Make sure the lesion is in a valid chromosome.  Otherwise, skip it.
                if not choppedUpLine[0] in acceptableChromosomes: continue

                # Make sure we have a valid read length
                readLength = int(choppedUpLine[2]) - int(choppedUpLine[1])
                if readLength in acceptableLengths:
                    
                    # Write the data
                    trimmedReadsFile.write('\t'.join((choppedUpLine[0],choppedUpLine[1],choppedUpLine[2],'NA','NA',choppedUpLine[5])) + '\n')


# Given a nucleotide sequence, do the conditions suggest a lesion is present?
# The search location indicates the position to look for the lesion relative to one end of the sequence.
#    A value of 1 indicates the first base at the 5' end, while a value of -1 indicates the first base at the 3' end.
# The function outputs the index of the lesion in the input sequence if it is found or "None" if it is not.
def searchForLesion(sequence, expectedLocation, acceptableBases):

    assert expectedLocation != 0
    assert abs(expectedLocation) <= len(sequence), ("Sequence of length {} does not have a {} position.".
                                                    format(len(sequence), expectedLocation))

    if expectedLocation > 0: expectedLocation -= 1

    if sequence[expectedLocation] not in acceptableBases: return None

    else:
        if expectedLocation >= 0: return expectedLocation
        else: return len(sequence) + expectedLocation


# From a given fasta file of reads, find likely BPDE-dG lesions and write their
# location to a bed file.
def writeLesions(fastaReadsFilePath, lesionsBedFilePath, expectedLocationsByLength, acceptableBasesByLength):

    with open(fastaReadsFilePath, 'r') as fastaReadsFile:
        with open(lesionsBedFilePath, 'w') as lesionsBedFile:

            # Look for lesions in each fasta entry and write them to the bed file.
            for fastaEntry in FastaFileIterator(fastaReadsFile):
                for expectedLocation in expectedLocationsByLength[len(fastaEntry.sequence)]:

                    lesionLocation = searchForLesion(fastaEntry.sequence, expectedLocation, acceptableBasesByLength[len(fastaEntry.sequence)])

                    if lesionLocation is not None:
                        
                        # IMPORTANT: If the sequence is on the minus strand, the location needs to be inverted with respect to the fragment
                        # because reverse complement and stuff.
                        if fastaEntry.strand == '-': lesionLocation = len(fastaEntry.sequence) - (lesionLocation + 1)

                        bedEntry = '\t'.join((fastaEntry.sequenceLocation[0],
                                                str(int(fastaEntry.sequenceLocation[1]) + lesionLocation),
                                                str(int(fastaEntry.sequenceLocation[1]) + lesionLocation + 1),
                                                fastaEntry.sequence[lesionLocation],
                                                "OTHER",
                                                fastaEntry.sequenceLocation[3])) + '\n'
                        lesionsBedFile.write(bedEntry)


# Given a file path, looks for the complementary file path (with the opposite strand designation) and returns both as a tuple.
def getFilePathPair(filePath: str) -> str:
    
    strandReverser = {'+': '-', '-': '+'}

    fileExtension = filePath.rsplit('.', 1)[-1]
    strandDesignation = filePath.rsplit('.', 1)[-2][-1]

    if not strandDesignation in strandReverser:
        raise ValueError("No strand designation (\'+\' or \'-\') immediately preceding file extension for file path: " + filePath)

    complementaryPath = (filePath.rsplit( '.'.join((strandDesignation,fileExtension)) , 1)[0] + 
                         '.'.join((strandReverser[strandDesignation],fileExtension)) )

    if not os.path.exists(complementaryPath):
            raise ValueError("Complementary file path not found at " + complementaryPath)

    return (filePath, complementaryPath)


# A class which stores all the relevant info for tXR-Seq data as it is parsed from a variety of input formats,
# and converts them to the output format.
class TXRSeqInputDataPipeline:

    def __init__(self, inputDataFilePath: str, callParamsFilePath: str, genomeFilePath: str, nucPosFilePath: str):

        self.inputDataFilePath = inputDataFilePath
        self.genomeFilePath = genomeFilePath
        self.nucPosFilePath = nucPosFilePath

        # Get the list of acceptable chromosomes
        self.acceptableChromosomes = getAcceptableChromosomes(self.genomeFilePath)

        # Read in information from the callParams File.
        self.callParamsFilePath = callParamsFilePath
        self.expectedLocationsByLength = dict()
        self.acceptableBasesByLength = dict()

        with open(self.callParamsFilePath, 'r') as callParamsFile:
            for line in callParamsFile:

                choppedUpLine: List[str] = line.split()

                sequenceLength = int(choppedUpLine[0])
                assert sequenceLength not in self.expectedLocationsByLength

                self.expectedLocationsByLength[sequenceLength] = [int(location) for location in choppedUpLine[1].split(',')]
                self.acceptableBasesByLength[sequenceLength] = choppedUpLine[2].split(',')

        # Initialize these values as empty until we know what form the input data is in.
        self.bigWigReadsFilePathPair = None
        self.bedGraphReadsFilePathPair = None
        self.bedInputFilePath = None

        # Initialize a variable to keep track of whether or not the input reads have been trimmed.
        self.readsHaveBeenTrimmed = False

        # Determine what form the input data is in and assign it accordingly.
        if self.inputDataFilePath.endswith(".bigWig"):
            self.bigWigReadsFilePathPair = getFilePathPair(self.inputDataFilePath)

        elif self.inputDataFilePath.endswith(".bedGraph"):
            self.bedGraphReadsFilePathPair = getFilePathPair(self.inputDataFilePath)

        elif self.inputDataFilePath.endswith(".bed"):
            self.bedInputFilePath = self.inputDataFilePath

        else: raise ValueError("Input data file path: " + self.inputDataFilePath + " is not in an acceptable format.  " +
                               "Expected a bigWig, bedGraph, or bed file type.")

        self.setUpFileSystem()


    # Prepares and stores data about the file system surrounding the input data for the rest of the pipeline.
    def setUpFileSystem(self):

        # Store useful paths and names.
        localRootDirectory = os.path.dirname(self.inputDataFilePath)
        dataGroupName = getIsolatedParentDir(self.inputDataFilePath)

        # Create the intermediate files directory if necessary
        intermediateFilesDirectory = os.path.join(localRootDirectory,"intermediate_files")
        if not os.path.exists(intermediateFilesDirectory):
            os.mkdir(intermediateFilesDirectory)

        # If the given input data is in bigWig form, use it to generate the file paths to the intermediate bedGraph files.
        if self.bigWigReadsFilePathPair is not None:
            self.bedGraphReadsFilePathPair = list()
            for bigWigReadsFilePath in self.bigWigReadsFilePathPair:
                self.bedGraphReadsFilePathPair.append(os.path.join(intermediateFilesDirectory,
                                                                   os.path.basename(bigWigReadsFilePath).rsplit('.',1)[0]+".bedGraph"))

        # Generate the trimmed reads output, the fasta output, and bed output file paths.
        self.trimmedReadsFilePath = os.path.join(intermediateFilesDirectory,dataGroupName+"_trimmed_reads.bed")
        self.fastaReadsFilePath = os.path.join(intermediateFilesDirectory,dataGroupName+"_trimmed_reads.fa")
        self.lesionsBedFilePath = generateFilePath(directory = localRootDirectory, dataGroup = dataGroupName,
                                                      context = "singlenuc", dataType = DataTypeStr.mutations, fileExtension = ".bed") 

        # Generate metadata
        generateMetadata(dataGroupName, getIsolatedParentDir(self.genomeFilePath), getIsolatedParentDir(self.nucPosFilePath),
                         os.path.basename(self.inputDataFilePath), InputFormat.tXRSeq, localRootDirectory,
                         callParamsFilePath = self.callParamsFilePath)

    
    # Generates a bed file of trimmed reads, a process which is dependent on the input format given in the class's initialization.
    def generateTrimmedReads(self):
        
        # If the input data is in bed format, generate the trimmed reads from that data.
        if self.bedInputFilePath is not None:
            trimBedTXRSeqData(self.bedInputFilePath, self.trimmedReadsFilePath, 
                              self.acceptableBasesByLength.keys(), self.acceptableChromosomes)

        # Otherwise, generate the trimmed reads from bedGraph data.
        else:

            # If we were given bigWig input data, convert it to bedGraph format first.
            if self.bigWigReadsFilePathPair is not None:
                for i in range(2):
                    subprocess.run(("bigWigToBedGraph",self.bigWigReadsFilePathPair[i],
                                    self.bedGraphReadsFilePathPair[i]), check = True)

            # Generate the trimmed reads.
            trimBedGraphTXRSeqData(self.bedGraphReadsFilePathPair, self.trimmedReadsFilePath, 
                                   getMinAdjustedCountValue(self.bedGraphReadsFilePathPair[0]),
                                   self.acceptableBasesByLength.keys(), self.acceptableChromosomes)

        self.readsHaveBeenTrimmed = True

    
    # Generates the output file that contains the locations of estimated lesions.
    def generateLesionsBedOutputFile(self):

        if not self.readsHaveBeenTrimmed: raise ValueError("Trying to generate final output without trimmed reads.")

        # Convert the trimmed file to fasta format to get the sequences associated with each read.
        bedToFasta(self.trimmedReadsFilePath, self.genomeFilePath, self.fastaReadsFilePath)

        # Find the lesions and write them to the final output file.
        writeLesions(self.fastaReadsFilePath, self.lesionsBedFilePath, 
                     self.expectedLocationsByLength, self.acceptableBasesByLength)

        # Sort the output file.
        subprocess.run(("sort","-k1,1","-k2,2n", self.lesionsBedFilePath,"-o", self.lesionsBedFilePath), check = True)    


def parseTXRSeq(inputDataFilePaths, callParamsFilePath, genomeFilePath, nucPosFilePath):

    # Use the TXRSeqInputDataPipeline object to convert the input data to the output data!
    for inputDataFilePath in inputDataFilePaths:

        print("\nWorking with:",os.path.basename(inputDataFilePath))
        tXRSeqInputDataPipeline = TXRSeqInputDataPipeline(inputDataFilePath, callParamsFilePath, genomeFilePath, nucPosFilePath)

        print("Generating trimmed reads bed file...")
        tXRSeqInputDataPipeline.generateTrimmedReads()

        print("Locating and writing lesion locations to final output file...")
        tXRSeqInputDataPipeline.generateLesionsBedOutputFile()
 

if __name__ == "__main__":

    # Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("XR-seq bigwig data (plus strand):",0,
                                      "+.bigWig",("BigWig Files",".bigWig"))
    dialog.createMultipleFileSelector("XR-seq bed data (alternative to bigwig):",1,"aligned_reads.bed",
                                      ("Bed Files",".bed"), additionalFileEndings=("Rep1.bed","Rep2.bed"))

    dialog.createFileSelector("Lesion Call Parameter File:", 2, ("Tab Seperated Values",".tsv"))

    dialog.createFileSelector("Genome Fasta File:",3,("Fasta Files",".fa"))
    dialog.createFileSelector("Strongly Positioned Nucleosome File:",4,("Bed Files",".bed"))

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    tXRSeqBigWigPlusReadsFilePaths: List[str] = list(selections.getFilePathGroups())[0]
    tXRSeqBedReadsFilePaths: List[str] = list(selections.getFilePathGroups())[1]
    callParamsFilePath = selections.getIndividualFilePaths()[0]
    genomeFilePath = selections.getIndividualFilePaths()[1]
    nucPosFilePath = selections.getIndividualFilePaths()[2]

    parseTXRSeq(tXRSeqBigWigPlusReadsFilePaths + tXRSeqBedReadsFilePaths, callParamsFilePath, genomeFilePath, nucPosFilePath)