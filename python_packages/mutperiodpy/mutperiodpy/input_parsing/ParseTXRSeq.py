# This script takes the repair map data from tXR-Seq and converts it to a trinucleotide context
# file of estimated lesion locations.

from typing import List
import os, subprocess
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog, Selections
from mutperiodpy.helper_scripts.UsefulBioinformaticsFunctions import bedToFasta, FastaFileIterator, baseChromosomes
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (getIsolatedParentDir, generateFilePath, dataDirectory,
                                                                  DataTypeStr, generateMetadata, InputFormat)


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
                           minReadLength = 26, maxReadLength = 28, roundingErrorLeniency = 0.01):
    
   with open(trimmedReadsFilePath, 'w') as trimmedReadsFile:

       # Trim/Modify all the entries from the file path pair.
       for tXRSeqBedGraphReadsFilePath in tXRSeqBedGraphReadsFilePathPair:

           plusOrMinus = tXRSeqBedGraphReadsFilePath.rsplit('.',1)[0][-1] # '+' or '-'

           with open(tXRSeqBedGraphReadsFilePath, 'r') as tXRSeqReadsFile:
               for line in tXRSeqReadsFile:

                    choppedUpLine = line.strip().split("\t")

                    # Make sure the lesion is in a valid chromosome.  Otherwise, skip it.
                    if not choppedUpLine[0] in baseChromosomes: continue

                    # Make sure we have a valid read length
                    readLength = int(choppedUpLine[2]) - int(choppedUpLine[1])
                    if not (readLength < minReadLength or readLength > maxReadLength):

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
def trimBedTXRSeqData(tXRSeqBedReadsFilePath: str, trimmedReadsFilePath, minReadLength = 26, maxReadLength = 28):
    
   with open(trimmedReadsFilePath, 'w') as trimmedReadsFile:

       # Trim/Modify all the entries from the file path.
        with open(tXRSeqBedReadsFilePath, 'r') as tXRSeqReadsFile:
            for line in tXRSeqReadsFile:

                choppedUpLine = line.strip().split("\t")

                # Make sure the lesion is in a valid chromosome.  Otherwise, skip it.
                if not choppedUpLine[0] in baseChromosomes: continue

                # Make sure we have a valid read length
                readLength = int(choppedUpLine[2]) - int(choppedUpLine[1])
                if not (readLength < minReadLength or readLength > maxReadLength):
                    
                    # Write the data
                    trimmedReadsFile.write('\t'.join((choppedUpLine[0],choppedUpLine[1],choppedUpLine[2],'NA','NA',choppedUpLine[5])) + '\n')


# Given a nucleotide sequence, find the most probable location of the BPDE-dG lesion.
# The lesion 3' shift value is the number of bases to shift from the 3' end to determine the expected lesion location.
#    A value of 1 indicates the first base at the 3' end.
# The exclusion values represent how many bases around the expected lesion location should be free of G's for the lesion 
#    to be considered valid.
# The function outputs the index of the lesion in the input sequence if it is found or "None" if it is not.
def findBPDEdGLesion(sequence, lesionThreePrimeShift = 7, upstreamExclusion = 0, downstreamExclusion = 0):

    if sequence[-lesionThreePrimeShift] != 'G': return None

    for i in range(-lesionThreePrimeShift - upstreamExclusion, -lesionThreePrimeShift):
        if sequence[i] == 'G': return None

    for i in range(-lesionThreePrimeShift + 1, -lesionThreePrimeShift + downstreamExclusion + 1):
        if sequence[i] == 'G': return None

    return len(sequence) - lesionThreePrimeShift


# From a given fasta file of reads, find likely BPDE-dG lesions and write their
# location and trinuc context to a bed file.
def writeTrinucLesions(fastaReadsFilePath, trinucLesionsFilePath):

    with open(fastaReadsFilePath, 'r') as fastaReadsFile:
        with open(trinucLesionsFilePath, 'w') as trinucLesionsFile:

            # Look for lesions in each fasta entry and write them to the trinuc file.
            for fastaEntry in FastaFileIterator(fastaReadsFile):

                lesionLocation = findBPDEdGLesion(fastaEntry.sequence)

                if lesionLocation is not None:
                    trinucEntry = '\t'.join((fastaEntry.sequenceLocation[0],
                                             str(int(fastaEntry.sequenceLocation[1]) + lesionLocation),
                                             str(int(fastaEntry.sequenceLocation[1]) + lesionLocation + 1),
                                             fastaEntry.sequence[lesionLocation-1:lesionLocation+2],
                                             "OTHER",
                                             fastaEntry.sequenceLocation[3])) + '\n'
                    trinucLesionsFile.write(trinucEntry)


def parseTXRSeq(tXRSeqBigWigReadsFilePathPairs, tXRSeqBedReadsFilePaths, genomeFilePath, nucPosFilePath):
    
    # Parse reads given in bigwig file pair format.
    for tXRSeqBigWigReadsFilePathPair in tXRSeqBigWigReadsFilePathPairs:

        print("\nWorking in:",os.path.basename(tXRSeqBigWigReadsFilePathPair[0]),
              "and",os.path.basename(tXRSeqBigWigReadsFilePathPair[1]))
        if (not os.path.basename(tXRSeqBigWigReadsFilePathPair[0]).endswith(".bigWig") or 
            not os.path.basename(tXRSeqBigWigReadsFilePathPair[1]).endswith(".bigWig")):
            raise ValueError("Error:  Expected bigWig file format.")

        ### Prep the file system for the output files.

        # Store useful paths and names.
        localRootDirectory = os.path.dirname(tXRSeqBigWigReadsFilePathPair[0])
        dataGroupName = getIsolatedParentDir(tXRSeqBigWigReadsFilePathPair[0])

        # Create the intermediate files directory if necessary
        intermediateFilesDirectory = os.path.join(localRootDirectory,"intermediate_files")
        if not os.path.exists(intermediateFilesDirectory):
            os.mkdir(intermediateFilesDirectory)

        # Generate the file paths to the intermediate bedGraph files.
        tXRSeqBedGraphReadsFilePathPair = list()
        for tXRSeqBigWigReadsFilePath in tXRSeqBigWigReadsFilePathPair:
            tXRSeqBedGraphReadsFilePathPair.append(os.path.join(intermediateFilesDirectory,
                                                                os.path.basename(tXRSeqBigWigReadsFilePath).rsplit('.',1)[0]+".bedGraph"))

        # Generate the trimmed reads output, the fasta output, and trinuc lesions output file paths.
        trimmedReadsFilePath = os.path.join(intermediateFilesDirectory,dataGroupName+"_trimmed_reads.bed")
        fastaReadsFilePath = os.path.join(intermediateFilesDirectory,dataGroupName+"_fasta_reads.fa")
        trinucLesionsFilePath = generateFilePath(directory = localRootDirectory, dataGroup = dataGroupName,
                                                 context = "trinuc", dataType = DataTypeStr.mutations, fileExtension = ".bed")
        
        # Generate metadata
        generateMetadata(dataGroupName, getIsolatedParentDir(genomeFilePath), getIsolatedParentDir(nucPosFilePath),
                         os.path.basename(tXRSeqBigWigReadsFilePathPair[0]), InputFormat.tXRSeqBigWig, localRootDirectory)

        # Convert from bigWig to bedGraph format.
        print("Converting to bedGraph...")
        for i in range(2):
            subprocess.run(" ".join(("bigWigToBedGraph",tXRSeqBigWigReadsFilePathPair[i],
                                     tXRSeqBedGraphReadsFilePathPair[i])), shell = True, check = True)

        # Trim the bedGraph file pair.
        print("Trimming bedGraph data and combining to one bed file...")
        trimBedGraphTXRSeqData(tXRSeqBedGraphReadsFilePathPair, trimmedReadsFilePath, 
                               getMinAdjustedCountValue(tXRSeqBedGraphReadsFilePathPair[0]))

        # Convert the trimmed file to fasta format to get the sequences associated with each read.
        print("Generating fasta file from trimmed bed file...")
        bedToFasta(trimmedReadsFilePath, genomeFilePath, fastaReadsFilePath)

        # Find the lesions and write them in their trinuc context to the final output file.
        print("Finding lesions and writing their trinuc context to final output file...")
        writeTrinucLesions(fastaReadsFilePath, trinucLesionsFilePath)

        # Sort the trinuc context file.
        print("Sorting trinuc context data...")
        subprocess.run(" ".join(("sort","-k1,1","-k2,2n",trinucLesionsFilePath,"-o",trinucLesionsFilePath)), 
                       shell = True, check = True)


    # Parse reads given in bed format.
    for tXRSeqBedReadsFilePath in tXRSeqBedReadsFilePaths:

        print("\nWorking in:",os.path.basename(tXRSeqBedReadsFilePath))
        if not os.path.basename(tXRSeqBedReadsFilePath).endswith(".bed"):
            raise ValueError("Error:  Expected bed file format.")

        ### Prep the file system for the output files.

        # Store useful paths and names.
        localRootDirectory = os.path.dirname(tXRSeqBedReadsFilePath)
        dataGroupName = getIsolatedParentDir(tXRSeqBedReadsFilePath)

        # Create the intermediate files directory if necessary
        intermediateFilesDirectory = os.path.join(localRootDirectory,"intermediate_files")
        if not os.path.exists(intermediateFilesDirectory):
            os.mkdir(intermediateFilesDirectory)

        # Generate the trimmed reads output, the fasta output, and trinuc lesions output file paths.
        trimmedReadsFilePath = os.path.join(intermediateFilesDirectory,dataGroupName+"_trimmed_reads.bed")
        fastaReadsFilePath = os.path.join(intermediateFilesDirectory,dataGroupName+"_fasta_reads.fa")
        trinucLesionsFilePath = generateFilePath(directory = localRootDirectory, dataGroup = dataGroupName,
                                                 context = "trinuc", dataType = DataTypeStr.mutations, fileExtension = ".bed")

        # Generate metadata.
        generateMetadata(dataGroupName, getIsolatedParentDir(genomeFilePath), getIsolatedParentDir(nucPosFilePath),
                         os.path.basename(tXRSeqBedReadsFilePath), InputFormat.tXRSeqBed, localRootDirectory)

        # Trim the bedGraph file pair.
        print("Trimming bed data...")
        trimBedTXRSeqData(tXRSeqBedReadsFilePath, trimmedReadsFilePath)

        # Convert the trimmed file to fasta format to get the sequences associated with each read.
        print("Generating fasta file from trimmed bed file...")
        bedToFasta(trimmedReadsFilePath, genomeFilePath, fastaReadsFilePath)

        # Find the lesions and write them in their trinuc context to the final output file.
        print("Finding lesions and writing their trinuc context to final output file...")
        writeTrinucLesions(fastaReadsFilePath, trinucLesionsFilePath)

        # Sort the trinuc context file.
        print("Sorting trinuc context data...")
        subprocess.run(" ".join(("sort","-k1,1","-k2,2n",trinucLesionsFilePath,"-o",trinucLesionsFilePath)), 
                       shell = True, check = True)

if __name__ == "__main__":

    # Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=dataDirectory)
    dialog.createMultipleFileSelector("tXR-seq bigwig data (plus strand):",0,
                                      "+.bigWig",("BigWig Files",".bigWig"))
    dialog.createMultipleFileSelector("tXR-seq bed data (alternative to bigwig):",1,
                                      "aligned_reads.bed",("Bed Files",".bed"))
    dialog.createFileSelector("Human Genome Fasta File:",2,("Fasta Files",".fa"))
    dialog.createFileSelector("Strongly Positioned Nucleosome File:",3,("Bed Files",".bed"))                    
    dialog.createExitButtons(4,0)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    tXRSeqBigWigPlusReadsFilePaths: List[str] = list(selections.getFilePathGroups())[0]
    tXRSeqBedReadsFilePaths: List[str] = list(selections.getFilePathGroups())[1]
    genomeFilePath = list(selections.getIndividualFilePaths())[0]
    nucPosFilePath = list(selections.getIndividualFilePaths())[1]

    # Pair up the plus and minus strands for each data set.
    tXRSeqBigWigReadsFilePathPairs = list()
    for tXRSeqBigWigPlusReadsFilePath in tXRSeqBigWigPlusReadsFilePaths:
        correspondingMinusPath = tXRSeqBigWigPlusReadsFilePath.rsplit("+.bigWig")[0] + "-.bigWig"
        if not os.path.exists(correspondingMinusPath):
            raise ValueError("Corresponding minus file path not found at " + correspondingMinusPath)
        else:
            tXRSeqBigWigReadsFilePathPairs.append((tXRSeqBigWigPlusReadsFilePath, correspondingMinusPath))

    parseTXRSeq(tXRSeqBigWigReadsFilePathPairs, tXRSeqBedReadsFilePaths, genomeFilePath, nucPosFilePath)