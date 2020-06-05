# This script takes the repair map data from tXR-Seq and converts it to a trinucleotide context
# file of estimated lesion locations.
from TkinterDialog import Selections, TkinterDialog
from UsefulBioinformaticsFunctions import bedToFasta, FastaFileIterator
from UsefulFileSystemFunctions import getIsolatedParentDir
from typing import List
import os, subprocess

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
def trimTXRSeqData(tXRSeqBedGraphReadsFilePathPair: List[str], trimmedReadsFilePath, minAdjustedCountValue, 
                   minReadLength = 26, maxReadLength = 28, roundingErrorLeniency = 0.01):
    
   with open(trimmedReadsFilePath, 'w') as trimmedReadsFile:

       # Trim/Modify all the entries from the file path pair.
       for tXRSeqBedGraphReadsFilePath in tXRSeqBedGraphReadsFilePathPair:

           plusOrMinus = tXRSeqBedGraphReadsFilePath.rsplit('.',1)[0][-1] # '+' or '-'

           with open(tXRSeqBedGraphReadsFilePath, 'r') as tXRSeqReadsFile:
               for line in tXRSeqReadsFile:

                    choppedUpLine = line.strip().split("\t")

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


# Given a nucleotide sequence, find the most probable location of the BPDE-dG lesion.
# The lesion 3' shift value is the number of bases to shift from the 3' end to determine the expected lesion location.
#    A value of 1 indicates the first base at the 3' end.
# The exclusion values represent how many bases around the expected lesion location should be free of G's for the lesion 
#    to be considered valid.
# The function outputs the index of the lesion in the input sequence if it is found or "None" if it is not.
def findBPDEdGLesion(sequence, lesionThreePrimeShift = 7, upstreamExclusion = 3, downstreamExclusion = 4):

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
                                             "NA",
                                             fastaEntry.sequenceLocation[3])) + '\n'
                    trinucLesionsFile.write(trinucEntry)


def parseTXRSeq(tXRSeqReadsFilePathPairs, humanGenomeFastaFilePath):
    
    for tXRSeqReadsFilePathPair in tXRSeqReadsFilePathPairs:

        print("\nWorking in:",os.path.basename(tXRSeqReadsFilePathPair[0]),
              "and",os.path.basename(tXRSeqReadsFilePathPair[1]))
        if (not os.path.basename(tXRSeqReadsFilePathPair[0]).endswith(".bigWig") or 
            not os.path.basename(tXRSeqReadsFilePathPair[1]).endswith(".bigWig")):
            raise ValueError("Error:  Expected bigWig file format.")

        ### Prep the file system for the output files.

        # Store useful paths and names.
        localRootDirectory = os.path.dirname(tXRSeqReadsFilePathPair[0])
        dataGroupName = getIsolatedParentDir(tXRSeqReadsFilePathPair[0])

        # Create the intermediate files directory if necessary
        intermediateFilesDirectory = os.path.join(localRootDirectory,"intermediate_files")
        if not os.path.exists(intermediateFilesDirectory):
            os.mkdir(intermediateFilesDirectory)

        # Generate the file paths to the intermediate bedGraph files.
        tXRSeqBedGraphReadsFilePathPair = list()
        for tXRSeqReadsFilePath in tXRSeqReadsFilePathPair:
            tXRSeqBedGraphReadsFilePathPair.append(os.path.join(intermediateFilesDirectory,
                                                                os.path.basename(tXRSeqReadsFilePath).rsplit('.',1)[0]+".bedGraph"))

        # Generate the trimmed reads output, the fasta output, and trinuc lesions output file paths.
        trimmedReadsFilePath = os.path.join(intermediateFilesDirectory,dataGroupName+"_trimmed_reads.bed")
        fastaReadsFilePath = os.path.join(intermediateFilesDirectory,dataGroupName+"_fasta_reads.fa")
        trinucLesionsFilePath = os.path.join(localRootDirectory,dataGroupName+"_trinuc_context.bed")

        # Convert from bigWig to bedGraph format.
        print("Converting to bedGraph...")
        for i in range(2):
            subprocess.run(" ".join(("bigWigToBedGraph",tXRSeqReadsFilePathPair[i],
                                     tXRSeqBedGraphReadsFilePathPair[i])), shell = True, check = True)

        # Trim the bedGraph file pair.
        print("Trimming bedGraph data and combining to one bed file...")
        trimTXRSeqData(tXRSeqBedGraphReadsFilePathPair, trimmedReadsFilePath, 
                       getMinAdjustedCountValue(tXRSeqBedGraphReadsFilePathPair[0]))

        # Convert the trimmed file to fasta format to get the sequences associated with each read.
        print("Generating fasta file from trimmed bed file...")
        bedToFasta(trimmedReadsFilePath, humanGenomeFastaFilePath, fastaReadsFilePath)

        # Find the lesions and write them in their trinuc context to the final output file.
        print("Finding lesions and writing their trinuc context to final output file...")
        writeTrinucLesions(fastaReadsFilePath, trinucLesionsFilePath)

        # Sort the trinuc context file.
        print("Sorting trinuc context data...")
        subprocess.run(" ".join(("sort","-k1,1","-k2,2n",trinucLesionsFilePath,"-o",trinucLesionsFilePath)), 
                        shell = True, check = True)


if __name__ == "__main__":

    # Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(__file__),"..","data"))
    dialog.createMultipleFileSelector("tXR-seq read data (plus strand):",0,
                                      "+.bigWig",("BigWig Files",".bigWig"))
    dialog.createFileSelector("Human Genome Fasta File:",1,("Fasta Files",".fa"))                    
    dialog.createReturnButton(2,0,2)
    dialog.createQuitButton(2,2,2)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    tXRSeqPlusReadsFilePaths: List[str] = list(selections.getFilePathGroups())[0]
    humanGenomeFastaFilePath = list(selections.getIndividualFilePaths())[0]

    # Pair up the plus and minus strands for each data set.
    tXRSeqReadsFilePathPairs = list()
    for tXRSeqPlusReadsFilePath in tXRSeqPlusReadsFilePaths:
        correspondingMinusPath = tXRSeqPlusReadsFilePath.rsplit("+.bigWig")[0] + "-.bigWig"
        if not os.path.exists(correspondingMinusPath):
            raise ValueError("Corresponding minus file path not found at " + correspondingMinusPath)
        else:
            tXRSeqReadsFilePathPairs.append((tXRSeqPlusReadsFilePath, correspondingMinusPath))

    parseTXRSeq(tXRSeqReadsFilePathPairs, humanGenomeFastaFilePath)