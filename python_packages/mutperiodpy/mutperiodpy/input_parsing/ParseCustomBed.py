# This script takes a pseudo-bed file of the following format:
#   Col 0: chromosome identifier, e.g. "chr1" NOT just a single character.
#   Col 1: 0 base start pos
#   Col 2: 1 base end pos
#   Col 3: The base in the reference genome at this position. (Will be auto-acquired if set to ".")
#   Col 4: The base that the position was mutated to, or one of the following identifiers: 
#          INS: An insertion at the given position.
#          DEL: A deltion of the given region.
#          OTHER: Any other lesion or feature
#   Col 5: The strand the mutation/alteration occurred in.  Single-base mutations are flipped so that they occur in 
#          the pyrimidine-containing strand if necessary.  (If set to ".", the strand is first determined from the genome file)
#   Col 6: The chort the tumor belongs to.  e.g. a donor ID or tumor type.  Optional, but required for stratifying 
#          data in future steps.  If any cohort designations are given, ALL entries must have designations.  A "." character
#          in this column can be used to avoid assigning an entry to another cohort without breaking this rule.
# The file is then converted to a format suitable for the rest of the package's analysis scripts.

import os, subprocess
from typing import List
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog, Selections
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (dataDirectory, getIsolatedParentDir, generateMetadata,
                                                                  generateFilePath, dataTypes, Metadata, checkDirs)
from mutperiodpy.helper_scripts.UsefulBioinformaticsFunctions import (bedToFasta, FastaFileIterator, baseChromosomes,
                                                                      isPurine, reverseCompliment)
from mutperiodpy.input_parsing.WriteManager import WriteManager
from mutperiodpy.input_parsing.IdentifyMSI import MSIIdentifier


# Checks for common errors in a line of input.
def checkForErrors(choppedUpLine: List[str], cohortDesignationPresent):

    if len(choppedUpLine) < 6 or len(choppedUpLine) > 8:
        raise ValueError("Entry given with invalid number of arguments (" + str(len(choppedUpLine)) + ").  " +
                         "Each entry should contain either 6 tab separated arguments or 7 (if a cohort designation is present.")

    if choppedUpLine[0] not in baseChromosomes:
        raise ValueError("Invalid chromosome identifier \"" + choppedUpLine[0] + "\" found.  " + 
                         "Expected chromosome in the set at: NEED TO LINK TO ACCEPTABLE CHROMOSOME FILE")

    if not (choppedUpLine[1].isnumeric() and choppedUpLine[2].isnumeric()) or int(choppedUpLine[1]) >= int(choppedUpLine[2]):
        raise ValueError("Base positions should be integers and specify a minimum range of 1 base.  " +
                         "The first base coordinate should be 0-based and the second should be 1-based (inclusive).  " +
                         "The given coordinate pair, " + choppedUpLine[1] + ", " + choppedUpLine[2] + ' ' +
                         "does not satisfy these conditions.")

    if choppedUpLine[3].upper() not in ('A','C','G','T','.'):
        raise ValueError("Invalid reference base: \"" + choppedUpLine[3] + "\".  Should be one of the four DNA bases or \".\" " + 
                         "to auto acquire it from the corresponding genome.")

    if choppedUpLine[4].upper() not in ('A','C','G','T',"INS","DEL","OTHER"):
        raise ValueError("Invalid altered base: \"" + choppedUpLine[4] + "\".  Should be one of the four DNA bases or " + 
                         "\"INS\", \"DEL\", or \"OTHER\".")

    if choppedUpLine[5] not in ('+','-','.'):
        raise ValueError("Invalid strand designation: \"" + choppedUpLine[5] + "\".  Should be \"+\" or \"-\" or \".\" " + 
                         "to auto acquire it based on pyrimidine containing strand (if single-base mutation) " + 
                         "or set it to \"+\" if not.")

    if cohortDesignationPresent and len(choppedUpLine) == 6:
        raise ValueError("A cohort designation was given for previous entries, but an entry was found without one.  " + 
                         "The designation should either be always present or always absent in the input file.  " +
                         "Use \".\" to denote an entry that does not belong to a cohort.")

    if cohortDesignationPresent and len(choppedUpLine[6].strip()) == 0:
        raise ValueError("Cohort designations should contain at least one non-whitespace character.  Use \".\" to denote " + 
                         "an entry that does not belong to a cohort.")

    if not cohortDesignationPresent and len(choppedUpLine) == 7:
        raise ValueError("No cohort designation was given for previous entries, but an entry was found with one.  " + 
                         "The designation should either be always present or always absent in the input file.  " +
                         "Use \".\" to denote an entry that does not belong to a cohort.")


def equivalentEntries(fastaEntry: FastaFileIterator.FastaEntry, choppedUpLine):

    return(fastaEntry.chromosome == choppedUpLine[0] and
           fastaEntry.startPos == choppedUpLine[1] and
           fastaEntry.endPos == choppedUpLine[2] and
           fastaEntry.strand == choppedUpLine[5])


# Converts the standard bed input into the singlenuc context format acceptable for further analysis.
# Auto-acquire sequences if necessary.
def convertToStandardInput(bedInputFilePath, genomeFilePath, autoAcquiredFilePath, writeManager: WriteManager, inputQAChecked = False):

    print("Converting custom bed file to standard bed input...")

    # To start, assume that no sequences need to be acquired, and do it on the fly if need be.
    autoAcquiring = False
    autoAcquireFastaIterator = None
    fastaEntry = None

    cohortDesignationPresent = None # Keeps track of whether cohort IDs are given in the data.
    currentCohort = None # Keeps track of the current cohort

    # Iterate through the input file one line at a time, converting each line to an acceptable format for the rest of the pipeline.
    with open(bedInputFilePath,'r') as bedInputFile:
        for line in bedInputFile:

            choppedUpLine = str(line).strip().split('\t')

            # If it isn't already, initialize the cohortDesignationPresent variable.
            if cohortDesignationPresent is None: cohortDesignationPresent = len(choppedUpLine) == 7

            # Check for some possible error states.
            if not inputQAChecked: checkForErrors(choppedUpLine, cohortDesignationPresent)

            # If this is the first entry requiring auto-acquiring, generate the required fasta file.
            if not autoAcquiring and (choppedUpLine[3] == '.' or choppedUpLine[5] == '.'):
                print("Found line with auto-acquire requested.  Generating fasta...")
                autoAcquiring = True
                bedToFasta(bedInputFilePath, genomeFilePath, autoAcquiredFilePath)
                autoAcquiredFile = open(autoAcquiredFilePath, 'r')
                autoAcquireFastaIterator = FastaFileIterator(autoAcquiredFile)
                fastaEntry = autoAcquireFastaIterator.readEntry()

            # Check for any base identities that need to be auto-acquired.
            if choppedUpLine[3] == '.':

                # Find the equivalent fasta entry.
                while not equivalentEntries(fastaEntry, choppedUpLine):
                    if autoAcquireFastaIterator.eof:
                        raise ValueError("Reached end of fasta file without finding a match for: ",' '.join(choppedUpLine))
                    fastaEntry = autoAcquireFastaIterator.readEntry()

                # Set the sequence.
                choppedUpLine[3] = fastaEntry.sequence

            # Check for any strand designations that need to be auto-acquired.
            if choppedUpLine[5] == '.':

                # Find the equivalent fasta entry.
                while not equivalentEntries(fastaEntry, choppedUpLine):
                    if autoAcquireFastaIterator.eof:
                        raise ValueError("Reached end of fasta file without finding a match for: ",' '.join(choppedUpLine))
                    fastaEntry = autoAcquireFastaIterator.readEntry()

                if fastaEntry.sequence == choppedUpLine[3]: choppedUpLine[5] = '+'
                elif fastaEntry.sequence == reverseCompliment(choppedUpLine[3]): choppedUpLine[5] = '-'
                else: raise ValueError("The given sequence " + choppedUpLine[3] + " for location " + fastaEntry.sequenceName + ' ' +
                                       "does not match the corresponding sequence in the given genome, or its reverse compliment.")

            # Is this an SNP from a purine?  If so, flip the strand to the pyrimidine containing strand.
            if len(choppedUpLine[4]) == 1 and len(choppedUpLine[3]) == 1 and isPurine(choppedUpLine[3]):

                choppedUpLine[3] = reverseCompliment(choppedUpLine[3])
                choppedUpLine[4] = reverseCompliment(choppedUpLine[4])
                if choppedUpLine[5] == '+': choppedUpLine[5] = '-'
                elif choppedUpLine[5] == '-': choppedUpLine[5] = '+'

            # Call on the write manager to handle the rest!
            if cohortDesignationPresent:
                writeManager.writeData(choppedUpLine[0], choppedUpLine[1], choppedUpLine[2],
                                       choppedUpLine[3], choppedUpLine[4], choppedUpLine[5],
                                       choppedUpLine[6])
            else:
                writeManager.writeData(choppedUpLine[0], choppedUpLine[1], choppedUpLine[2],
                                       choppedUpLine[3], choppedUpLine[4], choppedUpLine[5])


# Set up the WriteManager to stratify by microsatellite stability by idntifiying MSI cohorts.
def setUpForMSStratification(writeManager, bedInputFilePath, inputQAChecked):

    print("Prepping data for MSIseq...")

    # Get the MSIIdentifier from the write manager and complete its process.
    with writeManager.setUpForMSStratification() as myMSIIdentifier:
        with open(bedInputFilePath, 'r') as bedInputFile:

            for line in bedInputFile:

                choppedUpLine: List[str] = line.strip().split('\t')
                if not inputQAChecked: checkForErrors(choppedUpLine, True)

                if choppedUpLine[4] != "OTHER":

                    if len(choppedUpLine[4]) == 1 and len(choppedUpLine[3]) == 1:
                        choppedUpLine[4] = "SNP"

                else: continue

                myMSIIdentifier.addData(choppedUpLine[0], choppedUpLine[1], choppedUpLine[2],
                                        choppedUpLine[4], choppedUpLine[6])

        myMSIIdentifier.identifyMSICohorts()


# Handles the scripts main functionality.
def parseCustomBed(bedInputFilePaths, genomeFilePath, nucPosFilePath, stratifyByMS, separateIndividualCohorts):

    for bedInputFilePath in bedInputFilePaths:

        print("\nWorking in:",os.path.basename(bedInputFilePath))

        # Get some important file system paths for the rest of the function and generate metadata
        # If this is an intermediate file, keep in mind that it's not in the data group's root directory 
        # and metadata should already have been generated elsewhere
        if getIsolatedParentDir(bedInputFilePath) == "intermediate_files":
            dataDirectory = os.path.dirname(os.path.dirname(bedInputFilePath))
        else:
            dataDirectory = os.path.dirname(bedInputFilePath)
            generateMetadata(os.path.basename(dataDirectory), getIsolatedParentDir(genomeFilePath), getIsolatedParentDir(nucPosFilePath), 
                             os.path.basename(bedInputFilePath), os.path.dirname(bedInputFilePath))

        metadata = Metadata(dataDirectory)
        intermediateFilesDir = os.path.join(dataDirectory,"intermediate_files")
        checkDirs(intermediateFilesDir)
        autoAcquiredFilePath = os.path.join(intermediateFilesDir,"autoAcquire.fa")

        inputQAChecked = False # Whether or not the input bed file has been checked for common errors.

        # Create an instance of the WriteManager to handle writing.
        with WriteManager(dataDirectory) as writeManager:

            # Check to see if cohort designations are present to see if preparations need to be made.
            optionalArgument = tuple()
            with open(bedInputFilePath, 'r') as bedInputFile:
                line = bedInputFile.readline()

                # Is the cohort designation present?
                if len(line.strip().split('\t')) == 7: 

                    # Include in sort function
                    optionalArgument = ("-k7,7",)

                    # Prepare the write manager for individual cohorts if desired.
                    if separateIndividualCohorts: writeManager.setUpForIndividualCohorts()

                elif stratifyByMS: 
                    raise ValueError("Stratification by microsatellite stability requested, but no cohort designation given.")
                elif separateIndividualCohorts:
                    raise ValueError("Separation by individual cohorts requested, but no cohort designation given.")

            # Sort the input data (should also ensure that the output data is sorted)
            subprocess.run(" ".join(("sort",) + optionalArgument + 
                                    ("-k1,1","-k2,2n",bedInputFilePath,"-o",bedInputFilePath)), shell = True, check = True)

            # If requested, also prepare for stratification by microsatellite stability.
            if stratifyByMS:             
                setUpForMSStratification(writeManager, bedInputFilePath, inputQAChecked)
                inputQAChecked = True

            # Go, go, go!
            convertToStandardInput(bedInputFilePath, genomeFilePath, autoAcquiredFilePath, writeManager, inputQAChecked)


if __name__ == "__main__":

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=dataDirectory)
    dialog.createMultipleFileSelector("Custom bed Input Files:",0,"custom_input.bed",("bed files",".bed"))
    dialog.createFileSelector("Genome Fasta File:",1,("Fasta Files",".fa"))
    dialog.createFileSelector("Strongly Positioned Nucleosome File:",2,("Bed Files",".bed"))
    dialog.createCheckbox("Stratify data by microsatellite stability?", 3, 0)
    dialog.createCheckbox("Separate individual cohorts?", 3, 1)
    dialog.createExitButtons(4,0)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    bedInputFilePaths = list(selections.getFilePathGroups())[0]
    genomeFilePath = list(selections.getIndividualFilePaths())[0]
    nucPosFilePath = list(selections.getIndividualFilePaths())[1]
    stratifyByMS = list(selections.getToggleStates())[0]
    separateIndividualCohorts = list(selections.getToggleStates())[1]

    parseCustomBed(bedInputFilePaths, genomeFilePath, nucPosFilePath, stratifyByMS, separateIndividualCohorts)