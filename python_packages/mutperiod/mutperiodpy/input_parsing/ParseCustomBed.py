# This script takes a pseudo-bed file of the following format:
#   Col 0: chromosome identifier, e.g. "chr1" NOT just a single character.
#   Col 1: 0 base start pos
#   Col 2: 1 base end pos
#   Col 3: The base(s) in the reference genome at this position. 
#          Will be auto-acquired if set to "."
#          "*" indicates an insertion between the two given bases.
#   Col 4: The base(s) that the position(s) were mutated to.
#          "*" indicates a deletion of the given base
#          Use "OTHER" for any other lesion or feature ("." will be coerced to "OTHER")
#   Col 5: The strand the mutation/alteration occurred in.  Single-base mutations are flipped so that they occur in 
#          the pyrimidine-containing strand if necessary.  
#          If set to ".", the strand is first determined from the genome file, if possible (not an insertion).
#   Col 6: The chort the tumor belongs to.  e.g. a donor ID or tumor type.  Optional, but required for stratifying 
#          data in future steps.  If any cohort designations are given, ALL entries must have designations.  A "." character
#          in this column can be used to avoid assigning an entry to another cohort without breaking this rule.
# The file is then converted to a format suitable for the rest of the package's analysis scripts.

import os, subprocess, sys, shutil
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog, Selections
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (getDataDirectory, getIsolatedParentDir, generateMetadata, checkDirs, 
                                                                  DataTypeStr, InputFormat, getAcceptableChromosomes, generateFilePath,
                                                                  rScriptsDirectory)
from mutperiodpy.helper_scripts.CustomErrors import *
from benbiohelpers.FileSystemHandling.DirectoryHandling import getFilesInDirectory
from benbiohelpers.FileSystemHandling.BedToFasta import bedToFasta
from benbiohelpers.FileSystemHandling.FastaFileIterator import FastaFileIterator
from benbiohelpers.DNA_SequenceHandling import isPurine, reverseCompliment
from mutperiodpy.input_parsing.WriteManager import WriteManager
from benbiohelpers.CustomErrors import *


# Checks for common errors in a line of input.
def checkForErrors(choppedUpLine: List[str], cohortDesignationPresent, acceptableChromosomes, acceptableChromosomesFilePath):

    if len(choppedUpLine) < 6 or len(choppedUpLine) > 7:
        raise UserInputError("Entry given with invalid number of arguments (" + str(len(choppedUpLine)) + ").  "
                         "Each entry should contain either 6 tab separated arguments or 7 (if a cohort designation is present. "
                         f"The line in question is:\n{chr(9).join(choppedUpLine)}")

    if choppedUpLine[0] not in acceptableChromosomes:
        raise UserInputError("Invalid chromosome identifier \"" + choppedUpLine[0] + "\" found.  "
                         "Expected chromosome in the set at: " + acceptableChromosomesFilePath)

    try:
        float(choppedUpLine[1])
        float(choppedUpLine[2])
    except ValueError:
        raise UserInputError("Base positions should be numeric values.  "
                         "The given coordinate pair, " + choppedUpLine[1] + ", " + choppedUpLine[2] + ' '
                         "does not satisfy this condition.")

    if float(choppedUpLine[2]) - float(choppedUpLine[1]) < 1:
        raise UserInputError("Base positions should specify a minimum range of 1 base.  "
                         "The first base coordinate should be 0-based and the second should be 1-based.  " 
                         "The given coordinate pair, " + choppedUpLine[1] + ", " + choppedUpLine[2] + ' '
                         "does not satisfy these conditions.")

    if not {'A','C','G','T','N'}.issuperset(choppedUpLine[3]) and not choppedUpLine[3] in ('*','.'):
        raise UserInputError("Invalid reference base(s): \"" + choppedUpLine[3] + "\".  Should be a string made up of the four DNA bases "
                         "or \".\" to auto acquire the base(s) from the corresponding genome, or \"*\" to denote an insertion.")

    regionLength = float(choppedUpLine[2]) - float(choppedUpLine[1])

    # Cases where insertion...
    if choppedUpLine[3] == '*':

        if choppedUpLine[4] == '*':
            raise UserInputError("The reference base and mutation columns are both set to \"*\".  "
                             "(Entry cannot be insertion and deletion simultaneously)")

        if regionLength != 2:
            raise UserInputError("Reference column indicates an insertion, but the genome positions don't exactly flank the insertion site.  "
                             "(Indicated region should span exactly 2 bases.)")

    # Cases where not insertion...
    else:

        if choppedUpLine[3] != '.' and len(choppedUpLine[3]) < regionLength:
            raise UserInputError("References base(s): \"" + choppedUpLine[3] + "\" do not at least span the given region from "
                             "positions " + choppedUpLine[1] + " to " + choppedUpLine[2] + '.')

    if not {'A','C','G','T'}.issuperset(choppedUpLine[4]) and not choppedUpLine[4] in ('*',"OTHER", '.'):
        raise UserInputError("Invalid mutation designation: \"" + choppedUpLine[3] + "\".  Should be a string made up of the four DNA bases "
                         "or \"*\" to denote a deletion, or \"OTHER\" or \'.\' to denote an some other alteration.")

    if choppedUpLine[5] not in ('+','-','.'):
        raise UserInputError("Invalid strand designation: \"" + choppedUpLine[5] + "\".  Should be \"+\" or \"-\", or \".\" "
                         "to auto acquire the base if possible.")

    if cohortDesignationPresent and len(choppedUpLine) == 6:
        raise UserInputError("A cohort designation was given for previous entries, but an entry was found without one.  "
                         "The designation should either be always present or always absent in the input file.  "
                         "Use \".\" to denote an entry that does not belong to a cohort.")

    if cohortDesignationPresent and len(choppedUpLine[6].strip()) == 0:
        raise UserInputError("Cohort designations should contain at least one non-whitespace character.  Use \".\" to denote "
                         "an entry that does not belong to a cohort.")

    if not cohortDesignationPresent and len(choppedUpLine) == 7:
        raise UserInputError("No cohort designation was given for previous entries, but an entry was found with one.  "
                         "The designation should either be always present or always absent in the input file.  "
                         "Use \".\" to denote an entry that does not belong to a cohort.")


def equivalentEntries(fastaEntry: FastaFileIterator.FastaEntry, choppedUpLine):

    return(fastaEntry.chromosome == choppedUpLine[0] and
           fastaEntry.startPos == choppedUpLine[1] and
           fastaEntry.endPos == choppedUpLine[2] and
           fastaEntry.strand == choppedUpLine[5])


def isSingleBaseSubstitution(choppedUpLine):
    return float(choppedUpLine[2]) - float(choppedUpLine[1]) == 1 and choppedUpLine[4] in ('A','C','G','T')

# Checks each line for errors and auto acquire bases/strand designations where requested. 
# Overwrites the original bed file if auto-acquiring occurred.
# Also returns the numerical nucleotide context of the features.
def autoAcquireAndQACheck(bedInputFilePath: str, genomeFilePath, autoAcquiredFilePath, onlySingleBaseSubs, includeIndels):

    print("Checking custom bed file for formatting and auto-acquire requests...")

    # To start, assume that no sequences need to be acquired, and do it on the fly if need be.
    autoAcquiring = False
    autoAcquireFastaIterator = None
    fastaEntry = None
    cohortDesignationPresent = None

    # Unless indels are included, determine the context of the feqtures in the file.
    if includeIndels: context = 0
    else: context = None

    # Get the list of acceptable chromosomes
    acceptableChromosomes = getAcceptableChromosomes(genomeFilePath)
    acceptableChromosomesFilePath = getAcceptableChromosomes(genomeFilePath, True)

    # Create a temporary file to write the data to (potentially after auto-acquiring).  
    # Will replace original file at the end if auto-acquiring occurred.
    temporaryBedFilePath = bedInputFilePath + ".tmp"

    # Iterate through the input file one line at a time, checking the format of each entry and looking for auto-acquire requests.
    with open(bedInputFilePath, 'r') as bedInputFile:
        with open(temporaryBedFilePath, 'w')as temporaryBedFile:
            for line in bedInputFile:

                choppedUpLine = str(line).strip().split('\t')

                # If it isn't already, initialize the cohortDesignationPresent variable.
                if cohortDesignationPresent is None: cohortDesignationPresent = len(choppedUpLine) == 7

                # Check for possible error states.
                checkForErrors(choppedUpLine, cohortDesignationPresent, acceptableChromosomes, acceptableChromosomesFilePath)

                # If this is the first entry requiring auto-acquiring, generate the required fasta file.
                if (not autoAcquiring and (choppedUpLine[3] == '.' or choppedUpLine[4] == '.' or 
                    (choppedUpLine[5] == '.' and choppedUpLine[3] != '*'))):
                    print("Found line with auto-acquire requested.  Generating fasta...")
                    autoAcquiring = True
                    failedAutoAcquire = False
                    try:
                        bedToFasta(bedInputFilePath, genomeFilePath, autoAcquiredFilePath)
                    except subprocess.CalledProcessError:
                        failedAutoAcquire = True
                    if failedAutoAcquire:
                        raise UserInputError("The given custom bed file has lines with auto-acquire requested "
                                             "(symbolized by a '.' character), but a fasta file could not be generated from "
                                             "the given bed file.")
                    autoAcquiredFile = open(autoAcquiredFilePath, 'r')
                    autoAcquireFastaIterator = FastaFileIterator(autoAcquiredFile)
                    fastaEntry = autoAcquireFastaIterator.readEntry()
                    print("Continuing...")

                # Check for any base identities that need to be auto-acquired.
                if choppedUpLine[3] == '.':

                    # Find the equivalent fasta entry.
                    while not equivalentEntries(fastaEntry, choppedUpLine):
                        assert not autoAcquireFastaIterator.eof, (
                            "Reached end of fasta file without finding a match for: ",' '.join(choppedUpLine))
                        fastaEntry = autoAcquireFastaIterator.readEntry()

                    # Set the sequence.
                    choppedUpLine[3] = fastaEntry.sequence

                # Check for any strand designations that need to be auto-acquired.
                # Also, make sure this isn't an insertion, in which case the strand designation cannot be determined.
                if choppedUpLine[5] == '.' and choppedUpLine[3] != '*':

                    # Find the equivalent fasta entry.
                    while not equivalentEntries(fastaEntry, choppedUpLine):
                        assert not autoAcquireFastaIterator.eof, (
                            "Reached end of fasta file without finding a match for: ",' '.join(choppedUpLine))
                        fastaEntry = autoAcquireFastaIterator.readEntry()

                    # Determine which strand is represented.
                    if fastaEntry.sequence == choppedUpLine[3]: choppedUpLine[5] = '+'
                    elif fastaEntry.sequence == reverseCompliment(choppedUpLine[3]): choppedUpLine[5] = '-'
                    else: assert False, ("The given sequence " + choppedUpLine[3] + " for location " + fastaEntry.sequenceName + ' ' +
                                         "does not match the corresponding sequence in the given genome, or its reverse compliment.")

                # Change any '.' characters in the "altered to" column to "OTHER"
                if choppedUpLine[4] == '.': choppedUpLine[4] = "OTHER"

                # Determine the sequence context of the line and whether or not it matches the sequence context for other.
                # Skip this if the file is "mixed", this line is an indel, or only single base substitutions are allowed and this line isn't one.
                if (not context == 0 and 
                    not (choppedUpLine[3] == '*' or choppedUpLine[4] == '*') and 
                    (not onlySingleBaseSubs or isSingleBaseSubstitution(choppedUpLine))):

                    thisContext = len(choppedUpLine[3])
                    if context is None: context = thisContext
                    elif thisContext != context: context = 0

                # Write the current line to the temporary bed file.
                temporaryBedFile.write('\t'.join(choppedUpLine) + '\n')

    # If any lines were auto-acquired, replace the input bed file with the temporary bed file. (Which has auto-acquires)
    if autoAcquiring:
        print("Overwriting custom bed input with auto-acquired bases/strand designations.")
        os.replace(temporaryBedFilePath, bedInputFilePath)
    # Otherwise, just delete the temporary file.
    else: os.remove(temporaryBedFilePath)

    if context > 6: context = float("inf")
    return context
    

# Converts the custom bed input into the singlenuc context format acceptable for analysis further down the pipeline.
def convertToStandardInput(bedInputFilePath, writeManager: WriteManager, onlySingleBaseSubs, includeIndels):

    print("Converting custom bed file to standard bed input...")

    # Iterate through the input file one line at a time, converting each line to an acceptable format for the rest of the pipeline.
    with open(bedInputFilePath,'r') as bedInputFile:
        for line in bedInputFile:
            
            choppedUpLine = str(line).strip().split('\t')

            # Is this an SNP from a purine?  If so, flip the strand to the pyrimidine containing strand.
            if isPurine(choppedUpLine[3]) and choppedUpLine[4].upper() in ('A','C','G','T'):

                choppedUpLine[3] = reverseCompliment(choppedUpLine[3])
                choppedUpLine[4] = reverseCompliment(choppedUpLine[4])
                if choppedUpLine[5] == '+': choppedUpLine[5] = '-'
                elif choppedUpLine[5] == '-': choppedUpLine[5] = '+'

            # Is this an indel, and are those included?
            if choppedUpLine[3] == '*' or choppedUpLine[4] == '*':
                if not includeIndels: continue

            # Is this a single base substitution, and if not, should it even be included?
            if not isSingleBaseSubstitution(choppedUpLine):
                if onlySingleBaseSubs: continue

                # Center features greater than a single nucleotide so that they occur at a single nucleotide position (or half position)
                else: 
                    center = (float(choppedUpLine[1])+float(choppedUpLine[2])-1)/2
                    if int(center) == center: center = int(center) # Remove the decimal from the float if possible.
                    choppedUpLine[1] = str(center)
                    choppedUpLine[2] = str(center + 1)
            

            # Call on the write manager to handle the rest!
            if len(choppedUpLine) == 7:
                writeManager.writeData(choppedUpLine[0], choppedUpLine[1], choppedUpLine[2],
                                       choppedUpLine[3], choppedUpLine[4], choppedUpLine[5],
                                       choppedUpLine[6])
            else:
                writeManager.writeData(choppedUpLine[0], choppedUpLine[1], choppedUpLine[2],
                                       choppedUpLine[3], choppedUpLine[4], choppedUpLine[5])


# Set up the WriteManager to stratify by microsatellite stability by identifiying MSI cohorts.
def setUpForMSStratification(writeManager: WriteManager, bedInputFilePath):

    print("Prepping data for MSIseq...")
    
    # Get the MSIIdentifier from the write manager and complete its process.
    with writeManager.setUpForMSStratification() as myMSIIdentifier:
        with open(bedInputFilePath, 'r') as bedInputFile:

            mutType = None # Either SNP, INS, or DEL

            for line in bedInputFile:

                choppedUpLine: List[str] = line.strip().split('\t')

                if choppedUpLine[4] != "OTHER" and choppedUpLine[6] != '.':

                    if choppedUpLine[3] == '*':
                        mutType = "INS"
                    elif choppedUpLine[4] == '*':
                        mutType = "DEL"
                    else: mutType = "SNP"

                else: continue
                
                myMSIIdentifier.addData(choppedUpLine[0], str(int(choppedUpLine[1]) + 1), choppedUpLine[2],
                                        mutType, choppedUpLine[6])

        myMSIIdentifier.identifyMSICohorts()


# Set up the WriteManager to stratify by mutation signature by assigning mutation signatures to cohorts.
def setUpForMutSigStratification(writeManager: WriteManager, bedInputFilePath):

    print("Prepping data for deconstructSigs...")

    # Get the MutSigIdentifier from the write manager and complete its process.
    with writeManager.setUpForMutSigStratification() as mutSigIdentifier:
        with open(bedInputFilePath, 'r') as bedInputFile:

            for line in bedInputFile:

                choppedUpLine: List[str] = line.strip().split('\t')

                # Only use SNP's.  Skip this entry if it does not represent an SNP.
                if choppedUpLine[3] not in ('A','C','G','T') or choppedUpLine[4] not in ('A','C','G','T'): continue

                # Skip any entries that are independent of cohorts.
                if choppedUpLine[6] == '.': continue

                mutSigIdentifier.addData(choppedUpLine[6], choppedUpLine[0], choppedUpLine[2], 
                                         choppedUpLine[3], choppedUpLine[4])

        mutSigIdentifier.identifyMutSigs()


# Handles the scripts main functionality.
def parseCustomBed(bedInputFilePaths, genomeFilePath,
                   stratifyByMS = False, stratifyByMutSig = False, separateIndividualCohorts = False,
                   onlySingleBaseSubs = False, includeIndels = False,
                   simpleParsing = False):

    # This needs to be here to avoid a circular reference.
    from mutperiodpy.input_parsing.ParseStandardBed import parseStandardBed

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


    # If simple parsing was selected, do something kinda weird where the inputs get passed to
    # parseStandardBed first which eventually passes them back here.
    if simpleParsing: return parseStandardBed(bedInputFilePaths, genomeFilePath)

    if onlySingleBaseSubs and includeIndels: raise UserInputError("Indels are incompatible with single nucleotide substitutions.")
    if len(bedInputFilePaths) == 0: raise UserInputError("No bed files were found to parse.") 

    outputFilePaths = list()

    for bedInputFilePath in bedInputFilePaths:

        print("\nWorking in:",os.path.basename(bedInputFilePath))

        # Get some important file system paths for the rest of the function and generate metadata
        # If this is an intermediate file, keep in mind that it's not in the data group's root directory 
        # and metadata should already have been generated elsewhere
        if getIsolatedParentDir(bedInputFilePath) == "intermediate_files":
            dataDirectory = os.path.dirname(os.path.dirname(bedInputFilePath))
        else:
            dataDirectory = os.path.dirname(bedInputFilePath)
            generateMetadata(os.path.basename(dataDirectory), getIsolatedParentDir(genomeFilePath), 
                             os.path.basename(bedInputFilePath), InputFormat.customBed,  os.path.dirname(bedInputFilePath))

        intermediateFilesDir = os.path.join(dataDirectory,"intermediate_files")
        checkDirs(intermediateFilesDir)
        autoAcquiredFilePath = os.path.join(intermediateFilesDir,"auto_acquire.fa")

        context = autoAcquireAndQACheck(bedInputFilePath, genomeFilePath, autoAcquiredFilePath, onlySingleBaseSubs, includeIndels)

        # Make sure the input file is not named the same as what will become the output file.  If it is, it needs to be copied
        # to the intermediate_files directory so it is available to be read from as the new output file is being written.
        expectedOutputFilePath = generateFilePath(directory = dataDirectory, dataGroup = os.path.basename(dataDirectory),
                                                  context = context, dataType = DataTypeStr.mutations, fileExtension = ".bed")
        if bedInputFilePath == expectedOutputFilePath:
            inputFilePathCopy = os.path.join(intermediateFilesDir,os.path.basename(bedInputFilePath))
            print("Input file path is identical to generated output file path and will be overwritten. ",
                  "Creating a copy of the input file at:", inputFilePathCopy, "to use for reading.")
            shutil.copy2(bedInputFilePath, inputFilePathCopy)
            bedInputFilePath = inputFilePathCopy

        outputFilePaths.append(expectedOutputFilePath)

        # Create an instance of the WriteManager to handle writing.
        with WriteManager(dataDirectory, context) as writeManager:

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

                elif stratifyByMS or stratifyByMutSig: 
                    raise UserInputError("Additional stratification given, but no cohort designation given.")
                elif separateIndividualCohorts:
                    raise UserInputError("Separation by individual cohorts requested, but no cohort designation given.")

            # Sort the input data (should also ensure that the output data is sorted)
            subprocess.run(("sort",) + optionalArgument + 
                            ("-k1,1","-k2,2n", "-k3,3n",bedInputFilePath,"-s","-o",bedInputFilePath), check = True)

            # If requested, also prepare for stratification by microsatellite stability.
            if stratifyByMS:             
                setUpForMSStratification(writeManager, bedInputFilePath)

            if stratifyByMutSig:
                setUpForMutSigStratification(writeManager, bedInputFilePath)

            # Go, go, go!
            convertToStandardInput(bedInputFilePath, writeManager, onlySingleBaseSubs, includeIndels)

    return outputFilePaths


# Given a namespace resulting from an argparser object (constructed in mutperiodpy.Main),
# use the input to run this script.
def parseArgs(args):
    
    # If only the subcommand was given, run the UI.
    if len(sys.argv) == 2: 
        main(); return

    # Otherwise, check to make sure valid arguments were passed:
    if args.genome_file is None: raise UserInputError("No genome file was given.")
    checkIfPathExists(args.genome_file)
    genomeFilePath = os.path.abspath(args.genome_file)

    # If the coerce_bed argument has bee passed, make sure that no other incompatible arguments were passed.
    if args.coerce_bed and (args.stratify_by_microsatellite or args.stratify_by_mut_sigs or 
                            args.stratify_by_cohorts or args.only_sbs or args.include_indels):
        raise UserInputError("When coercing bed input to custom input, the following options cannot be used: "
                             "stratify-by-microsatellite, stratify-by-mut-sigs, stratify-by-cohorts, "
                             "only-sbs, include-indels")

    # Get the custom bed files from the given paths, searching directories if necessary.
    finalCustomBedPaths = list()
    for bedFilePath in args.bedFilePaths:
        if os.path.isdir(bedFilePath):
            finalCustomBedPaths += [os.path.abspath(filePath) for filePath in getFilesInDirectory(bedFilePath, "custom_input.bed")]
        elif checkIfPathExists(bedFilePath): finalCustomBedPaths.append(os.path.abspath(bedFilePath))

    # Run the parser.
    parseCustomBed(list(set(finalCustomBedPaths)), genomeFilePath, args.stratify_by_microsatellite, 
                   args.stratify_by_mut_sigs, args.stratify_by_cohorts, args.only_sbs, args.include_indels,
                   args.coerce_bed)


def main():

    FULL_CUSTOM = "FullCustom"

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory(), title = "Parse Custom Bed Data")
    dialog.createMultipleFileSelector("Custom bed Input Files:",0,"custom_input.bed",("bed files",".bed"))
    dialog.createFileSelector("Genome Fasta File:",1,("Fasta Files",".fa"))
    with dialog.createDynamicSelector(2, 0) as simplifiedBedDynSel:
        simplifiedBedDynSel.initCheckboxController("Input is bed3 (or bed6+ with strand info) and must be coereced to custom bed")
        fullCustomDialog = simplifiedBedDynSel.initDisplay(False, FULL_CUSTOM)
        fullCustomDialog.createCheckbox("Stratify data by microsatellite stability?", 0, 0)
        fullCustomDialog.createCheckbox("Stratify by mutation signature?", 0, 1)
        fullCustomDialog.createCheckbox("Separate individual cohorts?", 1, 0)
        fullCustomDialog.createCheckbox("Only use single nucleotide substitutions?", 2, 0)
        fullCustomDialog.createCheckbox("Include indels in output?", 2, 1)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    bedInputFilePaths = selections.getFilePathGroups()[0]
    genomeFilePath = selections.getIndividualFilePaths()[0]

    # If we have a simplified input, we don't need any other parameters and can start parsing now.
    if simplifiedBedDynSel.getControllerVar(): 
        parseCustomBed(bedInputFilePaths, genomeFilePath, simpleParsing = True)
    
    # Otherwise, set up the parser with the additional parameters.
    else:
        stratifyByMS = selections.getToggleStates(FULL_CUSTOM)[0]
        stratifyByMutSig = selections.getToggleStates(FULL_CUSTOM)[1]
        separateIndividualCohorts = selections.getToggleStates(FULL_CUSTOM)[2]
        onlySingleBaseSubs = selections.getToggleStates(FULL_CUSTOM)[3]
        includeIndels = selections.getToggleStates(FULL_CUSTOM)[4]

        parseCustomBed(bedInputFilePaths, genomeFilePath, stratifyByMS, 
                       stratifyByMutSig, separateIndividualCohorts, onlySingleBaseSubs, includeIndels)

if __name__ == "__main__": main()
