# This script is used when a file is already in an acceptable form of input for the rest of the mutperiod pipeline.
# Metadata is generated for the file and a few checks are performed to make sure the file is actually in the appropriate format.

import os, subprocess
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (getDataDirectory, getIsolatedParentDir, generateMetadata, Metadata, 
                                                                  InputFormat, DataTypeStr, getContext, getAcceptableChromosomes)
from mutperiodpy.input_parsing.ParseCustomBed import checkForErrors
from benbiohelpers.CustomErrors import *


def parsePreparedInput(inputFilePaths: List[str], genomeFilePath, checkEachLine = True):
    
    for inputFilePath in inputFilePaths:

        print("\nWorking in",os.path.basename(inputFilePath))

        # Perform some checks to make sure the input is formatted correctly.
        dataGroupName = getIsolatedParentDir(inputFilePath)
        inputFileBasename = os.path.basename(inputFilePath)
        inputFileContext = getContext(inputFilePath)

        if inputFileContext is None: raise UserInputError("No context is apparent from the given prepared input file.")
        if inputFileBasename.split('_'+inputFileContext)[0] != dataGroupName:
            raise InvalidPathError(inputFilePath, 
                                   "Prepared input file is not named as expected given the data group name generated from the "
                                   "parent directory.  Expected: \"" + dataGroupName + "\" immediately preceding the context definition "
                                   "but given file path is:")
        if not inputFileBasename.endswith(DataTypeStr.mutations + ".bed"):
            raise InvalidPathError(inputFilePath,
                                   "Prepared input file is not named properly to indicate the presence of mutation data.  "
                                   "Expected a file ending in \"" + DataTypeStr.mutations + ".bed\" but given path is:")

        acceptableChromosomes = getAcceptableChromosomes(genomeFilePath)
        acceptableChromosomesFilePath = getAcceptableChromosomes(genomeFilePath, True)

        # Perform QA with the checkForErrors function
        print("Checking for errors in line formatting...")
        with open(inputFilePath, 'r') as inputFile:
            choppedUpLine = inputFile.readline().strip().split('\t')
            cohortDesignationPresent = len(choppedUpLine) == 7
            checkForErrors(choppedUpLine, cohortDesignationPresent, acceptableChromosomes,
                            acceptableChromosomesFilePath)

            if checkEachLine:
                for line in inputFile:
                    choppedUpLine = line.strip().split('\t')
                    checkForErrors(choppedUpLine, cohortDesignationPresent, acceptableChromosomes,
                                    acceptableChromosomesFilePath)

        # If everything else looks good, generate the metadata.  This directory is now ready to go!
        print("Checks passed.  Generating metadata, including mutation counts using a call to wc -l")
        generateMetadata(dataGroupName, getIsolatedParentDir(genomeFilePath), 
                         os.path.basename(inputFilePath), InputFormat.prepared,  os.path.dirname(inputFilePath))
        featureCounts = int(subprocess.check_output(("wc", "-l", inputFilePath), encoding = "UTF-8").split()[0])
        Metadata(inputFilePath).addMetadata(Metadata.AddableKeys.mutCounts, featureCounts)


def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Prepared Input Files:",0,DataTypeStr.mutations + ".bed",("bed files",".bed"))
    dialog.createFileSelector("Genome Fasta File:",1,("Fasta Files",".fa"))
    dialog.createCheckbox("Skip formatting checks for all but the first line",2,0)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    parsePreparedInput(dialog.selections.getFilePathGroups()[0], dialog.selections.getIndividualFilePaths()[0],
                       not dialog.selections.getToggleStates()[0])

if __name__ == "__main__": main()