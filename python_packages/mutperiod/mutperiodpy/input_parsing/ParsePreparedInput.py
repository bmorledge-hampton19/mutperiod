# This script is used when a file is already in an acceptable form of input for the rest of the mutperiod pipeline.
# Metadata is generated for the file and a few checks are performed to make sure the file is actually in the appropriate format.

import os
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (getDataDirectory, getIsolatedParentDir, generateMetadata, checkDirs, 
                                                                  InputFormat, DataTypeStr, getContext, getAcceptableChromosomes, InputError)
from mutperiodpy.input_parsing.ParseCustomBed import checkForErrors


def parsePreparedInput(inputFilePaths: List[str], genomeFilePath, checkEachLine = True):
    
    for inputFilePath in inputFilePaths:

        print("\nWorking in",os.path.basename(inputFilePath))

        # Perform some checks to make sure the input is formatted correctly.
        dataGroupName = getIsolatedParentDir(inputFilePath)
        inputFileBasename = os.path.basename(inputFilePath)
        inputFileContext = getContext(inputFilePath)

        if inputFileContext is None: raise InputError("No context is apparent from the given prepared input file.")
        if inputFileBasename.split('_'+inputFileContext)[0] != dataGroupName:
            raise InputError("Prepared input file is not named as expected given the data group name generated from the "
                             "parent directory.  Expected: \"" + dataGroupName + "\" immediately preceding the context definition")
        if not inputFileBasename.endswith(DataTypeStr.mutations + ".bed"):
            raise InputError("Prepared input file is not named properly to indicate the presence of mutation data.  "
                             "Expected a file ending in \"" + DataTypeStr.mutations + ".bed\"")

        acceptableChromosomes = getAcceptableChromosomes(genomeFilePath)
        acceptableChromosomesFilePath = getAcceptableChromosomes(genomeFilePath, True)

        # Perform QA with the checkForErrors function
        print("Checking for errors in line formatting...")
        with open(inputFilePath, 'r') as inputFile:
            try:

                choppedUpLine = inputFile.readline().strip().split('\t')
                cohortDesignationPresent = len(choppedUpLine) == 7
                checkForErrors(choppedUpLine, cohortDesignationPresent, acceptableChromosomes,
                               acceptableChromosomesFilePath)

                if checkEachLine:
                    for line in inputFile:
                        choppedUpLine = line.strip().split('\t')
                        checkForErrors(choppedUpLine, cohortDesignationPresent, acceptableChromosomes,
                                       acceptableChromosomesFilePath)

            except AssertionError as assertionError:
                raise InputError(str(assertionError))

        # If everything else looks good, generate the metadata.  This directory is now ready to go!
        generateMetadata(dataGroupName, getIsolatedParentDir(genomeFilePath), 
                         os.path.basename(inputFilePath), InputFormat.prepared,  os.path.dirname(inputFilePath))


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