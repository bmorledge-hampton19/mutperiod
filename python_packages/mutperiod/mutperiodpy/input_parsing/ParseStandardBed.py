# This script takes standard bed formatted data and converts it to custom bed format.
# (Basically, it just replaces the 4th column with '.' and the 5th column with "OTHER" and also generates metadata.
# If a sixth column is present, it is kept as is, and the script assumes that it contains a strand designation.
# Otherwise, the strand is set to '+' in all rows.)
# The file is then passed along to ParseCustomBed to finish formatting for the mutperiod pipeline.

import os
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import Selections, TkinterDialog
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (getIsolatedParentDir, generateFilePath, getDataDirectory,
                                                                  DataTypeStr, generateMetadata, InputFormat, checkDirs,
                                                                  getAcceptableChromosomes)
from benbiohelpers.CustomErrors import *


def parseStandardBed(standardBedFilePaths: List[str], genomeFilePath):

    # This needs to be here to avoid a circular reference.
    from mutperiodpy.input_parsing.ParseCustomBed import parseCustomBed

    customBedOutputFilePaths = list() # The list of file paths to be passed to the custom bed parser.

    # Parse the given files into custom bed format.
    for standardBedFilePath in standardBedFilePaths:

        print("\nWorking in:",os.path.basename(standardBedFilePath))
        if not os.path.basename(standardBedFilePath).endswith(".bed"):
            raise InvalidPathError(standardBedFilePath, 
                                   "Given file does not appear to be in bed format. (missing \".bed\" extension)")

        # Store useful paths and names.
        localRootDirectory = os.path.dirname(standardBedFilePath)
        intermediateFilesDir = os.path.join(localRootDirectory,"intermediate_files")
        checkDirs(intermediateFilesDir)
        dataGroupName = getIsolatedParentDir(standardBedFilePath)

        # Generate the output file path and metadata
        customBedOutputFilePath = generateFilePath(directory = intermediateFilesDir, dataGroup = dataGroupName,
                                                   dataType = DataTypeStr.customInput, fileExtension = ".bed")
        customBedOutputFilePaths.append(customBedOutputFilePath)
        generateMetadata(dataGroupName, getIsolatedParentDir(genomeFilePath), 
                         os.path.basename(standardBedFilePath), InputFormat.standardBed, localRootDirectory)

        # Get the list of acceptable chromosomes.
        acceptableChromosomes = getAcceptableChromosomes(genomeFilePath)

        # Iterate through the standard bed file entries preparing them for custom-bed input.
        print("Converting entries for custom bed input...")
        with open(standardBedFilePath, 'r') as standardBedFile:
            with open(customBedOutputFilePath, 'w') as customBedOutputFile:

                for line in standardBedFile:
                    
                    choppedUpLine = line.strip().split("\t")

                    # Make sure we have at least a valid bed3 line
                    if len(choppedUpLine) < 3:
                        raise UserInputError(f"Found bed entry with less than 3 columns: \"{line.strip()}\"\n"
                                             "Bed entries must contain at least chromosome, start position, and end position.")

                    # Make sure the lesion is in a valid chromosome.  Otherwise, skip it.
                    if not choppedUpLine[0] in acceptableChromosomes: continue

                    # Ensure that the line has at least 6 columns by adding blank ones if necessary.
                    choppedUpLine += [''] * (6-len(choppedUpLine))

                    choppedUpLine[3] = '.'
                    choppedUpLine[4] = "OTHER"
                    if choppedUpLine[5] != '-': choppedUpLine[5] = '+'

                    customBedOutputFile.write('\t'.join(choppedUpLine[:6]) + '\n')


    # Pass the generated files to the custom bed parser.
    return parseCustomBed(customBedOutputFilePaths, genomeFilePath)


if __name__ == "__main__":

    # Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Standard Bed Data:",0,"dipy.bed",("Bed Files",".bed"),additionalFileEndings=("TA.bed",))    
    dialog.createFileSelector("Genome Fasta File:",1,("Fasta Files",".fa"))

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    standardBedFilePaths = list(selections.getFilePathGroups())[0]
    genomeFilePath = list(selections.getIndividualFilePaths())[0]

    parseStandardBed(standardBedFilePaths, genomeFilePath)
