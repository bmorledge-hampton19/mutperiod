# This script takes standard bed formatted data and converts it to custom bed format.
# (Basically, it just replaces the 4th column with the auto-acquire '.' and the 5th column with "OTHER" and also generates metadata)

import os
from mutperiodpy.Tkinter_scripts.TkinterDialog import Selections, TkinterDialog
from mutperiodpy.input_parsing.ParseCustomBed import parseCustomBed
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (getIsolatedParentDir, generateFilePath, getDataDirectory,
                                                                  DataTypeStr, generateMetadata, InputFormat, checkDirs,
                                                                  getAcceptableChromosomes)


def parseStandardBed(standardBedFilePaths, genomeFilePath, nucPosFilePath):

    customBedOutputFilePaths = list() # The list of file paths to be passed to the custom bed parser.

    # Parse the given files into custom bed format.
    for standardBedFilePath in standardBedFilePaths:

        print("\nWorking in:",os.path.basename(standardBedFilePath))
        if not os.path.basename(standardBedFilePath).endswith(".bed"):
            raise ValueError("Error:  Expected bed file format.")

        # Store useful paths and names.
        localRootDirectory = os.path.dirname(standardBedFilePath)
        intermediateFilesDir = os.path.join(localRootDirectory,"intermediate_files")
        checkDirs(intermediateFilesDir)
        dataGroupName = getIsolatedParentDir(standardBedFilePath)

        # Generate the output file path and metadata
        customBedOutputFilePath = generateFilePath(directory = intermediateFilesDir, dataGroup = dataGroupName,
                                                   dataType = DataTypeStr.customInput, fileExtension = ".bed")
        customBedOutputFilePaths.append(customBedOutputFilePath)
        generateMetadata(dataGroupName, getIsolatedParentDir(genomeFilePath), getIsolatedParentDir(nucPosFilePath), 
                         os.path.basename(standardBedFilePath), InputFormat.standardBed, localRootDirectory)

        # Get the list of acceptable chromosomes.
        acceptableChromosomes = getAcceptableChromosomes(genomeFilePath)

        # Iterate through the standard bed file entries preparing them for custom-bed input.
        print("Converting entries for custom bed input...")
        with open(standardBedFilePath, 'r') as standardBedFile:
            with open(customBedOutputFilePath, 'w') as customBedOutputFile:

                for line in standardBedFile:
                    
                    choppedUpLine = line.strip().split("\t")

                    # Make sure the lesion is in a valid chromosome.  Otherwise, skip it.
                    if not choppedUpLine[0] in acceptableChromosomes: continue

                    choppedUpLine[3] = '.'
                    choppedUpLine[4] = "OTHER"

                    customBedOutputFile.write('\t'.join(choppedUpLine[:6]) + '\n')


    # Pass the generated files to the custom bed parser.
    parseCustomBed(customBedOutputFilePaths, genomeFilePath, nucPosFilePath, False, False, False, False)




if __name__ == "__main__":

    # Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Standard Bed Data:",0,"dipy.bed",("Bed Files",".bed"),additionalFileEndings=("TA.bed",))    
    dialog.createFileSelector("Genome Fasta File:",1,("Fasta Files",".fa"))
    dialog.createFileSelector("Strongly Positioned Nucleosome File:",2,("Bed Files",".bed"))

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    standardBedFilePaths = list(selections.getFilePathGroups())[0]
    genomeFilePath = list(selections.getIndividualFilePaths())[0]
    nucPosFilePath = list(selections.getIndividualFilePaths())[1]

    parseStandardBed(standardBedFilePaths, genomeFilePath, nucPosFilePath)