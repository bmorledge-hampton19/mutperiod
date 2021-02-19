# This script takes data from the Alexandrov paper and parses it into a format acceptable for the rest of the pipeline.

import os
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog, Selections
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (generateFilePath, generateMetadata, DataTypeStr, InputFormat,
                                                                  checkDirs, getIsolatedParentDir, dataDirectory, getAcceptableChromosomes)
from mutperiodpy.input_parsing.ParseCustomBed import parseCustomBed


def parseAlexandrov (bedInputFilePaths, genomeFilePath, nucPosFilePath):

    outputBedFilePaths = list()

    for bedInputFilePath in bedInputFilePaths:

        print("\nWorking in:",os.path.basename(bedInputFilePath))

        # Get some important file system paths for the rest of the function and generate metadata.
        dataDirectory = os.path.dirname(bedInputFilePath)
        generateMetadata(os.path.basename(dataDirectory), getIsolatedParentDir(genomeFilePath), getIsolatedParentDir(nucPosFilePath), 
                            os.path.basename(bedInputFilePath), InputFormat.customBed, os.path.dirname(bedInputFilePath))

        intermediateFilesDir = os.path.join(dataDirectory,"intermediate_files")
        checkDirs(intermediateFilesDir)

        # Get the list of acceptable chromosomes
        acceptableChromosomes = getAcceptableChromosomes(genomeFilePath)

        # Generate the output file.
        outputBedFilePath = generateFilePath(directory = intermediateFilesDir, dataGroup = getIsolatedParentDir(bedInputFilePath),
                                             dataType = DataTypeStr.customInput, fileExtension = ".bed")

        # Write data to the output file.
        with open(bedInputFilePath, 'r') as bedInputFile:
            with open(outputBedFilePath, 'w') as outputBedFile:

                for line in bedInputFile:

                    choppedUpLine = str(line).strip().split('\t')

                    # Make sure we have a valid chromosome
                    if ("chr" + choppedUpLine[2]) in acceptableChromosomes and not '/' in choppedUpLine[5]:
                    
                        # Convert the line to custom bed format.
                        if choppedUpLine[5] == '-': choppedUpLine[5] = '*'
                        if choppedUpLine[6] == '-': choppedUpLine[6] = '*'
                        outputBedFile.write('\t'.join(("chr" + choppedUpLine[2], str(int(choppedUpLine[3])-1), choppedUpLine[4], 
                                                    choppedUpLine[5], choppedUpLine[6], '.', choppedUpLine[0])) + '\n')

        # Add the output file to the list.
        outputBedFilePaths.append(outputBedFilePath)

    # Pass the data to the custome bed parser.
    print("\nPassing data to custom bed parser.\n")
    parseCustomBed(outputBedFilePaths, genomeFilePath, nucPosFilePath, False, False, False)


if __name__ == "__main__":

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=dataDirectory)
    dialog.createMultipleFileSelector("Input Files:",0,"alexandrov.txt",("text files",".txt"))
    dialog.createFileSelector("Genome Fasta File:",1,("Fasta Files",".fa"))
    dialog.createFileSelector("Strongly Positioned Nucleosome File:",2,("Bed Files",".bed"))

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    bedInputFilePaths = list(selections.getFilePathGroups())[0]
    genomeFilePath = list(selections.getIndividualFilePaths())[0]
    nucPosFilePath = list(selections.getIndividualFilePaths())[1]

    parseAlexandrov(bedInputFilePaths, genomeFilePath, nucPosFilePath)