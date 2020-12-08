# This script takes the data obtained from mapping lesions cleaved by UVDE
# and converts it to a format suitable for downstream analysis.
# This is done by taking the 2 bp lesion and splitting it into 2 single base lesions.

import os
from nucperiodpy.Tkinter_scripts.TkinterDialog import Selections, TkinterDialog
from nucperiodpy.input_parsing.ParseCustomBed import parseCustomBed
from nucperiodpy.helper_scripts.UsefulFileSystemFunctions import (getIsolatedParentDir, generateFilePath, getDataDirectory,
                                                                  DataTypeStr, generateMetadata, InputFormat, checkDirs,
                                                                  getAcceptableChromosomes)


def parseUVDESeq(UVDESeqFilePaths, genomeFilePath, nucPosFilePath):

    customBedOutputFilePaths = list() # The list of file paths to be passed to the custom bed parser.

    # Parse the given reads into singlenuc context.
    for UVDESeqFilePath in UVDESeqFilePaths:

        print("\nWorking in:",os.path.basename(UVDESeqFilePath))
        if not os.path.basename(UVDESeqFilePath).endswith(".bed"):
            raise ValueError("Error:  Expected bed file format.")

        # Store useful paths and names.
        localRootDirectory = os.path.dirname(UVDESeqFilePath)
        intermediateFilesDir = os.path.join(localRootDirectory,"intermediate_files")
        checkDirs(intermediateFilesDir)
        dataGroupName = getIsolatedParentDir(UVDESeqFilePath)

        # Generate the output file path and metadata
        customBedOutputFilePath = generateFilePath(directory = intermediateFilesDir, dataGroup = dataGroupName,
                                                   dataType = DataTypeStr.customInput, fileExtension = ".bed")
        customBedOutputFilePaths.append(customBedOutputFilePath)
        generateMetadata(dataGroupName, getIsolatedParentDir(genomeFilePath), getIsolatedParentDir(nucPosFilePath), 
                         os.path.basename(UVDESeqFilePath), InputFormat.UVDESeq, localRootDirectory)

        # Get the list of acceptable chromosomes.
        acceptableChromosomes = getAcceptableChromosomes(genomeFilePath)

        # Iterate through the 2 bp lesions, adding 2 single base lesions to the singlenuc output file for each.
        print("Converting 2-bp lesions to 2 single base lesions...")
        with open(UVDESeqFilePath, 'r') as UVDESeqFile:
            with open(customBedOutputFilePath, 'w') as customBedOutputFile:

                for line in UVDESeqFile:
                    
                    choppedUpLine = line.strip().split("\t")

                    # Make sure the lesion is in a valid chromosome.  Otherwise, skip it.
                    if not choppedUpLine[0] in acceptableChromosomes: continue

                    # Extract the relevant data from the line.
                    chromosome = choppedUpLine[0]
                    startPos = int(choppedUpLine[1])
                    endPos = int(choppedUpLine[2]) - 1
                    plusOrMinus = choppedUpLine[5]

                    for i in range(2):

                        # Write the two single base lesions from the one 2 bp lesion.
                        customBedOutputFile.write('\t'.join((chromosome,str(startPos+i),str(endPos+i),".","OTHER",plusOrMinus)) + '\n')


    # Pass the generated files to the custom bed parser.
    parseCustomBed(customBedOutputFilePaths, genomeFilePath, nucPosFilePath, False, False, False)




if __name__ == "__main__":

    # Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("UVDE-seq data:",0,"dipy.bed",("Bed Files",".bed"),additionalFileEndings=("TA.bed",))    
    dialog.createFileSelector("Genome Fasta File:",1,("Fasta Files",".fa"))
    dialog.createFileSelector("Strongly Positioned Nucleosome File:",2,("Bed Files",".bed"))
    dialog.createExitButtons(3,0)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    UVDESeqFilePaths = list(selections.getFilePathGroups())[0]
    genomeFilePath = list(selections.getIndividualFilePaths())[0]
    nucPosFilePath = list(selections.getIndividualFilePaths())[1]

    parseUVDESeq(UVDESeqFilePaths, genomeFilePath, nucPosFilePath)