# This script takes the data obtained from mapping lesions cleaved by UVDE
# and converts it to a format suitable for downstream analysis.
# This is done by taking the 2 bp lesion and splitting it into 2 single base lesions.
from TkinterDialog import Selections, TkinterDialog
from UsefulFileSystemFunctions import getIsolatedParentDir, generateFilePath, dataTypes, generateMetadata
from UsefulBioinformaticsFunctions import baseChromosomes
import os, subprocess


def parseUVDESeq(UVDESeqFilePaths, genomeFilePath, nucPosFilePath):

    # Parse the given reads into singlenuc context.
    for UVDESeqFilePath in UVDESeqFilePaths:

        print("\nWorking in:",os.path.basename(UVDESeqFilePath))
        if not os.path.basename(UVDESeqFilePath).endswith(".bed"):
            raise ValueError("Error:  Expected bed file format.")

        # Store useful paths and names.
        localRootDirectory = os.path.dirname(UVDESeqFilePath)
        dataGroupName = getIsolatedParentDir(UVDESeqFilePath)

        # Generate the output file path and metadata
        singlenucOutputFilePath = generateFilePath(directory = localRootDirectory, dataGroup = dataGroupName,
                                                   context = "singlenuc", dataType = dataTypes.mutations,
                                                   fileExtension = ".bed")
        generateMetadata(dataGroupName, getIsolatedParentDir(genomeFilePath), getIsolatedParentDir(nucPosFilePath),
                         os.path.basename(UVDESeqFilePath), localRootDirectory)

        # Iterate through the 2 bp lesions, adding 2 single base lesions to the singlenuc output file for each.
        print("Converting 2-bp lesions to 2 single base lesions...")
        with open(UVDESeqFilePath, 'r') as UVDESeqFile:
            with open(singlenucOutputFilePath, 'w') as singlenucOutputFile:

                for line in UVDESeqFile:
                    
                    choppedUpLine = line.strip().split("\t")

                    # Make sure the lesion is in a valid chromosome.  Otherwise, skip it.
                    if not choppedUpLine[0] in baseChromosomes: continue

                    # Extract the relevant data from the line.
                    chromosome = choppedUpLine[0]
                    startPos = int(choppedUpLine[1])
                    endPos = int(choppedUpLine[2]) - 1
                    plusOrMinus = choppedUpLine[5]

                    for i in range(2):

                        # Write the two single base lesions from the one 2 bp lesion.
                        singlenucOutputFile.write('\t'.join((chromosome,str(startPos+i),str(endPos+i),"None","None",plusOrMinus)) + '\n')

        # Sort the output.
        print("Sorting output data...")
        subprocess.run(" ".join(("sort","-k1,1","-k2,2n",singlenucOutputFilePath,"-o",singlenucOutputFilePath)), 
                       shell = True, check = True)


if __name__ == "__main__":

    # Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(__file__),"..","data"))
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