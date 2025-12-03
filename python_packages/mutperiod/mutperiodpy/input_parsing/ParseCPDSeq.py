# This script takes CPD-seq data and converts it to custom bed format.
# Basically, it just replaces the 4th column with '.' and the 5th column with "OTHER" and also generates metadata.
# If a sixth column is present, it is kept as is, and the script assumes that it contains a strand designation.
# Otherwise, the strand is set to '+' in all rows. If the format is set to "split positions", each CPD is split
# into two single-base positions, otherwise the dinucleotide positions are maintained and will be coerced
# to half-base positions in ParseCustomBed.

import os
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import Selections, TkinterDialog
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (getIsolatedParentDir, generateFilePath, getDataDirectory,
                                                                  DataTypeStr, generateMetadata, InputFormat, checkDirs,
                                                                  getAcceptableChromosomes)
from benbiohelpers.CustomErrors import *

HALF_BASE_POSITION = 1
SPLIT_POSITIONS = 2

def parseCPDSeq(cpdSeqBedFilePaths: List[str], genomeFilePath, format: int):

    # This needs to be here to avoid a circular reference.
    from mutperiodpy.input_parsing.ParseCustomBed import parseCustomBed

    customBedOutputFilePaths = list() # The list of file paths to be passed to the custom bed parser.

    # Parse the given files into custom bed format.
    for cpdSeqBedFilePath in cpdSeqBedFilePaths:

        print("\nWorking in:",os.path.basename(cpdSeqBedFilePath))
        if not os.path.basename(cpdSeqBedFilePath).endswith(".bed"):
            raise InvalidPathError(cpdSeqBedFilePath, 
                                   "Given file does not appear to be in bed format. (missing \".bed\" extension)")

        # Store useful paths and names.
        localRootDirectory = os.path.dirname(cpdSeqBedFilePath)
        intermediateFilesDir = os.path.join(localRootDirectory,"intermediate_files")
        checkDirs(intermediateFilesDir)
        dataGroupName = getIsolatedParentDir(cpdSeqBedFilePath)

        # Generate the output file path and metadata
        customBedOutputFilePath = generateFilePath(directory = intermediateFilesDir, dataGroup = dataGroupName,
                                                   dataType = DataTypeStr.customInput, fileExtension = ".bed")
        customBedOutputFilePaths.append(customBedOutputFilePath)
        generateMetadata(dataGroupName, getIsolatedParentDir(genomeFilePath), 
                         os.path.basename(cpdSeqBedFilePath), InputFormat.standardBed, localRootDirectory)

        # Get the list of acceptable chromosomes.
        acceptableChromosomes = getAcceptableChromosomes(genomeFilePath)

        # Iterate through the standard bed file entries preparing them for custom-bed input.
        print("Converting entries for custom bed input...")
        with open(cpdSeqBedFilePath, 'r') as standardBedFile:
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

                    if format == HALF_BASE_POSITION:
                        customBedOutputFile.write('\t'.join(choppedUpLine[:6]) + '\n')
                    elif format == SPLIT_POSITIONS:
                        choppedUpLine[2] = str(int(choppedUpLine[2])-1)
                        customBedOutputFile.write('\t'.join(choppedUpLine[:6]) + '\n')
                        choppedUpLine[1] = str(int(choppedUpLine[1])+1)
                        choppedUpLine[2] = str(int(choppedUpLine[2])+1)
                        customBedOutputFile.write('\t'.join(choppedUpLine[:6]) + '\n')


    # Pass the generated files to the custom bed parser.
    return parseCustomBed(customBedOutputFilePaths, genomeFilePath)


if __name__ == "__main__":

    # Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory(), title = "Parse CPD-Seq")
    dialog.createMultipleFileSelector("CPDseq Bed Data:",0,"dipy.bed",("Bed Files",".bed"),additionalFileEndings=("TA.bed",))    
    dialog.createFileSelector("Genome Fasta File:",1,("Fasta Files",".fa"))
    dialog.createDropdown("Format:", 2, 0, ["Half-base position", "Split positions"])

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    cpdSeqBedFilePaths = list(selections.getFilePathGroups())[0]
    genomeFilePath = list(selections.getIndividualFilePaths())[0]
    format = HALF_BASE_POSITION if list(selections.getDropdownSelections())[0] == "Half-base position" else SPLIT_POSITIONS

    parseCPDSeq(cpdSeqBedFilePaths, genomeFilePath, format)
