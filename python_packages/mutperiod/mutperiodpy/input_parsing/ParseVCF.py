# This script takes VCF files and parses them into a format acceptable for the rest of the pipeline.

import os
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog, Selections
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (generateFilePath, generateMetadata, DataTypeStr, InputFormat,
                                                                  checkDirs, getIsolatedParentDir, getDataDirectory, 
                                                                  getAcceptableChromosomes)
from mutperiodpy.input_parsing.ParseCustomBed import parseCustomBed


def parseVCF (vcfInputFilePaths, genomeFilePath):

    outputBedFilePaths = list()

    for vcfInputFilePath in vcfInputFilePaths:

        print("\nWorking in:",os.path.basename(vcfInputFilePath))

        # Get some important file system paths for the rest of the function and generate metadata.
        dataDirectory = os.path.dirname(vcfInputFilePath)
        generateMetadata(os.path.basename(dataDirectory), getIsolatedParentDir(genomeFilePath), 
                         os.path.basename(vcfInputFilePath), InputFormat.vcf, os.path.dirname(vcfInputFilePath))

        intermediateFilesDir = os.path.join(dataDirectory,"intermediate_files")
        checkDirs(intermediateFilesDir)

        # Get the list of acceptable chromosomes
        acceptableChromosomes = getAcceptableChromosomes(genomeFilePath)

        # Generate the output file.
        outputBedFilePath = generateFilePath(directory = intermediateFilesDir, dataGroup = getIsolatedParentDir(vcfInputFilePath),
                                             dataType = DataTypeStr.customInput, fileExtension = ".bed")

        # Write data to the output file.
        with open(vcfInputFilePath, 'r') as vcfInputFile:
            with open(outputBedFilePath, 'w') as outputBedFile:

                for line in vcfInputFile:
                    
                    # Skip header lines.
                    if line.startswith("#"): continue

                    choppedUpLine = str(line).strip().split('\t')

                    # Make sure we have a valid chromosome.
                    if choppedUpLine[0] in acceptableChromosomes:
                    
                        # NOTE: Currently only knows how to handle SBSs. Will need to revise for more complicated VCF files.
                        if not (choppedUpLine[3].upper() in ('A','C','G','T') and choppedUpLine[4].upper() in ('A','C','G','T')):
                            print(f"Warning: line with non-SBS variant call encountered:\n{line}\nThe pipeline is not currently set up to handle this.")

                        # Convert the line to custom bed format.
                        outputBedFile.write('\t'.join((choppedUpLine[0], str(int(choppedUpLine[1])-1), choppedUpLine[1], 
                                                    choppedUpLine[3], choppedUpLine[4], '.')) + '\n')

        # Add the output file to the list.
        outputBedFilePaths.append(outputBedFilePath)

    # Pass the data to the custome bed parser.
    print("\nPassing data to custom bed parser.\n")
    return parseCustomBed(outputBedFilePaths, genomeFilePath, onlySingleBaseSubs = True)


if __name__ == "__main__":

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Input Files:",0,".vcf",("VCF files",".vcf"))
    dialog.createFileSelector("Genome Fasta File:",1,("Fasta Files",".fa"))

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    bedInputFilePaths = list(selections.getFilePathGroups())[0]
    genomeFilePath = list(selections.getIndividualFilePaths())[0]

    parseVCF(bedInputFilePaths, genomeFilePath)