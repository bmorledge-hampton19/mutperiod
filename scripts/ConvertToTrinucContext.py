# This script will, when given a bed formatted file with single base mutation entries and a corresponding fasta genome file,
# create a bed file with the trinucleotide context around the mutation.
# The input files will be selected using a tkinter interface for ease of use.

from TkinterDialog import TkinterDialog, Selections
from UsefulBioinformaticsFunctions import bedToFasta
import os, subprocess

# Expands the range of each mutation position in the original mutation file to encompass one extra base on either side.
def expandBedToTrinucRegion(singleBaseBedFilePath,trinucExpansionFilePath):
    "Expands the range of each mutation position in the original mutation file to encompass one extra base on either side."

    with open(trinucExpansionFilePath,'w') as trinucExpansionFile:
        with open(singleBaseBedFilePath, 'r') as singleBaseBedFile:

            print("Writing expanded mutation indicies to intermediate bed file...")
            for line in singleBaseBedFile:

                # Get a list of all the arguments for a single mutation in the bed file.
                choppedUpLine = line.strip().split('\t')

                # Make sure the line is in a single-nucleotide context.
                if not int(choppedUpLine[1]) - int(choppedUpLine[2]) == -1:
                    raise ValueError(line + " does not give a single nucleotide location to expand.")

                # Expand the position of the mutation to be one extra base on either side.
                choppedUpLine[1] = str(int(choppedUpLine[1]) - 1)
                choppedUpLine[2] = str(int(choppedUpLine[2]) + 1)

                # Write the results to the trinucExpansion file
                trinucExpansionFile.write("\t".join(choppedUpLine)+"\n")


# Uses the trinuc reads fasta file to create a new bed file with the trinuc mutational context.
def generateTrinucContext(singleBaseBedFilePath,trinucReadsFilePath,trinucContextFilePath):
    "Uses the trinuc reads fasta file to create a new bed file with the trinuc mutational context."

    print("Using fasta file to write trinuc context to new bed file...")
    # Open the singlenuc context bed file and the trinuc fasta reads that will be combined to create the trinuc context.
    with open(singleBaseBedFilePath, 'r') as singleBaseBedFile:
        with open(trinucReadsFilePath, 'r') as trinucReadsFile:
            with open(trinucContextFilePath, 'w') as trinucContextFile:

                # Work through the singlenuc context bed file one mutation at a time.
                for line in singleBaseBedFile:

                    # Split each line on tab characters.
                    choppedUpLine = line.strip().split("\t")

                    # Make sure that the corresponding line in the fasta file actually has the right nucleotide coordinates.
                    if not trinucReadsFile.readline().strip().find(choppedUpLine[3]) > 0:
                        raise ValueError(choppedUpLine[3] + " not found in expected position in corresponding fasta file")

                    # Replace the mutation's identifier with the trinuc context.
                    choppedUpLine[3] = trinucReadsFile.readline().strip()

                    # Write the result to the new trinuc context file.
                    trinucContextFile.write("\t".join(choppedUpLine)+"\n")


# Create the Tkinter dialog.
dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(__file__),"..","data"))
dialog.createMultipleFileSelector("Single-Base Bed File:",0,("Bed Files",".bed"))
dialog.createFileSelector("Human Genome Fasta File:",1,("Fasta Files",".fa"))
dialog.createReturnButton(2,0,2)
dialog.createQuitButton(2,2,2)

# Run the UI
dialog.mainloop()

# If no input was received (i.e. the UI was terminated prematurely), then quit!
if dialog.selections is None: quit()

# Get the user's input from the dialog.
selections: Selections = dialog.selections
singleBaseBedFilePaths = list(selections.getFilePathGroups())[0] # A list of paths to original bed mutation files
humanGenomeFastaFilePath = list(selections.getIndividualFilePaths())[0] # The path to the human genome fasta file

for singleBaseBedFilePath in singleBaseBedFilePaths:

    print("\nWorking in:",os.path.split(singleBaseBedFilePath)[1])
    if not "_singlenuc_context" in os.path.split(singleBaseBedFilePath)[1]:
        raise ValueError("Error:  Expected file with \"_singlenuc_context\" in the name.")

    # Get some information on the file system and generate file paths for convenience.
    workingDirectory = os.path.dirname(singleBaseBedFilePath) # The working directory for the current data group
    dataGroupName = os.path.split(singleBaseBedFilePath)[1].split("_singlenuc_context")[0] # The name of mutation data group
    trinucExpansionFilePath = os.path.join(workingDirectory,"intermediate_files",dataGroupName+"_trinuc_expansion.bed")
        # The path to an intermediate file fed to bedtools to generate the trinucleotide sequence.
    trinucReadsFilePath = os.path.join(workingDirectory,"intermediate_files",dataGroupName+"_trinuc_reads.fa")
        # The path to an intermediate fasta file that contains the trinuc reads generated from the locations in the trinucExpansion file.
    trinucContextFilePath = os.path.join(workingDirectory,dataGroupName+"_trinuc_context.bed") # The final output file.

    # Create a directory for intermediate files if it does not already exist...
    if not os.path.exists(os.path.join(workingDirectory,"intermediate_files")):
        os.mkdir(os.path.join(workingDirectory,"intermediate_files"))

    # Expand the nucleotide coordinates in the singlenuc context bed file to encompass a trinuc span.
    expandBedToTrinucRegion(singleBaseBedFilePath,trinucExpansionFilePath)

    # Convert the trinuc coordinates in the bed file to the referenced nucleotides in fasta format.
    bedToFasta(trinucExpansionFilePath,humanGenomeFastaFilePath,trinucReadsFilePath)

    # Using the newly generated fasta file, create a new bed file with the trinucleotide context.
    generateTrinucContext(singleBaseBedFilePath,trinucReadsFilePath,trinucContextFilePath)