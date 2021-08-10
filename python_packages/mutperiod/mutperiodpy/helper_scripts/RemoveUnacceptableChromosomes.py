# This script takes a bed formatted file and a genome fasta file and removes any chromosomes
# from the bed file that aren't in the fasta file.

import os
from typing import List
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getAcceptableChromosomes, getDataDirectory
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog


def removeUnacceptableChromosomes(bedFilePaths: List[str], genomeFastaFilepath):

    acceptableChromosomes = getAcceptableChromosomes(genomeFastaFilepath)

    for bedFilePath in bedFilePaths:

        print("\nWorking in",os.path.basename(bedFilePath))

        # Create an intermediate file to write acceptable entries to.
        acceptableChromBedFilePath = bedFilePath.rsplit('.', 1)[0] + "_acceptable_only.bed"

        removedEntries = 0
        unacceptableChromosomes = list()

        # Read through the bed file, writing only entries with acceptable chromosomes to the output file.
        with open(bedFilePath, 'r') as bedFile:
            with open(acceptableChromBedFilePath, 'w') as acceptableChromBedFile:

                for line in bedFile:
                    chromosome = line.split()[0]
                    if chromosome in acceptableChromosomes:
                        acceptableChromBedFile.write(line)
                    else:
                        removedEntries += 1
                        if chromosome not in unacceptableChromosomes:
                            unacceptableChromosomes.append(chromosome)
                            print("Unallowed chromosome found:",chromosome)

        # Rewrite the original file.
        print("Removed",removedEntries,"Entries.  Rewriting original file...")
        os.replace(acceptableChromBedFilePath, bedFilePath)


def main():

    # Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Bed Files:",0,"I_have_bad_chromosomes.bed",("Bed Files",".bed"))    
    dialog.createFileSelector("Genome Fasta File:",1,("Fasta Files",".fa"))

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    removeUnacceptableChromosomes(dialog.selections.getFilePathGroups()[0], dialog.selections.getIndividualFilePaths()[0])


if __name__ == "__main__": main()