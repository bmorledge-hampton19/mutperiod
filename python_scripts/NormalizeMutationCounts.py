# This script takes raw nucleosome mutation count files, and passes them an R script which normalizes the data.
import os, subprocess
from TkinterDialog import TkinterDialog, Selections
from typing import List

def normalizeMutationCounts(rawNucleosomeMutationCountsFilePaths: List[str]):

    normalizedNucleosomeMutationCountsFilePaths = list()

    for rawNucleosomeMutationCountsFilePath in rawNucleosomeMutationCountsFilePaths:

        print("\nWorking with",os.path.split(rawNucleosomeMutationCountsFilePath)[1])

        # Pass the path to the raw counts file to the R script to generate the normalized counts file.
        print("Calling R script to generate normalized counts...")
        subprocess.run(" ".join(("Rscript",os.path.join(os.path.dirname(__file__),"..","R_scripts","RunNucPeriod",
                                                        "NormalizeNucleosomeMutationCounts.R"),
                                 rawNucleosomeMutationCountsFilePath)), shell = True, check = True)

        # Generate the path to the normalized file.
        normalizedNucleosomeMutationCountsFilePath = rawNucleosomeMutationCountsFilePath.rsplit("mutation_counts.tsv")[0]
        normalizedNucleosomeMutationCountsFilePath += "mutation_counts_normalized.tsv"
        if not os.path.exists(normalizedNucleosomeMutationCountsFilePath):
            raise ValueError("Normalized counts file should have been generated at " + normalizedNucleosomeMutationCountsFilePath + 
                             " but the file does not exist.")
        normalizedNucleosomeMutationCountsFilePaths.append(normalizedNucleosomeMutationCountsFilePath)

    return normalizedNucleosomeMutationCountsFilePaths


if __name__ == "__main__":

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(__file__),"..","data"))
    dialog.createMultipleFileSelector("Raw Nucleosome Mutation Counts files:",0,
                                      "nucleosome_mutation_counts.tsv",("Tab Seperated Values Files",".tsv"))
    dialog.createReturnButton(1,0,2)
    dialog.createQuitButton(1,2,2)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    rawMutationCountsFilePaths = list(selections.getFilePathGroups())[0] # A list of mutation file paths

    normalizeMutationCounts(rawMutationCountsFilePaths)