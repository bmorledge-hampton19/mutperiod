# This script takes raw nucleosome mutation count files, and passes them an R script which normalizes the data.
import os, subprocess
from TkinterDialog import TkinterDialog, Selections
from typing import List
from UsefulFileSystemFunctions import getLinkerDNAAmount, getContext

# Pairs each background file path with its respective raw counts file path.
# Returns these pairings as a dictionary.
def getBackgroundRawPairs(rawCountsFilePaths, backgroundCountsFilePaths):

    # Create a dictionary matching each mutation group to its file path.
    rawCountsMutationGroups = dict()
    for rawCountsFilePath in rawCountsFilePaths:

        if not "raw_nucleosome_mutation_counts" in os.path.basename(rawCountsFilePath): 
            raise ValueError("Raw counts file should have \"raw_nucleosome_mutation_counts\" in the name.")

        dataGroupName = os.path.basename(rawCountsFilePath).split("_raw_nucleosome_mutation_counts")[0]

        if dataGroupName in rawCountsMutationGroups:
            raise ValueError("The mutation group name " + dataGroupName + " appears in multiple raw " +
                             "counts files but is supposed to be unique.")
        
        rawCountsMutationGroups[dataGroupName] = rawCountsFilePath

    # Match each background file path to its respective raw counts file path.
    backgroundRawPairs = dict()
    for backgroundCountsFilePath in backgroundCountsFilePaths:

        if not "nucleosome_mutation_background" in os.path.basename(backgroundCountsFilePath): 
            raise ValueError("Background counts file should have \"nucleosome_mutation_background\" in the name.")

        # Split the file name at the context identifier to determine the mutation group name.
        backgroundContext = getContext(backgroundCountsFilePath)
        dataGroupName = os.path.basename(backgroundCountsFilePath).split('_'+backgroundContext)[0]

        if not dataGroupName in rawCountsMutationGroups:
            raise ValueError("A background counts file was given for mutation group " + dataGroupName +
                             ", but no corresponding raw counts file was given.")

        backgroundRawPairs[backgroundCountsFilePath] = rawCountsMutationGroups[dataGroupName]

    return backgroundRawPairs



def normalizeCounts(rawCountsFilePaths: List[str], backgroundCountsFilePaths: List[str]):

    normalizedCountsFilePaths = list()

    backgroundRawPairs = getBackgroundRawPairs(rawCountsFilePaths, backgroundCountsFilePaths)

    for backgroundCountsFilePath in backgroundRawPairs:

        rawCountsFilePath = backgroundRawPairs[backgroundCountsFilePath]

        print("\nWorking with",os.path.basename(rawCountsFilePath),"and",os.path.basename(backgroundCountsFilePath))

        # Generate the path to the normalized file.
        extraLinker = getLinkerDNAAmount(backgroundCountsFilePath)
        normalizedCountsFilePath = rawCountsFilePath.rsplit("raw_nucleosome_mutation_counts",1)[0]
        normalizedCountsFilePath += getContext(backgroundCountsFilePath) + '_'
        if extraLinker > 0: normalizedCountsFilePath += str(extraLinker) + "linker+_"
        normalizedCountsFilePath += "normalized_nucleosome_mutation_counts.tsv"

        # Pass the path to the file paths to the R script to generate the normalized counts file.
        print("Calling R script to generate normalized counts...")
        subprocess.run(" ".join(("Rscript",os.path.join(os.path.dirname(__file__),"..","R_scripts","RunNucPeriod",
                                                        "NormalizeNucleosomeMutationCounts.R"),
                                 rawCountsFilePath,backgroundCountsFilePath,normalizedCountsFilePath)), 
                       shell = True, check = True)

        normalizedCountsFilePaths.append(normalizedCountsFilePath)

    return normalizedCountsFilePaths


if __name__ == "__main__":

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(__file__),"..","data"))
    dialog.createMultipleFileSelector("Raw Nucleosome Mutation Counts Files:",0,
                                      "raw_nucleosome_mutation_counts.tsv",("Tab Seperated Values Files",".tsv"))
    dialog.createMultipleFileSelector("Background NucleosomeMutation Counts Files:",1,
                                      "nucleosome_mutation_background.tsv",("Tab Seperated Values Files",".tsv"))
    dialog.createReturnButton(2,0,2)
    dialog.createQuitButton(2,2,2)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    rawCountsFilePaths = list(selections.getFilePathGroups())[0] # A list of raw mutation counts file paths
    backgroundCountsFilePaths = list(selections.getFilePathGroups())[1] # A list of background mutation counts file paths

    normalizeCounts(rawCountsFilePaths, backgroundCountsFilePaths)