# This script takes raw nucleosome mutation count files, and passes them an R script which normalizes the data.

import os, subprocess
from typing import List
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog, Selections
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (getLinkerOffset, getContext, dataDirectory, Metadata, 
                                                                  generateFilePath, dataTypes, RPackagesDirectory, checkForNucGroup)


# Pairs each background file path with its respective raw counts file path.
# Returns these pairings as a dictionary.
def getBackgroundRawPairs(backgroundCountsFilePaths):

    # Match each background file path to its respective raw counts file path.
    backgroundRawPairs = dict()
    for backgroundCountsFilePath in backgroundCountsFilePaths:

        if not dataTypes.nucMutBackground in os.path.basename(backgroundCountsFilePath): 
            raise ValueError("Background counts file should have \"" + dataTypes.nucMutBackground + "\" in the name.")

        # Generate the expected raw counts file path
        metadata = Metadata(backgroundCountsFilePath)
        rawCountsFilePath = generateFilePath(directory = metadata.directory, dataGroup = metadata.dataGroupName,
                                            linkerOffset = getLinkerOffset(backgroundCountsFilePath),
                                            usesNucGroup = checkForNucGroup(backgroundCountsFilePath),
                                            dataType = dataTypes.rawNucCounts, fileExtension = ".tsv")

        # Make sure it exists
        if not os.path.exists(rawCountsFilePath):
            raise ValueError("No raw counts file found to pair with " + backgroundCountsFilePath +
                             "\nExpected file with path: " + rawCountsFilePath)

        backgroundRawPairs[backgroundCountsFilePath] = rawCountsFilePath

    return backgroundRawPairs



def normalizeCounts(backgroundCountsFilePaths: List[str]):

    normalizedCountsFilePaths = list()

    backgroundRawPairs = getBackgroundRawPairs(backgroundCountsFilePaths)

    # Iterate through each background + raw counts pair
    for backgroundCountsFilePath in backgroundRawPairs:

        rawCountsFilePath = backgroundRawPairs[backgroundCountsFilePath]

        print("\nWorking with",os.path.basename(rawCountsFilePath),"and",os.path.basename(backgroundCountsFilePath))

        metadata = Metadata(backgroundCountsFilePath)

        # Generate the path to the normalized file.
        normalizedCountsFilePath = generateFilePath(directory = metadata.directory, dataGroup = metadata.dataGroupName,
                                                    context = getContext(backgroundCountsFilePath),
                                                    linkerOffset = getLinkerOffset(backgroundCountsFilePath),
                                                    usesNucGroup = checkForNucGroup(backgroundCountsFilePath),
                                                    dataType = dataTypes.normNucCounts, fileExtension = ".tsv")

        # Pass the path to the file paths to the R script to generate the normalized counts file.
        print("Calling R script to generate normalized counts...")
        subprocess.run(" ".join(("Rscript",os.path.join(RPackagesDirectory,"RunNucPeriod","NormalizeNucleosomeMutationCounts.R"),
                                 rawCountsFilePath,backgroundCountsFilePath,normalizedCountsFilePath)), 
                       shell = True, check = True)

        normalizedCountsFilePaths.append(normalizedCountsFilePath)

    return normalizedCountsFilePaths


if __name__ == "__main__":

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=dataDirectory)
    dialog.createMultipleFileSelector("Background Nucleosome Mutation Counts Files:",0,
                                      dataTypes.nucMutBackground + ".tsv",("Tab Seperated Values Files",".tsv"))
    dialog.createExitButtons(1,0)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    backgroundCountsFilePaths = list(selections.getFilePathGroups())[0] # A list of background mutation counts file paths

    normalizeCounts(backgroundCountsFilePaths)