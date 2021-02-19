# This script takes raw nucleosome mutation count files, and passes them an R script which normalizes the data.

import os, subprocess, datetime
from typing import List
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog, Selections
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (getLinkerOffset, getContext, getDataDirectory, Metadata, 
                                                                  generateFilePath, DataTypeStr, rScriptsDirectory, checkForNucGroup)


# Pairs each background file path with its respective raw counts file path.
# Returns these pairings as a dictionary.
def getBackgroundRawPairs(backgroundCountsFilePaths):

    # Match each background file path to its respective raw counts file path.
    backgroundRawPairs = dict()
    for backgroundCountsFilePath in backgroundCountsFilePaths:

        if not DataTypeStr.nucMutBackground in os.path.basename(backgroundCountsFilePath): 
            raise ValueError("Background counts file should have \"" + DataTypeStr.nucMutBackground + "\" in the name.")

        # Generate the expected raw counts file path
        metadata = Metadata(backgroundCountsFilePath)
        rawCountsFilePath = generateFilePath(directory = metadata.directory, dataGroup = metadata.dataGroupName,
                                            linkerOffset = getLinkerOffset(backgroundCountsFilePath),
                                            usesNucGroup = checkForNucGroup(backgroundCountsFilePath),
                                            dataType = DataTypeStr.rawNucCounts, fileExtension = ".tsv")

        # Make sure it exists
        if not os.path.exists(rawCountsFilePath):
            raise ValueError("No raw counts file found to pair with " + backgroundCountsFilePath +
                             "\nExpected file with path: " + rawCountsFilePath)

        backgroundRawPairs[backgroundCountsFilePath] = rawCountsFilePath

    return backgroundRawPairs


# Attempts to pair each raw counts file in the custom raw directory to a raw counts file in the custom background directory
def getCustomBackgroundRawPairs(customRawCountsFilePaths, customBackgroundCountsDir):

    customBackgroundRawPairs = dict()
    backgroundMetadata = Metadata(customBackgroundCountsDir)

    # For every raw counts file in the customRawCountsDir, try to match it to a raw counts file in the customBackgroundCountsDir.
    for customRawCountsFilePath in customRawCountsFilePaths:
        customBackgroundCountsFilePath = generateFilePath(
            directory = backgroundMetadata.directory, dataGroup = backgroundMetadata.dataGroupName,
            linkerOffset = getLinkerOffset(customRawCountsFilePath), 
            usesNucGroup = checkForNucGroup(customRawCountsFilePath),
            dataType = DataTypeStr.rawNucCounts, fileExtension = ".tsv")
        assert os.path.exists(customBackgroundCountsFilePath), (
            "No counts file found to use as custom background for " + customRawCountsFilePath +
            "\nExpected file at: " + customBackgroundCountsFilePath)
        customBackgroundRawPairs[customBackgroundCountsFilePath] = customRawCountsFilePath

    return customBackgroundRawPairs


def normalizeCounts(backgroundCountsFilePaths: List[str], customRawCountsFilePaths: List[str] = list(), customBackgroundCountsDir = None):

    normalizedCountsFilePaths = list()

    backgroundRawPairs = getBackgroundRawPairs(backgroundCountsFilePaths)

    # Get the background-raw pairs from the custom directories, if they were given.
    if customBackgroundCountsDir is not None:
        customBackgroundRawPairs = getCustomBackgroundRawPairs(customRawCountsFilePaths, customBackgroundCountsDir)
        for customBackgroundCountsFilePath in customBackgroundRawPairs:
            backgroundRawPairs[customBackgroundCountsFilePath] = customBackgroundRawPairs[customBackgroundCountsFilePath]

    # Iterate through each background + raw counts pair
    for backgroundCountsFilePath in backgroundRawPairs:

        rawCountsFilePath = backgroundRawPairs[backgroundCountsFilePath]

        print("\nWorking with",os.path.basename(rawCountsFilePath),"and",os.path.basename(backgroundCountsFilePath))

        metadata = Metadata(rawCountsFilePath)

        # Generate the path to the normalized file.
        if DataTypeStr.rawNucCounts in backgroundCountsFilePath: context = "custom_context"
        else: context = getContext(backgroundCountsFilePath)
        normalizedCountsFilePath = generateFilePath(directory = metadata.directory, dataGroup = metadata.dataGroupName,
                                                    context = context, linkerOffset = getLinkerOffset(backgroundCountsFilePath),
                                                    usesNucGroup = checkForNucGroup(backgroundCountsFilePath),
                                                    dataType = DataTypeStr.normNucCounts, fileExtension = ".tsv")

        # Pass the file paths to the R script to generate the normalized counts file.
        print("Calling R script to generate normalized counts...")
        subprocess.run(" ".join(("Rscript",os.path.join(rScriptsDirectory,"NormalizeNucleosomeMutationCounts.R"),
                                 rawCountsFilePath,backgroundCountsFilePath,normalizedCountsFilePath)), 
                       shell = True, check = True)

        normalizedCountsFilePaths.append(normalizedCountsFilePath)

    # Document where the custom background counts came from in each relevant directory.
    if customBackgroundCountsDir is not None:
        for customRawCountsDir in set([os.path.dirname(customRawCountsFilePath) for customRawCountsFilePath in customRawCountsFilePaths]):
            metadata = Metadata(customRawCountsDir)
            customBackgroundInfoFilePath = generateFilePath(directory = metadata.directory, dataGroup = metadata.dataGroupName,
                                                            dataType = DataTypeStr.customBackgroundInfo, fileExtension = ".txt")
            with open(customBackgroundInfoFilePath, 'w') as customBackgroundInfoFile:
                customBackgroundInfoFile.write("Custom background directory: " + customBackgroundCountsDir + '\n')
                customBackgroundInfoFile.write("Last date used: " + str(datetime.datetime.now()).rsplit(':',1)[0] + '\n')

    return normalizedCountsFilePaths


def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Background Nucleosome Mutation Counts Files:",0,
                                      DataTypeStr.nucMutBackground + ".tsv",("Tab Seperated Values Files",".tsv"))

    customBackgroundSelector = dialog.createDynamicSelector(1, 0)
    customBackgroundSelector.initCheckboxController("Use another data set as custom background.")
    customBackgroundFileSelector = customBackgroundSelector.initDisplay(True, "customBackground")
    customBackgroundFileSelector.createFileSelector("Data Directory for nucleosome counts to be used as raw", 0, directory = True)
    customBackgroundFileSelector.createFileSelector("Data Directory for nucleosome counts to be used as background", 1, directory = True)
    customBackgroundSelector.initDisplayState()

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    backgroundCountsFilePaths = selections.getFilePathGroups()[0] # A list of background mutation counts file paths

    if customBackgroundSelector.getControllerVar():
        customRawCountsDir = selections.getFilePaths("customBackground")[0]
        customBackgroundCountsDir = selections.getFilePaths("customBackground")[1]
    else:
        customRawCountsDir = None
        customBackgroundCountsDir = None

    normalizeCounts(backgroundCountsFilePaths, customRawCountsDir, customBackgroundCountsDir)

if __name__ == "__main__": main()