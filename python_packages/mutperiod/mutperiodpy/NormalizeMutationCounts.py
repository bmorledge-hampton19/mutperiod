# This script takes raw nucleosome mutation count files, and passes them an R script which normalizes the data.

import os, subprocess, datetime
from typing import List, Dict
from benbiohelpers.CustomErrors import UserInputError, InvalidPathError
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog, Selections
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (getLinkerOffset, getContext, getDataDirectory, Metadata, 
                                                                  generateFilePath, DataTypeStr, rScriptsDirectory, checkForNucGroup)


# Pairs each background file path with its respective raw counts file path.
# Returns these pairings as a dictionary.
def getBackgroundRawPairs(backgroundCountsFilePaths):

    # Match each background file path to its respective raw counts file path.
    backgroundRawPairs: Dict[str, List[str]] = dict()
    for backgroundCountsFilePath in backgroundCountsFilePaths:

        if not DataTypeStr.nucMutBackground in os.path.basename(backgroundCountsFilePath): 
            raise InvalidPathError("Background counts file should have \"" + DataTypeStr.nucMutBackground + "\" in the name.  Given:")

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

        if backgroundCountsFilePath not in backgroundRawPairs:
            backgroundRawPairs[backgroundCountsFilePath] = list()
        backgroundRawPairs[backgroundCountsFilePath].append(rawCountsFilePath)

    return backgroundRawPairs


# Attempts to pair each raw counts file in the custom raw directory to a raw counts file in the custom background directory
def getCustomBackgroundRawPairs(customRawCountsFilePaths, customBackgroundCountsDir):

    customBackgroundRawPairs: Dict[str, List[str]] = dict()
    

    # For every raw counts file given, try to match it to a raw counts file in the customBackgroundCountsDir.
    for customRawCountsFilePath in customRawCountsFilePaths:

        rawMetadata = Metadata(customRawCountsFilePath)
        backgroundDir = os.path.join(customBackgroundCountsDir,rawMetadata.nucPosName)
        if not os.path.exists(backgroundDir):
            raise UserInputError ("Expected a directory at " + backgroundDir + " to contain the background for " +
                              customRawCountsFilePath + " but the directory does not exist.  Have you forgotten to run "
                              "the analysis for the related nucleosome map?")
        backgroundMetadata = Metadata(backgroundDir)

        customBackgroundCountsFilePath = generateFilePath(
            directory = backgroundMetadata.directory, dataGroup = backgroundMetadata.dataGroupName,
            linkerOffset = getLinkerOffset(customRawCountsFilePath), 
            usesNucGroup = checkForNucGroup(customRawCountsFilePath),
            dataType = DataTypeStr.rawNucCounts, fileExtension = ".tsv")
        if not os.path.exists(customBackgroundCountsFilePath):
            raise UserInputError("Expected file at " + customBackgroundCountsFilePath + " to use as custom background for "
                             + customRawCountsFilePath + " but this file does not exist.  Have you forgotten to "
                             "run the relevant analysis to generate it?")
        if customBackgroundCountsFilePath not in customBackgroundRawPairs:
            customBackgroundRawPairs[customBackgroundCountsFilePath] = list()
        customBackgroundRawPairs[customBackgroundCountsFilePath].append(customRawCountsFilePath)

    return customBackgroundRawPairs


# Given a file path, use its metadata to determine how many 
# features (mutations, repair reads, etc.) are present in its parent data set.
def getParentDataFeatureCounts(filePath):
    counts = Metadata(os.path.dirname(Metadata(filePath).parentDataFilePath)).mutationCounts
    assert counts is not None, "Feature counts were never recorded for the parent data from: " + filePath
    return(counts)


def normalizeCounts(backgroundCountsFilePaths: List[str], customRawCountsFilePaths: List[str] = list(), 
                    customBackgroundCountsDir = None, includeAlternativeScaling = False):

    normalizedCountsFilePaths = list()

    backgroundRawPairs = getBackgroundRawPairs(backgroundCountsFilePaths)

    # Get the background-raw pairs from the custom directories, if they were given.
    if customBackgroundCountsDir is not None:
        customBackgroundRawPairs = getCustomBackgroundRawPairs(customRawCountsFilePaths, customBackgroundCountsDir)
        for customBackgroundCountsFilePath in customBackgroundRawPairs:
            assert customBackgroundCountsFilePath not in backgroundRawPairs, "Unexpected intersection!"
            backgroundRawPairs[customBackgroundCountsFilePath] = customBackgroundRawPairs[customBackgroundCountsFilePath]

    # Iterate through each background + raw counts pair
    for backgroundCountsFilePath in backgroundRawPairs:
        for rawCountsFilePath in backgroundRawPairs[backgroundCountsFilePath]:

            print("\nWorking with",os.path.basename(rawCountsFilePath),"and",os.path.basename(backgroundCountsFilePath))

            metadata = Metadata(rawCountsFilePath)

            # Generate the path to the normalized file.
            if DataTypeStr.rawNucCounts in backgroundCountsFilePath: context = "custom_context"
            else: context = getContext(backgroundCountsFilePath)
            normalizedCountsFilePath = generateFilePath(directory = metadata.directory, dataGroup = metadata.dataGroupName,
                                                        context = context, linkerOffset = getLinkerOffset(backgroundCountsFilePath),
                                                        usesNucGroup = checkForNucGroup(backgroundCountsFilePath),
                                                        dataType = DataTypeStr.normNucCounts, fileExtension = ".tsv")

            # Prepare the arguments to the subprocess call.
            args = ["Rscript",os.path.join(rScriptsDirectory,"NormalizeNucleosomeMutationCounts.R"),
                    rawCountsFilePath,backgroundCountsFilePath,normalizedCountsFilePath]

            # If alternative scaling is requested, determine the appropriate scaling factor and add it to the arguments
            if includeAlternativeScaling:

                # If we are normalizing by sequence context, just revert the automatic scaling.
                if customBackgroundCountsDir is None: args.append(1)

                # If we are normalizing by a custom context, scale based on the relative sizes of the parent background and raw data sets.
                else:
                    args.append(str(getParentDataFeatureCounts(backgroundCountsFilePath) /
                                    getParentDataFeatureCounts(rawCountsFilePath)))                    

            # Pass the file paths to the R script to generate the normalized counts file.
            print("Calling R script to generate normalized counts...")
            subprocess.run(args, check = True)

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
    customBackgroundSelector.initCheckboxController("Include files with another data set as custom background.")
    customBackgroundFileSelector = customBackgroundSelector.initDisplay(True, "customBackground")
    customBackgroundFileSelector.createMultipleFileSelector("Raw nucleosome counts to be normalized", 0, 
                                                            DataTypeStr.rawNucCounts + ".tsv", ("TSV Files", ".tsv"))
    customBackgroundFileSelector.createFileSelector("Data Directory for nucleosome counts to be used as background", 1, directory = True)
    customBackgroundSelector.initDisplayState()

    dialog.createCheckbox("Include alternative scaling factor indepedent of nucleosome map.", 2, 0)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    backgroundCountsFilePaths = selections.getFilePathGroups()[0] # A list of background mutation counts file paths

    if customBackgroundSelector.getControllerVar():
        customRawCountsFilePaths = selections.getFilePathGroups("customBackground")[0]
        customBackgroundCountsDir = selections.getIndividualFilePaths("customBackground")[0]
    else:
        customRawCountsFilePaths = None
        customBackgroundCountsDir = None

    includeAlternativeScaling = selections.getToggleStates()[0]

    normalizeCounts(backgroundCountsFilePaths, customRawCountsFilePaths, customBackgroundCountsDir, includeAlternativeScaling)

if __name__ == "__main__": main()