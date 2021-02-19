# This script takes normalized nucleosome mutation counts files and passes them to an R script
# which outputs relevant data about them such as periodicity snr, assymetry, and differences between MSI and MSS data.

import os, subprocess, sys
from typing import List
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (DataTypeStr, getDataDirectory, Metadata,
                                                                  rScriptsDirectory, getContext, checkForNucGroup,
                                                                  getFilesInDirectory)
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog, Selections


# Given a list of file paths pointing to nucleosome mutation data, returns the paths that fit the given specifications
# normalizationMethods: A list of integers designating which normalization methods are allowed.
#   0: Raw, 1: Singlenuc, 3: trinuc, 5: pentanuc -1: custom
# singleNuc and nucGroup: Boolean values specifying whether the relevant nucleosome radii are allowed.
# MSS and MSI: Boolean values specifying whether the relevant microsatellite staibility status is allowed.
# acceptableCohorts: a list of cohorts.  The mutation group must belong to at least one of these cohorts to be accepted.
def getFilePathGroup(potentialFilePaths, normalizationMethods: List[int], singleNuc, nucGroup, acceptableCohorts: List[str]):
    
    filePathGroup = list() # The file paths to be returned.

    for potentialFilePath in potentialFilePaths:

        potentialFileName = os.path.basename(potentialFilePath)

        # Does it satisfy the normalization methods qualification? 
        # (Also ensure that we have nucleosome counts, whether raw or normalized.)
        if len(normalizationMethods) != 0:
            if DataTypeStr.rawNucCounts in potentialFileName and 0 in normalizationMethods:
                passed = True
            elif DataTypeStr.normNucCounts in potentialFileName and getContext(potentialFilePath, True) in normalizationMethods:
                passed = True
            else: continue

        # Does it satisfy the nucleosome radius qualification?
        if singleNuc or nucGroup:
            if checkForNucGroup(potentialFilePath) and nucGroup:
                passed = True
            elif not checkForNucGroup(potentialFilePath) and singleNuc:
                passed = True
            else: continue

        # Does it belong to one of the acceptable cohorts?
        if len(acceptableCohorts) != 0:
            filePathCohortDesignations = Metadata(potentialFilePath).cohorts
            acceptableCohortFound = False
            for cohort in filePathCohortDesignations:
                if cohort in acceptableCohorts:
                    acceptableCohortFound = True
                    continue
            if not acceptableCohortFound: continue

        # If we've made it this far, add the file path to the return group!
        filePathGroup.append(potentialFilePath)        
        
    return filePathGroup


def runNucleosomeMutationAnalysis(nucleosomeMutationCountsFilePaths: List[str], outputFilePath: str, 
                                  filePathGroup1: List[str] = list(), filePathGroup2: List[str] = list()):

    # Check for valid input.
    assert (len(filePathGroup1) == 0) == (len(filePathGroup2) == 0), (
        "One file path group contains file paths, but the other is empty.")
    assert len(nucleosomeMutationCountsFilePaths) > 0, (
        "No normalized counts files given.")
    assert outputFilePath.endswith(".rda") or outputFilePath.endswith(".tsv"), (
        "Output file should end with \".rda\" or \".tsv\".")

    # Write the inputs to a temporary file to be read by the R script
    inputsFilePath = os.path.join(os.getenv("HOME"), ".mutperiod","R_inputs.txt")

    with open(inputsFilePath, 'w') as inputsFile:
        if (len(filePathGroup1) == 0 and len(filePathGroup2) == 0):
            print("Generating inputs to run analysis without grouped comparison...")
            inputsFile.write('\n'.join(('$'.join(nucleosomeMutationCountsFilePaths), outputFilePath)) + '\n')
            
        else:
            print("Generating inputs to run analysis with grouped comparison...")
            inputsFile.write('\n'.join(('$'.join(nucleosomeMutationCountsFilePaths), outputFilePath,
                                        '$'.join(filePathGroup1), '$'.join(filePathGroup2))) + '\n')

    # Call the R script
    print("Calling R script...")
    subprocess.run(" ".join(("Rscript",os.path.join(rScriptsDirectory,"RunNucleosomeMutationAnalysis.R"),inputsFilePath)),
                   shell = True, check = True)

    print("Results can be found at",outputFilePath)


def parseArgs(args):
    
    # If only the subcommand was given, run the UI.
    if len(sys.argv) == 2: 
        main(); return

    # Determine what files were passed to each argument.
    filePathGroups = list()
    for i in range(3): filePathGroups.append(list())

    # Pulls all the relevant file paths out of the three groups that could have been passed as arguments.
    for i,filePaths in enumerate((args.nucleosomeMutationFilePaths, args.group_1, args.group_2)):
        if filePaths is not None:

            for filePath in filePaths:

                if os.path.isdir(filePath):
                    filePathGroups[i] += getFilesInDirectory(filePath, DataTypeStr.generalNucCounts + ".tsv")
                else: filePathGroups[i].append(filePath)

    # Make sure that any file paths passed to group 1 or group 2 are present in the default group.
    for i in range(1,3):
        for filePath in filePathGroups[i]:
            if filePath not in filePathGroups[0]: filePathGroups[0].append(filePath)

    runNucleosomeMutationAnalysis(filePathGroups[0], args.output_file_path,
                                  filePathGroups[1], filePathGroups[2])


def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory(), scrollable=True)
    dialog.createMultipleFileSelector("Nucleosome Mutation Counts files:",0,
                                      DataTypeStr.normNucCounts + ".tsv",("Tab Seperated Values Files",".tsv"),
                                      additionalFileEndings = (DataTypeStr.rawNucCounts + ".tsv",))
    dialog.createFileSelector("Output File", 1, ("R Data File", ".rda"), ("Tab Separated Values File", ".tsv"), newFile = True)

    mainGroupSearchRefine = dialog.createDynamicSelector(2, 0)
    mainGroupSearchRefine.initCheckboxController("Filter counts files")
    mainGroupSearchRefine.initDisplay(True).createNucMutGroupSubDialog("MainGroup", 0)
    mainGroupSearchRefine.initDisplayState()

    periodicityComparison = dialog.createDynamicSelector(3, 0)
    periodicityComparison.initCheckboxController("Compare periodicities between two groups")
    periodicityGroups = periodicityComparison.initDisplay(True,"periodicityGroups")

    # Create two "sub-dialogs" for each of the groups, allowing the user to specify the make-up of that group.
    for i, dialogID in enumerate(("Group1", "Group2")):
        periodicityGroups.createNucMutGroupSubDialog(dialogID, i+1)

    periodicityComparison.initDisplayState()

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    nucleosomeMutationCountsFilePaths = list(selections.getFilePathGroups())[0]
    outputFilePath = list(selections.getIndividualFilePaths())[0]
    
    # If group comparisons were requested, get the respective groups.
    filePathGroups = list()
    for i in range(3): filePathGroups.append(list())

    groups = ["MainGroup"]
    if periodicityComparison.getControllerVar(): groups += ("Group1", "Group2")

    for i, dialogID in enumerate(groups):

        # Was filtering even requested for the main group?
        if i == 0 and not mainGroupSearchRefine.getControllerVar():
            filePathGroups[i] = nucleosomeMutationCountsFilePaths
            continue

        # Determine what normalization methods were requested
        normalizationMethods = list()
        normalizationSelections = selections.getToggleStates(dialogID)[:5]
        if normalizationSelections[0]: normalizationMethods.append(1)
        if normalizationSelections[1]: normalizationMethods.append(3)
        if normalizationSelections[2]: normalizationMethods.append(5)
        if normalizationSelections[3]: normalizationMethods.append(-1)
        if normalizationSelections[4]: normalizationMethods.append(0)

        # Determine what microsatellite stability states were requested.
        acceptableCohorts = dict()
        if selections.getToggleStates(dialogID)[7]:
            MSSelection = selections.getDropdownSelections(dialogID+"MS")[0]
            if MSSelection == "MSS": acceptableCohorts["MSS"] = None
            elif MSSelection == "MSI": acceptableCohorts["MSI"] = None
            else: acceptableCohorts["MSS"] = None; acceptableCohorts["MSI"] = None

        # Check for custom cohort input
        if selections.getToggleStates(dialogID)[8]:

            customCohortsFilePath = selections.getIndividualFilePaths(dialogID + "CustomCohorts")[0]
            with open(customCohortsFilePath, 'r') as customCohortsFile:

                for line in customCohortsFile: acceptableCohorts[line.strip()] = None

        # Get the file paths associated with the given parameters.
        filePathGroups[i] += getFilePathGroup(nucleosomeMutationCountsFilePaths, normalizationMethods, 
                                              selections.getToggleStates(dialogID)[5], selections.getToggleStates(dialogID)[6],
                                              acceptableCohorts)

        #If this is the first pass through the loop, set the file paths list to the newly filtered list.
        if i == 0: nucleosomeMutationCountsFilePaths = filePathGroups[0]

    runNucleosomeMutationAnalysis(filePathGroups[0], outputFilePath,
                                  filePathGroups[1], filePathGroups[2])

if __name__ == "__main__": main()