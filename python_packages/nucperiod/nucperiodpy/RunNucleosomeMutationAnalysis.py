# This script takes normalized nucleosome mutation counts files and passes them to an R script
# which outputs relevant data about them such as periodicity snr, assymetry, and differences between MSI and MSS data.

import os, subprocess
from typing import List
from nucperiodpy.helper_scripts.UsefulFileSystemFunctions import (DataTypeStr, getDataDirectory, Metadata,
                                                                  rScriptsDirectory, getContext, checkForNucGroup)
from nucperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog, Selections


# Given a list of file paths pointing to nucleosome mutation data, returns the paths that fit the given specifications
# normalizationMethods: A list of integers designating which normalization methods are allowed.
#   0: Raw, 1: Singlenuc, 3: trinuc, 5: pentanuc
# singleNuc and nucGroup: Boolean values specifying whether the relevant nucleosome radii are allowed.
# MSS and MSI: Boolean values specifying whether the relevant microsatellite staibility status is allowed.
def getFilePathGroup(potentialFilePaths, normalizationMethods: List[int], singleNuc, nucGroup, MSS, MSI):
    
    filePathGroup = list() # The file paths to be returned.

    for potentialFilePath in potentialFilePaths:

        potentialFileName = os.path.basename(potentialFilePath)

        # Does it satisfy the normalization methods qualification? 
        # (Also ensure that we have nucleosome counts, whether raw or normalized.)
        if DataTypeStr.rawNucCounts in potentialFileName and 0 in normalizationMethods:
            passed = True
        elif DataTypeStr.normNucCounts in potentialFileName and getContext(potentialFilePath, True) in normalizationMethods:
            passed = True
        else: continue

        # Does it satisfy the nucleosome radius qualification?
        if checkForNucGroup(potentialFilePath) and nucGroup:
            passed = True
        elif not checkForNucGroup(potentialFilePath) and singleNuc:
            passed = True
        else: continue

        # Does it satisfy the microsatellite stability qualifications?
        cohortDesignations = Metadata(potentialFilePath).cohorts
        if "MSI" in cohortDesignations and not MSI: continue
        elif "MSS" in cohortDesignations and not MSS: continue

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
    assert outputFilePath.endswith(".rda"), (
        "Output file should end with \".rda\".")

    # Write the inputs to a temporary file to be read by the R script
    inputsFilePath = os.path.join(rScriptsDirectory,"inputs.txt")

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


def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Normalized Nucleosome Mutation Counts files:",0,
                                      DataTypeStr.normNucCounts + ".tsv",("Tab Seperated Values Files",".tsv"))
    dialog.createFileSelector("Output File", 1, ("R Data File", ".rda"), newFile = True)

    dialog.createNucMutGroupSubDialog("MainGroup", 2)

    periodicityComparison = dialog.createDynamicSelector(3, 0)
    periodicityComparison.initCheckboxController("Compare periodicities between two groups")
    periodicityGroups = periodicityComparison.initDisplay(True,"periodicityGroups")
    periodicityGroups.createLabel("",0,0)

    # Create two "sub-dialogs" for each of the groups, allowing the user to specify the make-up of that group.
    for i, dialogID in enumerate(("Group1", "Group2")):
        periodicityGroups.createNucMutGroupSubDialog(dialogID, i+1)

    periodicityComparison.initDisplayState()

    dialog.createExitButtons(4,0)

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
    filePathGroups.append(list())
    filePathGroups.append(list())
    filePathGroups.append(list())

    groups = ["MainGroup"]
    if periodicityComparison.getControllerVar(): groups += ("Group1", "Group2")

    for i, dialogID in enumerate(groups):

        # Determine what normalization methods were requested
        normalizationMethods = list()
        normalizationSelections = selections.getToggleStates(dialogID)[:5]
        if normalizationSelections[0]: normalizationMethods.append(0)
        if normalizationSelections[1]: normalizationMethods.append(1)
        if normalizationSelections[2]: normalizationMethods.append(3)
        if normalizationSelections[3]: normalizationMethods.append(5)
        if normalizationSelections[4]: normalizationMethods.append(-1)

        # Ensure valid input was given
        assert len(normalizationMethods) > 0, (
            "No normalization method chosen for " + dialogID)
        assert selections.getToggleStates(dialogID)[5] or selections.getToggleStates(dialogID)[6], (
            "No nucleosome radius given for " + dialogID)

        # Determine what microsatellite stability states were requested.
        MSSelection = selections.getDropdownSelections(dialogID)[0]
        MSS = True
        MSI = True
        if MSSelection == "MSS": MSI = False
        elif MSSelection == "MSI": MSS = False

        # Get the file paths associated with the given parameters.
        filePathGroups[i] += getFilePathGroup(nucleosomeMutationCountsFilePaths, normalizationMethods, 
                                                selections.getToggleStates(dialogID)[5], selections.getToggleStates(dialogID)[6],
                                                MSS, MSI)

    runNucleosomeMutationAnalysis(filePathGroups[0], outputFilePath,
                                  filePathGroups[1], filePathGroups[2])

if __name__ == "__main__": main()