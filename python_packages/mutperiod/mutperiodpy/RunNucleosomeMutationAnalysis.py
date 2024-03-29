# This script takes normalized nucleosome mutation counts files and passes them to an R script
# which outputs relevant data about them such as periodicity snr, assymetry, and differences between MSI and MSS data.

import os, subprocess, sys
from typing import List

from benbiohelpers.CustomErrors import UserInputError, InvalidPathError, checkIfPathExists
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (DataTypeStr, getDataDirectory, Metadata,
                                                                  rScriptsDirectory, getContext, checkForNucGroup, getExpectedPeriod)
from benbiohelpers.FileSystemHandling.DirectoryHandling import getFilesInDirectory
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog, Selections


# Given a list of file paths pointing to nucleosome mutation data, returns the paths that fit the given specifications
# normalizationMethods: A list of integers designating which normalization methods are allowed.
#   0: Raw, 1/2: Singlenuc or dinuc, 3/4: trinuc or quadrunuc, 5/6: pentanuc or hexanuc, -1: custom
# singleNuc and nucGroup: Boolean values specifying whether the relevant nucleosome radii are allowed.
# MSS and MSI: Boolean values specifying whether the relevant microsatellite staibility status is allowed.
# acceptableCohorts: a list of cohorts.  The mutation group must belong to at least one of these cohorts to be accepted.
def getFilePathGroup(potentialFilePaths, normalizationMethods: List[int], singleNuc, nucGroup,
                     acceptableMSCohorts: List[str], acceptableMutSigCohorts: List[str], acceptableCustomCohorts: List[str],
                     acceptableNucleosomeMaps: List[str]):

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

        # Does it belong to one of the acceptable cohorts in each category?
        invalidCohortGroup = False
        for acceptableCohortsGroup in (acceptableMSCohorts, acceptableMutSigCohorts, acceptableCustomCohorts):

            if len(acceptableCohortsGroup) != 0:

                filePathCohortDesignations = Metadata(potentialFilePath).cohorts
                acceptableCohortFound = False
                for cohort in filePathCohortDesignations:
                    if cohort in acceptableCohortsGroup:
                        acceptableCohortFound = True
                        continue

                if not acceptableCohortFound: 
                    invalidCohortGroup = True
                    continue

        if invalidCohortGroup: continue

        # Does it belong to one of the acceptable nucleosome maps given?
        if len(acceptableNucleosomeMaps) != 0:
            filePathNucleosomeMap = Metadata(potentialFilePath).nucPosName
            if not filePathNucleosomeMap in acceptableNucleosomeMaps: continue


        # If we've made it this far, add the file path to the return group!
        filePathGroup.append(potentialFilePath)        
        
    return filePathGroup


def runNucleosomeMutationAnalysis(nucleosomeMutationCountsFilePaths: List[str], outputFilePath: str, overridePeakPeriodicityWithExpected,
                                  alignStrands, filePathGroup1: List[str] = list(), filePathGroup2: List[str] = list()):

    # Check for valid input.
    if (len(filePathGroup1) == 0) != (len(filePathGroup2) == 0):
        raise UserInputError("One file path group contains file paths, but the other is empty.")
    if len(nucleosomeMutationCountsFilePaths) == 0:
        raise UserInputError("No nucleosome counts files given.")
    if not (outputFilePath.endswith(".rda") or outputFilePath.endswith(".tsv")):
        raise InvalidPathError(outputFilePath, "Given output file does not end with \".rda\" or \".tsv\":")
    try:
        outputFile = open(outputFilePath, 'w' )
        outputFile.close()
    except IOError:
        raise InvalidPathError(outputFilePath, "Given output file path is not writeable: ")

    # Retrieve the expected periods for each of the given counts files.
    expectedPeriods = [str(getExpectedPeriod(nucleosomeMutationCountsFilePath)) for nucleosomeMutationCountsFilePath in nucleosomeMutationCountsFilePaths]

    # Write the inputs to a temporary file to be read by the R script
    inputsFilePath = os.path.join(os.getenv("HOME"), ".mutperiod","R_inputs.txt")

    with open(inputsFilePath, 'w') as inputsFile:
        if (len(filePathGroup1) == 0 and len(filePathGroup2) == 0):
            print("Generating inputs to run analysis without grouped comparison...")
            inputsFile.write('\n'.join(('$'.join(nucleosomeMutationCountsFilePaths), outputFilePath, 
                                        str(overridePeakPeriodicityWithExpected), 
                                        '$'.join(expectedPeriods),
                                        str(alignStrands))) + '\n')
            
        else:
            print("Generating inputs to run analysis with grouped comparison...")
            inputsFile.write('\n'.join(('$'.join(nucleosomeMutationCountsFilePaths), outputFilePath,
                                        '$'.join(filePathGroup1), '$'.join(filePathGroup2), 
                                        str(overridePeakPeriodicityWithExpected),
                                        '$'.join(expectedPeriods),
                                        str(alignStrands))) + '\n')

    # Call the R script
    print("Calling R script...")
    subprocess.run(("Rscript",os.path.join(rScriptsDirectory,"RunNucleosomeMutationAnalysis.R"),inputsFilePath), check = True)

    print("Results can be found at",outputFilePath)


def parseArgs(args):
    
    # If only the subcommand was given, run the UI.
    if len(sys.argv) == 2: 
        main(); return

    # Make sure an output file path was given.
    if args.output_file_path is None: raise UserInputError("No output file path was given.")

    # Determine what files were passed to each argument.
    filePathGroups = list()
    for i in range(3): filePathGroups.append(list())

    # Pulls all the relevant file paths out of the three groups that could have been passed as arguments.
    for i,filePaths in enumerate((args.nucleosomeMutationFilePaths, args.group_1, args.group_2)):
        if filePaths is not None:

            for filePath in filePaths:

                checkIfPathExists(filePath)
                if os.path.isdir(filePath):
                    filePathGroups[i] += [os.path.abspath(filePath) for filePath in getFilesInDirectory(filePath, DataTypeStr.generalNucCounts + ".tsv")]
                else: filePathGroups[i].append(os.path.abspath(filePath))
            
        filePathGroups[i] = set(filePathGroups[i])

    # Make sure that any file paths passed to group 1 or group 2 are present in the default group.
    filePathGroups[0] = filePathGroups[0] | filePathGroups[1] | filePathGroups[2]

    runNucleosomeMutationAnalysis(list(filePathGroups[0]), args.output_file_path, args.use_expected_periodicity, args.align_strands,
                                  list(filePathGroups[1]), list(filePathGroups[2]))


def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory(), scrollable=True, title = "Nucleosome Mutation Analysis")
    dialog.createMultipleFileSelector("Nucleosome Mutation Counts files:",0,
                                      DataTypeStr.normNucCounts + ".tsv",("Tab Seperated Values Files",".tsv"),
                                      additionalFileEndings = (DataTypeStr.rawNucCounts + ".tsv",))
    dialog.createFileSelector("Output File", 1, ("R Data File", ".rda"), ("Tab Separated Values File", ".tsv"), newFile = True)
    
    dialog.createCheckbox("Use expected periodicity from nucleosome maps instead of peak periodicity", 2, 0)
    dialog.createCheckbox("Align both DNA strands to run 5' to 3' before running the analysis", 3, 0)
    dialog.createLabel('', 4, 0)

    mainGroupSearchRefine = dialog.createDynamicSelector(5, 0)
    mainGroupSearchRefine.initCheckboxController("Filter counts files")
    mainGroupSearchRefine.initDisplay(True).createNucMutGroupSubDialog("MainGroup", 0)
    mainGroupSearchRefine.initDisplayState()

    periodicityComparison = dialog.createDynamicSelector(6, 0)
    periodicityComparison.initCheckboxController("Compare periodicities between two groups")
    periodicityGroupType = periodicityComparison.initDisplay(True,"periodicityGroupType")

    periodicityGroupTypeSelector = periodicityGroupType.createDynamicSelector(0,0)
    periodicityGroupTypeSelector.initDropdownController("Compare periodicities...", ("Within original selection", "Against a newly selected group"))
    periodicityGroupsWithin = periodicityGroupTypeSelector.initDisplay("Within original selection", "withinGroup")
    secondaryPeriodicityGroup = periodicityGroupTypeSelector.initDisplay("Against a newly selected group", "secondaryGroupFilePaths")

    # Create two "sub-dialogs" for each of the groups, allowing the user to specify the make-up of that group for the "withinGroup" dialog.
    for i, dialogID in enumerate(("Sub-Group 1", "Sub-Group 2")):
        periodicityGroupsWithin.createNucMutGroupSubDialog(dialogID, i+1)

    # Create one multiple file selector and one sub-dialog for the "secondaryGroup" dialog
    secondaryPeriodicityGroup.createMultipleFileSelector("Nucleosome Mutation Counts files:",0,
                                                        DataTypeStr.normNucCounts + ".tsv",("Tab Seperated Values Files",".tsv"),
                                                        additionalFileEndings = (DataTypeStr.rawNucCounts + ".tsv",))

    secondaryPeriodicityGroupSearchRefine = secondaryPeriodicityGroup.createDynamicSelector(2, 0)
    secondaryPeriodicityGroupSearchRefine.initCheckboxController("Filter counts files")
    secondaryPeriodicityGroupSearchRefine.initDisplay(True).createNucMutGroupSubDialog("Secondary Group", 0)
    secondaryPeriodicityGroupSearchRefine.initDisplayState()

    periodicityGroupTypeSelector.initDisplayState()
    periodicityComparison.initDisplayState()

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    nucleosomeMutationCountsFilePaths = selections.getFilePathGroups()[0]
    if periodicityGroupTypeSelector.getControllerVar() == "Against a newly selected group":
        secondaryNucMutCountsFilePaths = selections.getFilePathGroups("secondaryGroupFilePaths")[0]
    outputFilePath = list(selections.getIndividualFilePaths())[0]

    # Get the default periodicity value, testing the string to see if it is a valid float, if necessary.
    overridePeakPeriodWithExpected = bool(selections.getToggleStates()[0])
    alignStrands = bool(selections.getToggleStates()[1])

    # If group comparisons were requested, get the respective groups.
    filePathGroups:List[list] = list()
    for i in range(3): filePathGroups.append(list())

    groups = ["MainGroup"]
    if periodicityComparison.getControllerVar(): 
        if periodicityGroupTypeSelector.getControllerVar() == "Within original selection": groups += ("Sub-Group 1", "Sub-Group 2")
        else: groups.append("Secondary Group")

    for i, dialogID in enumerate(groups):

        # Was filtering even requested for the main group?
        if i == 0 and not mainGroupSearchRefine.getControllerVar():
            filePathGroups[0] = nucleosomeMutationCountsFilePaths
            continue

        # If we are examining the secondary group, check to see if filtering was even requested.
        if dialogID == "Secondary Group" and not secondaryPeriodicityGroupSearchRefine.getControllerVar():
            assert i == 1, "Secondary group encountered on unexpected iteration of for loop: " + str(i)
            filePathGroups[2] = secondaryNucMutCountsFilePaths
            filePathGroups[1] = filePathGroups[0].copy()
            filePathGroups[0] += filePathGroups[2]
            continue

        # Determine what normalization methods were requested
        normalizationMethods = list()
        normalizationSelections = selections.getToggleStates(dialogID)[:5]
        if normalizationSelections[0]: normalizationMethods += (1,2)
        if normalizationSelections[1]: normalizationMethods += (3,4)
        if normalizationSelections[2]: normalizationMethods += (5,6)
        if normalizationSelections[3]: normalizationMethods.append(-1)
        if normalizationSelections[4]: normalizationMethods.append(0)

        # Determine what microsatellite stability states were requested.
        acceptableMSCohorts = dict()
        if selections.getToggleStates(dialogID)[7]:
            MSSelection = selections.getDropdownSelections(dialogID+"MS")[0]
            if MSSelection == "MSS": acceptableMSCohorts["MSS"] = None
            elif MSSelection == "MSI": acceptableMSCohorts["MSI"] = None
            else: acceptableMSCohorts["MSS"] = None; acceptableMSCohorts["MSI"] = None

        # Determine what mutation signature states were requested.
        acceptableMutSigCohorts = dict()
        if selections.getToggleStates(dialogID)[8]:
            with open(selections.getIndividualFilePaths(dialogID+"MutSig")[0], 'r') as mutSigsFile:
                for line in mutSigsFile: acceptableMutSigCohorts["mut_sig_" + line.strip()] = None

        # Check for custom cohort input
        acceptableCustomCohorts = dict()
        if selections.getToggleStates(dialogID)[9]:

            customCohortsFilePath = selections.getIndividualFilePaths(dialogID + "CustomCohorts")[0]
            with open(customCohortsFilePath, 'r') as customCohortsFile:

                for line in customCohortsFile: acceptableCustomCohorts[line.strip()] = None

        # Check for nucleosome map input
        acceptableNucleosomeMaps = dict()
        if selections.getToggleStates(dialogID)[10]:

            acceptableNucMapsFilePath = selections.getIndividualFilePaths(dialogID + "NucleosomeMaps")[0]
            with open(acceptableNucMapsFilePath, 'r') as acceptableNucMapsFile:

                for line in acceptableNucMapsFile: acceptableNucleosomeMaps[line.strip()] = None

        # Get the file paths associated with the given parameters.
        if dialogID != "Secondary Group":
            filePathGroups[i] += getFilePathGroup(nucleosomeMutationCountsFilePaths, normalizationMethods, 
                                                  selections.getToggleStates(dialogID)[5], selections.getToggleStates(dialogID)[6],
                                                  acceptableMSCohorts, acceptableMutSigCohorts, acceptableCustomCohorts,
                                                  acceptableNucleosomeMaps)
        else: 
            assert i == 1, "Secondary group encountered on unexpected iteration of for loop: " + str(i)
            filePathGroups[2] += getFilePathGroup(secondaryNucMutCountsFilePaths, normalizationMethods, 
                                                  selections.getToggleStates(dialogID)[5], selections.getToggleStates(dialogID)[6],
                                                  acceptableMSCohorts, acceptableMutSigCohorts, acceptableCustomCohorts,
                                                  acceptableNucleosomeMaps)
            filePathGroups[1] = filePathGroups[0].copy()
            filePathGroups[0] += filePathGroups[2]

        #If this is the first pass through the loop, set the file paths list to the newly filtered list.
        if i == 0: nucleosomeMutationCountsFilePaths = filePathGroups[0]

    runNucleosomeMutationAnalysis(filePathGroups[0], outputFilePath, overridePeakPeriodWithExpected, alignStrands,
                                  filePathGroups[1], filePathGroups[2])

if __name__ == "__main__": main()