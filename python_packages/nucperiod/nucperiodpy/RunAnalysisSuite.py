# This script runs a suite of scripts from this project to take a singlenuc context (or trinuc if it's already available) bed file
# and produce the normalized dyad position counts, along with all the relevant intermediate files.

from typing import List
import os, sys
from nucperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog, Selections
from nucperiodpy.helper_scripts.UsefulFileSystemFunctions import DataTypeStr, getDataDirectory, getFilesInDirectory
from nucperiodpy.ExpandContext import expandContext
from nucperiodpy.GenerateMutationBackground import generateMutationBackground
from nucperiodpy.GenerateNucleosomeMutationBackground import generateNucleosomeMutationBackground
from nucperiodpy.CountNucleosomePositionMutations import countNucleosomePositionMutations
from nucperiodpy.NormalizeMutationCounts import normalizeCounts


def parseArgs(args):
    
    # If only the subcommand was given, run the UI.
    if len(sys.argv) == 2: 
        main(); return

    # Get the bed mutation files from the given paths, searching directories if necessary.
    finalBedMutationPaths = list()
    for mutationFilePath in args.mutationFilePaths:
        if os.path.isdir(mutationFilePath):
            finalBedMutationPaths += getFilesInDirectory(mutationFilePath, DataTypeStr.mutations + ".bed")
        else: finalBedMutationPaths.append(mutationFilePath)

    assert len(finalBedMutationPaths) > 0, "No bed mutation files were found."

    # Determine what normalization method was selected.
    normalizationMethod = "No Normalization"
    if args.context_normalization == 1: normalizationMethod = "Singlenuc"
    elif args.context_normalization == 3: normalizationMethod = "Trinuc"
    elif args.context_normalization == 5: normalizationMethod = "Pentanuc"
    elif args.background is not None: normalizationMethod = "Custom Background"

    runAnalysisSuite(finalBedMutationPaths, normalizationMethod, args.background, 
                     args.singlenuc_radius, args.add_linker, args.nuc_group_radius)


def main():

    # Create the Tkinter dialog.
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Bed Mutation Files:",0,DataTypeStr.mutations + ".bed",("Bed Files",".bed"))

    normalizationSelector = dialog.createDynamicSelector(1, 0)
    normalizationSelector.initDropdownController("Normalization Method",("Trinuc", "Singlenuc", "Pentanuc", "Custom Background", "No Normalization"))
    customBackgroundFileSelector = normalizationSelector.initDisplay("Custom Background", "customBackground")
    customBackgroundFileSelector.createFileSelector("Custom Background Directory:", 0, ("Bed Files", ".bed"), directory = True)
    normalizationSelector.initDisplayState()

    selectNucleosomeDyadRadius = dialog.createDynamicSelector(2,0)
    selectNucleosomeDyadRadius.initCheckboxController("Run analysis with a single nucleosome dyad radius (73 bp)")
    linkerSelectionDialog = selectNucleosomeDyadRadius.initDisplay(1, "singleNuc")
    linkerSelectionDialog.createCheckbox("Include 30 bp linker DNA on either side of single nucleosome dyad radius.",0,0)
    selectNucleosomeDyadRadius.initDisplayState()

    dialog.createCheckbox("Count with a nucleosome group radius (1000 bp)", 3, 0)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    mutationFilePaths: List[str] = selections.getFilePathGroups()[0] # A list of paths to bed mutation files
    normalizationMethod = normalizationSelector.getControllerVar() # The normalization method to be used.
    if normalizationMethod == "Custom Background":
        customBackgroundDir = selections.getFilePaths("customBackground")[0] # Where to find raw counts files to use as custom background
    else: customBackgroundDir = None
    useSingleNucRadius = selectNucleosomeDyadRadius.getControllerVar() # Whether or not to generate data with a 73 bp single nuc dyad radius
    if useSingleNucRadius: 
        includeLinker = selections.getToggleStates("singleNuc")[0] # Whether or not to include 30 bp linker DNA in nucleosome dyad positions
    else: includeLinker = False
    useNucGroupRadius = selections.getToggleStates()[0] # Whether or not to generate data with a 1000 bp nuc group dyad radius

    runAnalysisSuite(mutationFilePaths, normalizationMethod, customBackgroundDir, useSingleNucRadius, 
                     includeLinker, useNucGroupRadius)


def runAnalysisSuite(mutationFilePaths: List[str], normalizationMethod, customBackgroundDir, useSingleNucRadius,
                     includeLinker, useNucGroupRadius):

    # Make sure at least one radius was selected.
    if not useNucGroupRadius and not useSingleNucRadius:
        raise ValueError("Must select at least one radius.")

    # Make sure at least one input file was found.
    assert len(mutationFilePaths) > 0, "No valid input files given."

    # Convert background context to int
    if normalizationMethod == "Singlenuc":
        normalizationMethodNum = 1
    elif normalizationMethod == "Trinuc":
        normalizationMethodNum = 3
    elif normalizationMethod == "Pentanuc":
        normalizationMethodNum = 5
    elif normalizationMethod in ("No Normalization", "Custom Background"):
        normalizationMethodNum = None
    else: raise ValueError("Matching strings is hard.")

    # Set the linker offset
    if includeLinker: linkerOffset = 30
    else: linkerOffset = 0

    ### Ensure that every mutation file has a context sufficient for the requested background.

    # Returns the number associated with the context of the given mutation file.
    def determineMutationFileContext(mutationFilePath: str):

        if mutationFilePath.endswith("singlenuc_context_mutations.bed"): return 1
        elif mutationFilePath.endswith("trinuc_context_mutations.bed"): return 3
        elif mutationFilePath.endswith("pentanuc_context_mutations.bed"): return 5
        else: raise ValueError("Unexpected file ending for " + os.path.basename(mutationFilePath))

    # create a new list of mutation file paths, replacing any with contexts that are too low.
    if normalizationMethodNum is not None:
        print("\nExpanding file context where necessary...\n")
        updatedMutationFilePaths = list()
        for mutationFilePath in mutationFilePaths:
            if determineMutationFileContext(mutationFilePath) < normalizationMethodNum:
                updatedMutationFilePaths += expandContext((mutationFilePath,),normalizationMethodNum)
            else: updatedMutationFilePaths.append(mutationFilePath)
    else: updatedMutationFilePaths = mutationFilePaths

    ### Run the rest of the analysis.

    print("\nCounting mutations at each dyad position...\n")                                                                             
    nucleosomeMutationCountsFilePaths = countNucleosomePositionMutations(updatedMutationFilePaths, useSingleNucRadius,
                                                                         useNucGroupRadius, linkerOffset)

    if normalizationMethodNum is not None:

        print("\nGenerating genome-wide mutation background...\n")
        mutationBackgroundFilePaths = generateMutationBackground(updatedMutationFilePaths,normalizationMethodNum)

        print("\nGenerating nucleosome mutation background...\n")
        nucleosomeMutationBackgroundFilePaths = generateNucleosomeMutationBackground(mutationBackgroundFilePaths, useSingleNucRadius, 
                                                                                     useNucGroupRadius, linkerOffset)

        print("\nNormalizing counts with nucleosome background data...\n")
        normalizeCounts(nucleosomeMutationBackgroundFilePaths)

    elif normalizationMethod == "Custom Background":
        print("\nNormalizing counts using custom background data...\n")
        normalizeCounts(list(), nucleosomeMutationCountsFilePaths, customBackgroundDir)

if __name__ == "__main__": main()