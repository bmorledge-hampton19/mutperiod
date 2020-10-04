# This script runs a suite of scripts from this project to take a singlenuc context (or trinuc if it's already available) bed file
# and produce the normalized dyad position counts, along with all the relevant intermediate files.

from typing import List
import os
from nucperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog, Selections
from nucperiodpy.helper_scripts.UsefulFileSystemFunctions import DataTypeStr, dataDirectory
from nucperiodpy.ExpandContext import expandContext
from nucperiodpy.GenerateMutationBackground import generateMutationBackground
from nucperiodpy.GenerateNucleosomeMutationBackground import generateNucleosomeMutationBackground
from nucperiodpy.CountNucleosomePositionMutations import countNucleosomePositionMutations
from nucperiodpy.NormalizeMutationCounts import normalizeCounts


def main():

    # Create the Tkinter dialog.
    dialog = TkinterDialog(workingDirectory=dataDirectory)
    dialog.createMultipleFileSelector("Bed Mutation Files:",0,DataTypeStr.mutations + ".bed",("Bed Files",".bed"))
    dialog.createDropdown("Background Context",1,0,("Trinuc","Singlenuc", "Pentanuc"))

    selectNucleosomeDyadRadius = dialog.createDynamicSelector(2,0)
    selectNucleosomeDyadRadius.initCheckboxController("Run analysis with a single nucleosome dyad radius (73 bp)")
    linkerSelectionDialog = selectNucleosomeDyadRadius.initDisplay(1, "singleNuc")
    linkerSelectionDialog.createCheckbox("Include 30 bp linker DNA on either side of single nucleosome dyad radius.",0,0)
    selectNucleosomeDyadRadius.initDisplay(0)
    selectNucleosomeDyadRadius.initDisplayState()

    dialog.createCheckbox("Count with a nucleosome group radius (1000 bp)", 3, 0)

    dialog.createExitButtons(4,0)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    mutationFilePaths: List[str] = list(selections.getFilePathGroups())[0] # A list of paths to bed mutation files
    backgroundContext = list(selections.getDropdownSelections())[0] # The context to be used when generating the background
    useSingleNucRadius = selectNucleosomeDyadRadius.getControllerVar() # Whether or not to generate data with a 73 bp single nuc dyad radius
    if useSingleNucRadius: 
        includeLinker = list(selections.getToggleStates("singleNuc"))[0] # Whether or not to include 30 bp linker DNA in nucleosome dyad positions
    useNucGroupRadius = selections.getToggleStates()[0] # Whether or not to generate data with a 1000 bp nuc group dyad radius

    # Make sure at least one radius was selected.
    if not useNucGroupRadius and not useSingleNucRadius:
        raise ValueError("Must select at least one radius.")

    # Convert background context to int
    if backgroundContext == "Singlenuc":
        backgroundContextNum = 1
    elif backgroundContext == "Trinuc":
        backgroundContextNum = 3
    elif backgroundContext == "Pentanuc":
        backgroundContextNum = 5
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
    print("\nExpanding file context where necessary...\n")
    updatedMutationFilePaths = list()
    for mutationFilePath in mutationFilePaths:
        if determineMutationFileContext(mutationFilePath) < backgroundContextNum:
            updatedMutationFilePaths += expandContext((mutationFilePath,),backgroundContextNum)
        else: updatedMutationFilePaths.append(mutationFilePath)

    ### Run the rest of the analysis.

    print("\nGenerating genome-wide mutation background...\n")
    mutationBackgroundFilePaths = generateMutationBackground(updatedMutationFilePaths,backgroundContextNum)

    print("\nGenerating nucleosome mutation background...\n")
    nucleosomeMutationBackgroundFilePaths = generateNucleosomeMutationBackground(mutationBackgroundFilePaths, useSingleNucRadius, 
                                                                                useNucGroupRadius, linkerOffset)

    print("\nCounting mutations at each dyad position...\n")                                                                             
    nucleosomeMutationCountsFilePaths = countNucleosomePositionMutations(updatedMutationFilePaths, useSingleNucRadius,
                                                                        useNucGroupRadius, linkerOffset)

    print("\nNormalizing counts with nucleosome background data...\n")
    normalizedNucleosomeMutationCountsFilePaths = normalizeCounts(nucleosomeMutationBackgroundFilePaths)

if __name__ == "__main__": main()