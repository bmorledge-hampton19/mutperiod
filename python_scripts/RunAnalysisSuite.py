# This script runs a suite of scripts from this project to take a singlenuc context (or trinuc if it's already available) bed file
# and produce the normalized dyad position counts, along with all the relevant intermediate files.

from TkinterDialog import TkinterDialog, Selections
from typing import List
import os
from UsefulFileSystemFunctions import dataTypes
from ExpandContext import expandContext
from GenerateMutationBackground import generateMutationBackground
from GenerateNucleosomeMutationBackground import generateNucleosomeMutationBackground
from CountNucleosomePositionMutations import countNucleosomePositionMutations
from NormalizeMutationCounts import normalizeCounts

# Create the Tkinter dialog.
dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(__file__),"..","data"))
dialog.createMultipleFileSelector("Bed Mutation Files:",0,dataTypes.mutations + ".bed",("Bed Files",".bed"))
dialog.createDropdown("Background Context",1,0,("Trinuc","Singlenuc", "Pentanuc"))
dialog.createCheckbox("Include 30 bp linker DNA",1,2)
dialog.createExitButtons(2,0)

# Run the UI
dialog.mainloop()

# If no input was received (i.e. the UI was terminated prematurely), then quit!
if dialog.selections is None: quit()

# Get the user's input from the dialog.
selections: Selections = dialog.selections
mutationFilePaths: List[str] = list(selections.getFilePathGroups())[0] # A list of paths to bed mutation files
backgroundContext = list(selections.getDropdownSelections())[0] # The context to be used when generating the background
includeLinker = list(selections.getToggleStates())[0] # Whether or not to include 30 bp linker DNA in nucleosome dyad positions

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
nucleosomeMutationBackgroundFilePaths = generateNucleosomeMutationBackground(mutationBackgroundFilePaths, linkerOffset)

print("\nCounting mutations at each dyad position...\n")                                                                             
nucleosomeMutationCountsFilePaths = countNucleosomePositionMutations(updatedMutationFilePaths, linkerOffset)

print("\nNormalizing counts with nucleosome background data...\n")
normalizedNucleosomeMutationCountsFilePaths = normalizeCounts(nucleosomeMutationBackgroundFilePaths)