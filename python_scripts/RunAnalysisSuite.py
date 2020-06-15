# This script runs a suite of scripts from this project to take a singlenuc context (or trinuc if it's already available) bed file
# and produce the normalized dyad position counts, along with all the relevant intermediate files.

from TkinterDialog import TkinterDialog, Selections
from typing import List
import os
from ConvertToTrinucContext import convertToTrinucContext
from GenerateMutationBackground import generateMutationBackground
from GenerateNucleosomeMutationBackground import generateNucleosomeMutationBackground
from CountNucleosomePositionMutations import countNucleosomePositionMutations
from NormalizeMutationCounts import normalizeMutationCounts

# Create the Tkinter dialog.
dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(__file__),"..","data"))
dialog.createMultipleFileSelector("Singlenuc/Trinuc Bed Files:",0,"singlenuc_context.bed",("Bed Files",".bed"),
                                  additionalFileEndings=("trinuc_context.bed",))
dialog.createFileSelector("Genome Fasta File:",1,("Fasta Files",".fa"))
dialog.createFileSelector("Strongly Positioned Nucleosome File:",2,("Bed Files",".bed"))
dialog.createReturnButton(3,0,2)
dialog.createQuitButton(3,2,2)

# Run the UI
dialog.mainloop()

# If no input was received (i.e. the UI was terminated prematurely), then quit!
if dialog.selections is None: quit()

# Get the user's input from the dialog.
selections: Selections = dialog.selections
mutationFilePaths: List[str] = list(selections.getFilePathGroups())[0] # A list of paths to bed mutation files
genomeFastaFilePath = list(selections.getIndividualFilePaths())[0] # The path to the human genome fasta file
strongPosNucleosomeFilePath = list(selections.getIndividualFilePaths())[1] # The path to the bed file of nucleosome positions

### Ensure that for every singlenuc file, there is a corresponding trinuc file in the same directory, creating one if necessary.

singlenucMutationFilePaths = list()
trinucMutationFilePaths = list()

# Sort the given mutation file paths by single/trinuc context.
for mutationFilePath in mutationFilePaths:
    if mutationFilePath.endswith("singlenuc_context.bed"): singlenucMutationFilePaths.append(mutationFilePath)
    elif mutationFilePath.endswith("trinuc_context.bed"): trinucMutationFilePaths.append(mutationFilePath)
    else: raise ValueError("Error: mutation file " + os.path.basename(mutationFilePath) + " does not conform to naming standards.\n\
                            The file should end in either \"trinuc_context.bed\" or \"singlenuc_context.bed\".")

# Convert to trinuc context where necessary.
print("\nConverting singlenuc mutation files to trinuc context...\n")
for singlenucMutationFilePath in singlenucMutationFilePaths:
    correspondingTrinucPath = singlenucMutationFilePath.rsplit("singlenuc_context.bed",1)[0] + "trinuc_context.bed"
    if not correspondingTrinucPath in trinucMutationFilePaths: 
        trinucMutationFilePaths += convertToTrinucContext((singlenucMutationFilePath,),genomeFastaFilePath)

### Run the rest of the analysis.

print("\nGenerating genome-wide mutation background...\n")
mutationBackgroundFilePaths = generateMutationBackground(trinucMutationFilePaths,genomeFastaFilePath)

print("\nGenerating nucleosome mutation background...\n")
nucleosomeMutationBackgroundFilePaths = generateNucleosomeMutationBackground(mutationBackgroundFilePaths, genomeFastaFilePath,
                                                                             strongPosNucleosomeFilePath)

print("\nCounting mutations at each dyad position...\n")                                                                             
nucleosomeMutationCountsFilePaths = countNucleosomePositionMutations(trinucMutationFilePaths, strongPosNucleosomeFilePath, False)

print("\nNormalizing counts with nucleosome background data...\n")
normalizedNucleosomeMutationCountsFilePaths = normalizeMutationCounts(nucleosomeMutationCountsFilePaths)