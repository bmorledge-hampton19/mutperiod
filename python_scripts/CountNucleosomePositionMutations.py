# This script takes a mutation file and the coordinates for DNA bases around nucleosomes and calculates how many
# mutations occured at each dyad position for each (and both) strands.
# NOTE:  Both input files must be sorted for this script to run properly. 
#        (Sorted first by chromosome (string) and then by nucleotide position (numeric))
import os
from TkinterDialog import TkinterDialog, Selections
from UsefulBioinformaticsFunctions import baseChromosomes
from UsefulFileSystemFunctions import Metadata, generateFilePath
from typing import IO


# Contains data on a single mutation position obtained by reading the next available line in a given file.
class MutationData:

    def __init__(self,file: IO):
        
        self.file = file # The file containing mutation data.

        # Read in the first line and check that it isn't empty.
        choppedUpLine = file.readline().strip().split()
        if len(choppedUpLine) == 0: raise ValueError("Error:  Empty mutation file given.")

        # Initialize class variables
        self.isEmpty = False # This variable keeps track of when EOF is reached.
        self.chromosome = choppedUpLine[0] # The chromosome that houses the mutation.
        self.position = int(choppedUpLine[1]) # The position of the mutation in its chromosome.
        self.strand = choppedUpLine[5] # Either '+' or '-' depending on which strand houses the mutation.

        # Make sure the mutation is in a valid chromosome.  Otherwise, something is wrong..
        if not self.chromosome in baseChromosomes:
            raise ValueError(choppedUpLine[0] + " is not a valid chromosome for the mutation trinuc file.")


    def readNextMutation(self):
        
        # Read in the next line.
        choppedUpLine = self.file.readline().strip().split()

        # Check if EOF has been reached.
        if len(choppedUpLine) == 0: 
            self.isEmpty = True
            self.chromosome = None
            self.position = None
            self.strand = None
            return

        # Assign variables
        self.chromosome = choppedUpLine[0]
        self.position = int(choppedUpLine[1])
        self.strand = choppedUpLine[5]

        # Make sure the mutation is in a valid chromosome.  Otherwise, something is wrong..
        if not self.chromosome in baseChromosomes:
            raise ValueError(choppedUpLine[0] + " is not a valid chromosome for the mutation trinuc file.")
        

# Contains data on a single nucleosome position obtained by reading the next available line in a given file.
class NucleosomeData:

    def __init__(self,file: IO):
        
        self.file = file # The file containing nucleosome data.

        # Read in the first line and check that it isn't empty.
        choppedUpLine = file.readline().strip().split()
        if len(choppedUpLine) == 0: raise ValueError("Error:  Empty nucleosome file given.")

        # Initialize class variables
        self.isEmpty = False # This variable keeps track of when EOF is reached.
        self.chromosome = choppedUpLine[0] # The chromosome that houses the nucleosome.
        self.dyadPosNegative74 = int(choppedUpLine[1]) # The position of the base pair at dyad position -74 for the current nucleosome.
            
    
    def readNextNucleosome(self):
        
        # Read in the next line.
        choppedUpLine = self.file.readline().strip().split()

        # Check if EOF has been reached.
        if len(choppedUpLine) == 0: 
            self.isEmpty = True
            self.chromosome = None
            self.dyadPosNegative74 = None
            return

        # Assign variables
        self.chromosome = choppedUpLine[0]
        self.dyadPosNegative74 = int(choppedUpLine[1])-74


# Takes a mutation object and nucleosome object which have unequal chromosomes and read through data until they are equal.
def reconcileChromosomes(mutation: MutationData, nucleosome: NucleosomeData):
    
    # Until the chromosomes are the same for both mutations and 
    while not mutation.chromosome == nucleosome.chromosome and not (mutation.isEmpty or nucleosome.isEmpty):
        if mutation.chromosome < nucleosome.chromosome: mutation.readNextMutation()
        else: nucleosome.readNextNucleosome()

    if not (nucleosome.isEmpty or mutation.isEmpty): print("Counting in",nucleosome.chromosome)


# Determines whether or not the given mutation is past the range of the given nucleosome.
def isMutationPastNucleosome(mutation: MutationData, nucleosome: NucleosomeData, linkerOffset):

    if mutation.isEmpty:
        return True
    elif mutation.position-nucleosome.dyadPosNegative74 > 147 + linkerOffset:
        return True
    elif not mutation.chromosome == nucleosome.chromosome:
        return True
    else: 
        return False


# Uses the given nucleosome position file and mutation file to count the number of mutations 
# at each dyad position (-73 to 73) for each strand.  Generates a new file to store these results.
# Optionally, some amount of surrounding linker DNA may be included in the analysis as "linker offset"
# NOTE:  It is VITAL that both files are sorted, first by chromosome number and then by starting coordinate.
#        (Sorted first by chromosome (string) and then by nucleotide position (numeric))
#        This code is pretty slick, but it will crash and burn and give you a heap of garbage as output if the inputs aren't sorted.
def generateCountsFile(mutationFilePath, nucPosFilePath, nucleosomeMutationCountsFilePath, linkerOffset):
    print("Counting mutations at each nucleosome position...")

    # Dictionaries holding the number of mutations found at each dyad position from -73 to 73 for each strand.
    minusStrandMutationCounts = dict() 
    plusStrandMutationCounts = dict()
    for i in range(-73-linkerOffset,74+linkerOffset):
        minusStrandMutationCounts[i] = 0
        plusStrandMutationCounts[i] = 0

    # Keeps track of mutations in the current and previous nucleosomes to catch mutations that span overlapping nucleosomes.
    # Each muation is a tuple containing the mutation position and the strand it is housed on.
    mutationsInCurrentNucleosome = list()
    mutationsInPreviousNucleosome = list()

    # A function which checks to see if a given mutation falls within the given nucleosome and if it does, adds it to the given list.
    def addMutationIfInNucleosome(mutationPosition, mutationStrand, dyadPosNegative74, linkerOffset):

        if mutationPosition - dyadPosNegative74 > 0 - linkerOffset:
            if mutationStrand == '+': plusStrandMutationCounts[mutationPosition-dyadPosNegative74 - 74] += 1
            elif mutationStrand == '-': minusStrandMutationCounts[mutationPosition-dyadPosNegative74 - 74] += 1
            else:  raise ValueError("Error:  No strandedness found for mutation")
            
            # Add the mutation to the list of mutations in the current nucleosome
            mutationsInCurrentNucleosome.append((mutationPosition,mutationStrand))

    # Open the mutation and nucleosome files to compare against one another.
    with open(mutationFilePath, 'r') as mutationFile:
        with open(nucPosFilePath,'r') as nucPosFile:

            #Get data on the first mutation and nucleosome and reconcile their chromosomes if necessary to start things off.
            currentMutation = MutationData(mutationFile)            
            currentNucleosome = NucleosomeData(nucPosFile) 
            if currentNucleosome.chromosome == currentMutation.chromosome: print("Counting in",currentNucleosome.chromosome)
            else: reconcileChromosomes(currentMutation, currentNucleosome)

            # The core loop goes through each nucleosome one at a time and checks mutation positions against it until 
            # one exceeds its rightmost position or is on a different chromosome.  
            # Then, the next nucleosome is checked until either no mutations or no nucleosomes remain.
            while not (currentNucleosome.isEmpty or currentMutation.isEmpty):
                    
                # Check to see if any of the mutations in the last nucleosome are present in the current nucleosome due to overlap.
                mutationsInPreviousNucleosome = mutationsInCurrentNucleosome.copy()
                mutationsInCurrentNucleosome.clear()
                for mutation in mutationsInPreviousNucleosome:
                    addMutationIfInNucleosome(mutation[0],mutation[1],currentNucleosome.dyadPosNegative74,linkerOffset)

                # Read mutations until the mutation is past the range of the current nucleosome.
                while not isMutationPastNucleosome(currentMutation,currentNucleosome,linkerOffset):

                    # Check and see if we need to add the mutation to our lists.
                    addMutationIfInNucleosome(currentMutation.position,currentMutation.strand,currentNucleosome.dyadPosNegative74,linkerOffset)

                    #Get data on the next mutation.
                    currentMutation.readNextMutation()

                # If no recnociliation is necessary at the moment, read in a new nucleosome.
                if currentMutation.chromosome == currentNucleosome.chromosome: currentNucleosome.readNextNucleosome()

                # Reconcile the mutation data and chromosome data if necessary to be sure
                # that they are looking at the same chromosome for the next iteration
                if not currentMutation.chromosome == currentNucleosome.chromosome:
                    reconcileChromosomes(currentMutation, currentNucleosome)
                    mutationsInCurrentNucleosome.clear() # Nucleosomes can't overlap on different chromosomes.

    # Write the results to the output file.
    with open(nucleosomeMutationCountsFilePath,'w') as nucleosomeMutationCountsFile:
        
        # Write the headers to the file.
        nucleosomeMutationCountsFile.write('\t'.join(("Dyad_Pos","Plus_Strand_Counts",
                                            "Minus_Strand_Counts","Both_Strands_Counts")) + '\n')
        
        # Write the data.
        for i in range(-73-linkerOffset,74+linkerOffset):
            nucleosomeMutationCountsFile.write('\t'.join((str(i), str(plusStrandMutationCounts[i]), str(minusStrandMutationCounts[i]), 
                                                str(plusStrandMutationCounts[i] + minusStrandMutationCounts[i]))) + '\n')


def countNucleosomePositionMutations(mutationFilePaths, linkerOffset):

    nucleosomeMutationCountsFilePaths = list() # A list of paths to the output files generated by the function

    # Loop through each given mutation file path, creating a corresponding nucleosome mutation count file for each.
    for mutationFilePath in mutationFilePaths:

        print("\nWorking with",os.path.split(mutationFilePath)[1])

        # Make sure we have the expected file type.
        if not "context_mutations" in os.path.basename(mutationFilePath): 
            raise ValueError("Mutation file should have \"context_mutations\" in the name.")

        # Get metadata and use it to generate a path to the nucleosome positions file.
        metadata = Metadata(mutationFilePath)

        # Generate the output file path
        nucleosomeMutationCountsFilePath = generateFilePath(directory = metadata.directory,
                                                            dataGroup = metadata.dataGroupName,
                                                            linkerOffset = linkerOffset, fileExtension = ".tsv",
                                                            dataType = "raw_nucleosome_mutation_counts")

        # Ready, set, go!
        generateCountsFile(mutationFilePath, metadata.baseNucPosFilePath, 
                           nucleosomeMutationCountsFilePath, linkerOffset)

        nucleosomeMutationCountsFilePaths.append(nucleosomeMutationCountsFilePath)

    return nucleosomeMutationCountsFilePaths


if __name__ == "__main__":

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(__file__),"..","data"))
    dialog.createMultipleFileSelector("Mutation Files:",0,"_context_mutations.bed",("Bed Files",".bed"))
    dialog.createCheckbox("Include 30 bp linker DNA on either side.",1, 0, 2)
    dialog.createReturnButton(2,0,2)
    dialog.createQuitButton(2,2,2)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    mutationFilePaths = list(selections.getFilePathGroups())[0] # A list of mutation file paths
    includeLinker = list(selections.getToggleStates())[0]

    if includeLinker: linkerOffset = 30
    else: linkerOffset = 0

    countNucleosomePositionMutations(mutationFilePaths, linkerOffset)