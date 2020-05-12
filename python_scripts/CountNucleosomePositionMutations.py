# This script takes a mutation file and the coordinates for DNA bases around nucleosomes and calculates how many
# mutations occured at each dyad position for each (and both) strands.
# NOTE:  Both input files must be sorted for this script to run properly. 
#        (Sorted first by chromosome (string) and then by nucleotide position (numeric))
import os
import time
from TkinterDialog import TkinterDialog, Selections
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


    def readNextMutation(self):
        
        # Read in the next line.
        choppedUpLine = self.file.readline().strip().split()

        # Check if EOF has been reached.
        if len(choppedUpLine) == 0: 
            self.isEmpty = True
            return

        # Assign variables
        self.chromosome = choppedUpLine[0]
        self.position = int(choppedUpLine[1])
        self.strand = choppedUpLine[5]
        

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
            return

        # Assign variables
        self.chromosome = choppedUpLine[0]
        self.dyadPosNegative74 = int(choppedUpLine[1])


# This function takes a bed file of strongly positioned nucleosomes and expands their coordinates to encompass
# 74 bases on either side of the dyad. (in order to get trinucleotide sequences for positions -73 to 73)
# Returns the file path to the newly expanded file.
def expandNucleosomeCoordinates(strongPosNucleosomeFilePath):

    # Generate a file path for the expanded nucleosome coordinate file.
    strongPosNucleosomeExpansionFilePath = strongPosNucleosomeFilePath.rsplit('.',1)[0] + "_expansion.bed"

    # Check to see if the expansion file already exists.
    if os.path.exists(strongPosNucleosomeExpansionFilePath):
        print ("Expanded nucleosome coordinate file already exists.  Not re-expanding.")
        return strongPosNucleosomeExpansionFilePath

    print("Expanding nucleosome coordinates...")

    # Open the files.
    with open(strongPosNucleosomeFilePath,'r') as strongPosNucleosomeFile:
        with open(strongPosNucleosomeExpansionFilePath, 'w') as strongPosNucleosomeExpansionFile:

            # Write the expanded positions to the new file, one line at a time.
            for line in strongPosNucleosomeFile:
                choppedUpLine = line.strip().split('\t')
                choppedUpLine[1] = str(int(choppedUpLine[1]) - 74)
                choppedUpLine[2] = str(int(choppedUpLine[2]) + 74)

                strongPosNucleosomeExpansionFile.write('\t'.join(choppedUpLine) + '\n')

    # Return the file path to the newly expanded file.
    return strongPosNucleosomeExpansionFilePath


# Takes a mutation object and nucleosome object which have unequal chromosomes and read through data until they are equal.
def reconcileChromosomes(mutation: MutationData, nucleosome: NucleosomeData):
    
    # Until the chromosomes are the same for both mutations and 
    while not mutation.chromosome == nucleosome.chromosome and not (mutation.isEmpty or nucleosome.isEmpty):
        if mutation.chromosome < nucleosome.chromosome: mutation.readNextMutation()
        else: nucleosome.readNextNucleosome()

    if not (nucleosome.isEmpty or mutation.isEmpty): print("Counting in",nucleosome.chromosome)


# Determines whether or not the given mutation is past the range of the given nucleosome.
def isMutationPastNucleosome(mutation: MutationData, nucleosome: NucleosomeData):
    if mutation.position-nucleosome.dyadPosNegative74 > 147:
        return True
    elif not mutation.chromosome == nucleosome.chromosome:
        return True
    else: 
        return False


# Checks to see whether or not the given nucleosome positioning data needs to be expanded, 
# and calls the necessary function to expand it if needed.
def parseStrongPosNucleosomeData(strongPosNucleosomeFilePath):

    # Check if the given bed (hopefully) file needs to be expanded, and do it if it does.
    with open(strongPosNucleosomeFilePath, 'r') as strongPosNucleosomeFile:
        strongPosNucleosomeFile.readline()
        choppedUpLine = strongPosNucleosomeFile.readline().strip().split('\t')

        # If it is not actually expanded, expand it!
        if int(choppedUpLine[2]) - int(choppedUpLine[1]) == 1:
            print("Unexpanded nucleosome coordinates were given.")
            strongPosNucleosomeFilePath = expandNucleosomeCoordinates(strongPosNucleosomeFilePath)
        # If it isn't "not" expanded, is it expanded properly?
        elif not int(choppedUpLine[2]) - int(choppedUpLine[1]) == 149:
            raise ValueError("Error: Strongly positioned nucleosome data is not in the expected format.\n" +  
                "Each coordinate should contain the central base pair in the dyad only, or in addition, exactly 74 bp on either side.")
        else:
            print("Given nucleosome coordinates are already expanded.")

    print("Using nucleosome coordinate file", os.path.split(strongPosNucleosomeFilePath)[1])
    return strongPosNucleosomeFilePath


# Uses the given nucleosome position file and mutation file to count the number of mutations 
# at each dyad position (-73 to 73) for each strand.  Generates a new file to store these results.
# NOTE:  It is VITAL that both files are sorted, first by chromosome number and then by starting coordinate.
#        (Sorted first by chromosome (string) and then by nucleotide position (numeric))
#        This code is pretty slick, but it will crash and burn and give you a heap of garbage as output if the inputs aren't sorted.
def countNucleosomePositionMutations(mutationFilePath, strongPosNucleosomeFilePath, nucleosomeMutationCountsFilePath):
    print("Counting mutations at each nucleosome position...")

    # Dictionaries holding the number of mutations found at each dyad position from -73 to 73 for each strand.
    minusStrandMutationCounts = dict() 
    plusStrandMutationCounts = dict()
    for i in range(-73,74):
        minusStrandMutationCounts[i] = 0
        plusStrandMutationCounts[i] = 0

    # Keeps track of mutations in the current and previous nucleosomes to catch mutations that span overlapping nucleosomes.
    # Each muation is a tuple containing the mutation position and the strand it is housed on.
    mutationsInCurrentNucleosome = list()
    mutationsInPreviousNucleosome = list()

    # A function which checks to see if a given mutation falls within the given nucleosome and if it does, adds it to the given list.
    def addMutationIfInNucleosome(mutationPosition, mutationStrand, dyadPosNegative74):

        if mutationPosition - dyadPosNegative74 > 0:
            if mutationStrand == '+': plusStrandMutationCounts[mutationPosition-dyadPosNegative74 - 74] += 1
            elif mutationStrand == '-': minusStrandMutationCounts[mutationPosition-dyadPosNegative74 - 74] += 1
            else:  raise ValueError("Error:  No strandedness found for mutation")
            
            # Add the mutation to the list of mutations in the current nucleosome
            mutationsInCurrentNucleosome.append((mutationPosition,mutationStrand))

    # Open the mutation and nucleosome files to compare against one another.
    with open(mutationFilePath, 'r') as mutationFile:
        with open(strongPosNucleosomeFilePath,'r') as strongPosNucleosomeFile:

            #Get data on the first mutation and nucleosome and reconcile their chromosomes if necessary to start things off.
            currentMutation = MutationData(mutationFile)            
            currentNucleosome = NucleosomeData(strongPosNucleosomeFile) 
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
                    addMutationIfInNucleosome(mutation[0],mutation[1],currentNucleosome.dyadPosNegative74)

                # Read mutations until the mutation is past the range of the current nucleosome.
                while not isMutationPastNucleosome(currentMutation,currentNucleosome):

                    # Check and see if we need to add the mutation to our lists.
                    addMutationIfInNucleosome(currentMutation.position,currentMutation.strand,currentNucleosome.dyadPosNegative74)

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
        for i in range(-73,74):
            nucleosomeMutationCountsFile.write('\t'.join((str(i), str(plusStrandMutationCounts[i]), str(minusStrandMutationCounts[i]), 
                                                str(plusStrandMutationCounts[i] + minusStrandMutationCounts[i]))) + '\n')


#Create the Tkinter UI
dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(__file__),"..","data"))
dialog.createMultipleFileSelector("Mutation Files:",0,"trinuc_context.bed",("Bed Files",".bed"))
dialog.createFileSelector("Strongly Positioned Nucleosome File:",1,("Bed Files",".bed"))
dialog.createReturnButton(2,0,2)
dialog.createQuitButton(2,2,2)

# Run the UI
dialog.mainloop()

# If no input was received (i.e. the UI was terminated prematurely), then quit!
if dialog.selections is None: quit()

# Get the user's input from the dialog.
selections: Selections = dialog.selections
mutationFilePaths = list(selections.getFilePathGroups())[0] # A list of mutation file paths
# The path to the file containing strongly positioned nucleosomes.
strongPosNucleosomeFilePath: str = list(selections.getIndividualFilePaths())[0] 

# Make sure we have a path the expanded nucleosome position file.
strongPosNucleosomeFilePath = parseStrongPosNucleosomeData(strongPosNucleosomeFilePath)

# Loop through each given mutation file path, creating a corresponding nucleosome mutation count file for each.
for mutationFilePath in mutationFilePaths:

    print("\nWorking with",os.path.split(mutationFilePath)[1])

    # Make the file path for the output file
    workingDirectory = os.path.split(mutationFilePath)[0]
    mutationGroupName = os.path.split(mutationFilePath)[1].rsplit("_trinuc",1)[0].rsplit("_singlenuc",1)[0]
    # Double check that the mutation group name was generated correctly.
    if '.' in mutationGroupName: raise ValueError("Error, expected mutation file with \"trinuc\" or \"singlenuc\" in the name.")
    nucleosomeMutationCountsFilePath = os.path.join(workingDirectory,mutationGroupName+"_nucleosome_mutation_counts.tsv")

    # Ready, set, go!
    countNucleosomePositionMutations(mutationFilePath, strongPosNucleosomeFilePath, nucleosomeMutationCountsFilePath)