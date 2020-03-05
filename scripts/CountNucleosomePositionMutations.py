# This script takes a mutation file and the coordinates for DNA bases around nucleosomes and calculates how many
# mutations occured at each dyad position for each (and both) strands.
# NOTE:  Both input files must be sorted for this script to run properly.
import os
import time
from TkinterDialog import TkinterDialog, Selections

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
#        This code is pretty slick, but it will crash and burn and give you a heap of garbage as output if the inputs aren't sorted.
# ANOTHER NOTE:  This function got out of hand fast.  I probably should have used some classes to clean up implementation,
#                but I didn't.  Sorry future me (or whoever is reading this). 
#                You should probably implement Mutation and Nucleosome classes that initialize by reading a line from a given file.
def countNucleosomePositionMutations(mutationFilePath, strongPosNucleosomeFilePath, nucleosomeMutationCountsFilePath):
    print("Counting mutations at each nucleosome position...")

    # Dictionaries holding the number of mutations found at each dyad position from -73 to 73 for each strand.
    minusStrandMutationCounts = dict() 
    plusStrandMutationCounts = dict()
    for i in range(-73,74):
        minusStrandMutationCounts[i] = 0
        plusStrandMutationCounts[i] = 0

    mutationChromosome = '' # The chromosome that houses the mutation being checked against the current nucleosome.
    mutationPos = 0 # The position of the mutation being checked against the current nucleosome.
    strand = '' # Either '+' or '-' depending on which strand houses the mutation.

    nucleosomeChromosome = '' # The chromosome that houses the nucleosome currently being checked against.
    dyadPosNegative74 = 0 # The position of the base pair at dyad position -74 for the current nucleosome.

    chromosomesWithMutations = list() # A list of all the chromosomes containing mutations.

    # Keeps track of mutations in the current and previous nucleosomes to catch mutations that span overlapping nucleosomes.
    # Each item in the list is a tuple containing the mutation position followed by its strand.
    mutationsInCurrentNucleosome = list()
    mutationsInPreviousNucleosome = list()
    
    # Record which chromosomes have at least one mutation in them.
    with open(mutationFilePath, 'r') as mutationFile:
        for line in mutationFile:
            if line.split()[0] not in chromosomesWithMutations: chromosomesWithMutations.append(line.split()[0])

    # A function to check whether or not the given mutation could be in a nucleosome further down the sorted list.
    def isMutationPastNucleosome():

        # First, check to see if we are on the same chromosome or not.
        if not mutationChromosome == nucleosomeChromosome:
            return True
        # Then, check to see if we are past dyad position 73.
        elif mutationPos-dyadPosNegative74 > 147:
            return True
        else:
            return False

    with open(mutationFilePath, 'r') as mutationFile:
        with open(strongPosNucleosomeFilePath,'r') as strongPosNucleosomeFile:

            # A function which retrieves data from the next line in the mutation data file.
            def readMutationData():

                choppedUpLine = mutationFile.readline().strip().split('\t')

                # Set the value for the chromosome housing the mutation, or set it as an empty string if EOF.
                nonlocal mutationChromosome
                mutationChromosome = choppedUpLine[0]
                if mutationChromosome == '': return

                # Set the values for the position of the mutation and the strand it is on.
                nonlocal mutationPos
                mutationPos = int(choppedUpLine[1])
                nonlocal strand
                strand = choppedUpLine[5]

            # A function which retrieves data from the next line in the nucleosome data file.
            def readNucleosomeData():

                choppedUpLine = strongPosNucleosomeFile.readline().strip().split('\t')
                
                nonlocal nucleosomeChromosome
                nonlocal dyadPosNegative74
                nonlocal mutationsInCurrentNucleosome

                # Check and see if we have moved on to a new chromosome.
                # If so, we need to find the next chromosome that has both recorded nucleosome positions and mutations!
                if not choppedUpLine[0] == nucleosomeChromosome:
                    nucleosomeChromosome = choppedUpLine[0]
                    if nucleosomeChromosome == '': return

                    print("Searching for mutations in nucleosomes in",nucleosomeChromosome)

                    # If we are switching chromosomes, overlap is impossible between this and the last nucleosome.
                    mutationsInCurrentNucleosome.clear()

                    # Update dyadPos so it does not overwrite the values set by later function calls which may be invoked
                    # as this function searches for the next chromosome with both mutations and nucleosome positions
                    dyadPosNegative74 = int(choppedUpLine[1])

                    # If there are no mutations on the new chromosome, keep going until we get to a chromosome that has some.
                    while not nucleosomeChromosome in chromosomesWithMutations: readNucleosomeData()

                    # Read down the list of mutations until we get to the first one in the current chromosome, or reach EOF.
                    while not mutationChromosome == nucleosomeChromosome and not mutationChromosome == '': readMutationData()

                # We don't want to overwrite the global dyad position after just changing chromosomes as another
                # call of this function may have updated it to a more current value.
                else: dyadPosNegative74 = int(choppedUpLine[1])


            #Get data on the first mutation and nucleosome to start things off.
            readMutationData()             
            readNucleosomeData()   

            # The core loop goes through each nucleosome one at a time and checks mutation positions against it until 
            # one exceeds its rightmost position or is on a different chromosome.  Then, the next nucleosome is checked, etc.
            while not nucleosomeChromosome == '':

                # Make sure that none of the mutations in the last nucleosome are present in the current nucleosome due to overlap.
                mutationsInPreviousNucleosome = mutationsInCurrentNucleosome.copy()
                mutationsInCurrentNucleosome.clear()
                for mutation in mutationsInPreviousNucleosome:
                    if mutation[0] - dyadPosNegative74 > 0:
                        if mutation[1] == '+': plusStrandMutationCounts[mutation[0]-dyadPosNegative74 - 74] += 1
                        elif mutation[1] == '-': minusStrandMutationCounts[mutation[0]-dyadPosNegative74 - 74] += 1
                        else:  raise ValueError("Error:  No strandedness found for mutation")
                        mutationsInCurrentNucleosome.append(mutation)

                # Make sure the mutation isn't past the range of the current nucleosome.
                while not isMutationPastNucleosome():

                    # Look for mutations that fall within dyad positions -73 to 73, and add them to the counts!
                    if mutationPos - dyadPosNegative74 > 0:
                        if strand == '+': plusStrandMutationCounts[mutationPos-dyadPosNegative74 - 74] += 1
                        elif strand == '-': minusStrandMutationCounts[mutationPos-dyadPosNegative74 - 74] += 1
                        else:  raise ValueError("Error:  No strandedness found for mutation")
                        
                        # Add the mutation to the list of mutations in the current nucleosome
                        mutationsInCurrentNucleosome.append((mutationPos,strand))

                        # print("Mutation at position",mutationPos,"in",mutationChromosome,"on the",strand,"strand",
                        #         " was found at dyad position", str(mutationPos-dyadPosNegative74-74))

                    #Get data on the next mutation.
                    readMutationData()

                # Read in the next nucleosome
                readNucleosomeData()

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
dialog.createMultipleFileSelector("Mutation Files:",0,("Bed Files",".bed"))
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

    print("Working with",os.path.split(mutationFilePath)[1])

    # Make the file path for the output file
    workingDirectory = os.path.split(mutationFilePath)[0]
    mutationGroupName = os.path.split(mutationFilePath)[1].rsplit("_trinuc",1)[0].rsplit("_singlenuc",1)[0]
    # Double check that the mutation group name was generated correctly.
    if '.' in mutationGroupName: raise ValueError("Error, expected mutation file with \"trinuc\" or \"singlenuc\" in the name.")
    nucleosomeMutationCountsFilePath = os.path.join(workingDirectory,mutationGroupName+"_nucleosome_mutation_counts.txt")

    # Ready, set, go!
    countNucleosomePositionMutations(mutationFilePath, strongPosNucleosomeFilePath, nucleosomeMutationCountsFilePath)