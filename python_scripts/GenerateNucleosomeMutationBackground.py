# This script, when given a mutation background file, a genome file, and a file of strongly positioned nucleosome coordinates,
# generates a background file with the expected mutations at each dyad position from -73 to 73 (inclusive).

from TkinterDialog import TkinterDialog, Selections
import os, subprocess
from typing import Dict
from UsefulBioinformaticsFunctions import bedToFasta, reverseCompliment

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


# The user could have given a nucleosome positioning file in one of three reasonable formats... Parse away!
def parseStrongPosNucleosomeData(strongPosNucleosomeFilePath,strongPosNucleosomeFastaFilePath):

    # The fasta format should already be formatted correctly, so no further formatting should be necessary.
    if strongPosNucleosomeFilePath.endswith(".fa"):
        strongPosNucleosomeFastaFilePath = strongPosNucleosomeFilePath
        print("Nucleosome positioning data given in fasta format.")

    elif os.path.exists(strongPosNucleosomeFastaFilePath):
        print("Found fasta file associated with given nucleosome positioning file!") 

    # It's not a fasta file, so check if the given bed (hopefully) file needs to be expanded, then generate the fasta file.
    else:
        with open(strongPosNucleosomeFilePath, 'r') as strongPosNucleosomeFile:
            strongPosNucleosomeFile.readline()
            choppedUpLine = strongPosNucleosomeFile.readline().strip().split('\t')

            strongPosNucleosomeExpandedFilePath = strongPosNucleosomeFilePath # Assume the file is expanded for now.

            # If it is not actually expanded, expand it!
            if int(choppedUpLine[2]) - int(choppedUpLine[1]) == 1:
                print("Found unexpanded nucleosome coordinates.")
                strongPosNucleosomeExpandedFilePath = expandNucleosomeCoordinates(strongPosNucleosomeFilePath)
            # If it is expanded, we're already ready to convert to fasta.
            elif not int(choppedUpLine[2]) - int(choppedUpLine[1]) == 149:
                raise ValueError("Error: Strongly positioned nucleosome data is not in the expected format.\n" +  
                    "Each coordinate should contain the central base pair in the dyad only, " + 
                    "or in addition, exactly 74 bp on either side.")

            # Convert the file to fasta format!
            print("Converting strongly positioned nucleosome coordinates to fasta file...")
            bedToFasta(strongPosNucleosomeExpandedFilePath,genomeFilePath,strongPosNucleosomeFastaFilePath)


# This function returns a dictionary of mutation rates for all trinucleotide contexts in the associated genome.
def getTrinucMutationRate(mutationBackgroundFilePath):

    trinucMutationRate = dict()

    # Open the file and read the lines into the dictionary.
    with open(mutationBackgroundFilePath, 'r') as mutationBackgroundFile:
        mutationBackgroundFile.readline() # Skip the line with headers.
        for line in mutationBackgroundFile:
            choppedUpLine = line.strip().split('\t')
            if len(choppedUpLine[0]) != 3:
                raise ValueError("Error: found \"" + choppedUpLine[0] + "\" in " + mutationBackgroundFilePath + " but expected a trinucleotide sequence.")
            trinucMutationRate[choppedUpLine[0]] = float(choppedUpLine[1])
    
    return trinucMutationRate


# This function generates a file of trinuc counts in the genome for each dyad position.
def generateDyadPosTrinucCounts(strongPosNucleosomeFastaFilePath, dyadPosTrinucCountsFilePath):

    # Dictionary of trinuc counts for every dyad position from -73 to 73 on the plus strand.
    plusStrandNucleosomeDyadPosTrinucCounts = dict() 

    # Hash table of observed trinucleotide sequences.
    observedTrinucs = dict()

    # Initialize the dictionary
    for dyadPos in range(-73,74): 
        plusStrandNucleosomeDyadPosTrinucCounts[dyadPos] = dict()

    # Initialize the counter for dyad position
    dyadPos = 74

    # Read through the file, counging trinucs for every dyad position to the running total in the dictionary
    with open(strongPosNucleosomeFastaFilePath, 'r') as strongPosNucleosomeFastaFile:

        for line in strongPosNucleosomeFastaFile:

            # Go away, whitespace! >:(
            line = line.strip()

            # Check and see if we're starting a new nucleosome with this line.
            if line.startswith(">"):

                # Make sure we didn't end the last nucleosome prematurely.
                if not dyadPos == 74:
                    raise ValueError("Error in nucleosome positioning file.  " +  
                                    "Last nucleosome did not give trinucleotides up to position 73.  "+
                                    "Last dyad position was " + str(dyadPos-1))

                # Reset line leftovers and dyad position
                lineLeftovers = ''
                dyadPos = -73

            else:

                # Add any leftovers to the current line.
                currentSequence = lineLeftovers + line

                # Count all available trinuc sequences.
                for i in range(0,len(currentSequence)-2):

                    trinuc = currentSequence[i:i+3]
                    if trinuc not in observedTrinucs: observedTrinucs[trinuc] = None
                    
                    # Add the trinuc to the running total in the counts dictionary.
                    plusStrandNucleosomeDyadPosTrinucCounts[dyadPos][trinuc] = plusStrandNucleosomeDyadPosTrinucCounts[dyadPos].setdefault(trinuc, 0) + 1

                    # Increment the dyadPos counter
                    dyadPos += 1

                # Don't forget about the leftovers for this line!
                lineLeftovers = line[-2:]

    # Write the trinuc counts for every dyad position on the plus strand in the output file
    with open(dyadPosTrinucCountsFilePath, 'w') as dyadPosTrinucCountsFile:

        # Write the header for the file
        dyadPosTrinucCountsFile.write("Dyad_Pos\t" + '\t'.join(observedTrinucs.keys()) + '\n')

        for dyadPos in plusStrandNucleosomeDyadPosTrinucCounts.keys():
            
            dyadPosTrinucCountsFile.write(str(dyadPos))

            for trinuc in observedTrinucs.keys():
                if trinuc in plusStrandNucleosomeDyadPosTrinucCounts[dyadPos]:
                    dyadPosTrinucCountsFile.write('\t' + str(plusStrandNucleosomeDyadPosTrinucCounts[dyadPos][trinuc]))
                else: dyadPosTrinucCountsFile.write('\t0')

            dyadPosTrinucCountsFile.write('\n')
            
    
# This function retrieves the trinuc counts for each dyad position in a genome from a given file.
# The data is returned as a dictionary of dictionaries, with the first key being dyad position and the second
# being a trinucleotide sequence.
def getDyadPosTrinucCounts(dyadPosTrinucCountsFilePath):
    with open(dyadPosTrinucCountsFilePath, 'r') as dyadPosTrinucCountsFile:
        
        trinucSequences = list()
        headersHaveBeenRead = False
        dyadPosTrinucCounts= dict()

        for line in dyadPosTrinucCountsFile:

            if not headersHaveBeenRead:
                headersHaveBeenRead = True
                trinucSequences = line.strip().split('\t')[1:]
    
            else: 
                choppedUpLine = line.strip().split('\t')
                dyadPos = int(choppedUpLine[0])
                dyadPosTrinucCounts[dyadPos] = dict()
                for i,trinuc in enumerate(trinucSequences):
                    dyadPosTrinucCounts[dyadPos][trinuc] = int(choppedUpLine[i+1])

    return dyadPosTrinucCounts


# This function generates a nucleosome mutation background file from a general mutation background file
# and a file of strongly positioned nucleosome coordinates.
def generateNucleosomeMutationBackground(dyadPosTrinucCountsFilePath, mutationBackgroundFilePath,
    nucleosomeMutationBackgroundFilePath):
    print("Generating nucleosome mutation background file...")

    # Dictionaries of expected mutations for every dyad position from -73 to 73, one for each strand.
    plusStrandNucleosomeMutationBackground = dict() 
    minusStrandNucleosomeMutationBackground = dict()

    # Initialize the dictionary
    for dyadPos in range(-73,74): 
        plusStrandNucleosomeMutationBackground[dyadPos] = 0
        minusStrandNucleosomeMutationBackground[dyadPos] = 0

    # Get the trinuc mutation rate dictionary.
    trinucMutationRate = getTrinucMutationRate(mutationBackgroundFilePath)
    # Get the dyad position trinuc counts dictionary.
    dyadPosTrinucCounts = getDyadPosTrinucCounts(dyadPosTrinucCountsFilePath)

    # Calculate the expected mutation rates for each dyad position based on the trinuc counts at that position and that trinuc's mutation rate
    for dyadPos in dyadPosTrinucCounts:

        for trinuc in dyadPosTrinucCounts[dyadPos]:

            reverseTrinuc = reverseCompliment(trinuc)
            
            # Add the trinuc's mutation rate to the running total in the background dictionaries.
            plusStrandNucleosomeMutationBackground[dyadPos] += trinucMutationRate[trinuc] * dyadPosTrinucCounts[dyadPos][trinuc]
            minusStrandNucleosomeMutationBackground[dyadPos] += trinucMutationRate[reverseTrinuc] * dyadPosTrinucCounts[dyadPos][trinuc]

    # Write the results of the dictionary to the nucleosome mutation background file.
    with open(nucleosomeMutationBackgroundFilePath, 'w') as nucleosomeMutationBackgroundFile:

        # Write the headers for the data.
        headers = '\t'.join(("Dyad_Position","Expected_Mutations_Plus_Strand",
        "Expected_Mutations_Minus_Strand","Expected_Mutations_Both_Strands"))

        nucleosomeMutationBackgroundFile.write(headers + '\n')
        
        # Write the data for each dyad position.
        for dyadPos in range(-73,74):

            dataRow = '\t'.join((str(dyadPos),str(plusStrandNucleosomeMutationBackground[dyadPos]),
            str(minusStrandNucleosomeMutationBackground[dyadPos]),
            str(plusStrandNucleosomeMutationBackground[dyadPos] + minusStrandNucleosomeMutationBackground[dyadPos])))

            nucleosomeMutationBackgroundFile.write(dataRow + '\n')


#Create the Tkinter UI
dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(__file__),"..","data"))
dialog.createMultipleFileSelector("Mutation Background Files:",0,"mutation_background.tsv",("Tab Seperated Values Files",".tsv"))
dialog.createFileSelector("Genome Fasta File:",1,("Fasta Files",".fa"))
dialog.createFileSelector("Strongly Positioned Nucleosome File:",2,("Bed or Fasta Files",".bed .fa"))
dialog.createReturnButton(3,0,2)
dialog.createQuitButton(3,2,2)

# Run the UI
dialog.mainloop()

# If no input was received (i.e. the UI was terminated prematurely), then quit!
if dialog.selections is None: quit()

# Get the user's input from the dialog.
selections: Selections = dialog.selections
filePaths = list(selections.getFilePaths())
mutationBackgroundFilePaths = list(selections.getFilePathGroups())[0] # A list of mutation background file paths
genomeFilePath = filePaths[0] # The path to the genome fasta file
strongPosNucleosomeFilePath: str = filePaths[1] # The path to the file containing strongly positioned nucleosomes.

# Generate a path to the fasta file of strongly positioned nucleosome coordinates. 
if strongPosNucleosomeFilePath.rsplit('.',1)[0].endswith("_expansion"):
    strongPosNucleosomeFastaFilePath = strongPosNucleosomeFilePath.rsplit("_expansion",1)[0] + ".fa"
else:
    strongPosNucleosomeFastaFilePath = strongPosNucleosomeFilePath.rsplit('.',1)[0] + ".fa"

# Generate the path to the tsv file of dyad position trinuc counts
dyadPosTrinucCountsFilePath = strongPosNucleosomeFastaFilePath.rsplit(".",1)[0] + "_dyad_pos_trinuc_counts.tsv"

# Make sure we have a fasta file for strongly positioned nucleosome coordinates
parseStrongPosNucleosomeData(strongPosNucleosomeFilePath,strongPosNucleosomeFastaFilePath)

# Make sure we have a tsv file with trinuc counts at each dyad position.
if not os.path.exists(dyadPosTrinucCountsFilePath): 
    print("Dyad position trinuc counts file not found at",dyadPosTrinucCountsFilePath)
    print("Generating genome wide dyad position trinuc counts file...")
    generateDyadPosTrinucCounts(strongPosNucleosomeFastaFilePath, dyadPosTrinucCountsFilePath)

# Loop through each given mutation background file path, creating a corresponding nucleosome mutation background for each.
for mutationBackgroundFilePath in mutationBackgroundFilePaths:

    mutationBackgroundFileName = os.path.split(mutationBackgroundFilePath)[1]
    print("\nWorking with file:",mutationBackgroundFileName)
    if not "_mutation_background" in mutationBackgroundFileName: 
        raise ValueError("Error, expected file with \"_mutation_background\" in the name.")

    # Get some information on the file system and generate file paths for convenience.
    workingDirectory = os.path.dirname(mutationBackgroundFilePath) # The working directory for the current data "group"

    # The name of the mutation data set the background pertains to.
    mutationDataGroup = mutationBackgroundFileName.split("_mutation_background",1)[0]

    # A path to the final output file.
    nucleosomeMutationBackgroundFilePath = os.path.join(workingDirectory,mutationDataGroup+"_nucleosome_mutation_background.tsv")

    # Generate the nucleosome mutation background file!
    generateNucleosomeMutationBackground(dyadPosTrinucCountsFilePath,mutationBackgroundFilePath,nucleosomeMutationBackgroundFilePath)