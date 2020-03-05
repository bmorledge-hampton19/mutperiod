# This script, when given a mutation background file, a genome file, and a file of strongly positioned nucleosome coordinates,
# generates a background file with the expected mutations at each dyad position from -73 to 73 (inclusive).

from TkinterDialog import TkinterDialog, Selections
import os, subprocess
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
            trinucMutationRate[choppedUpLine[0]] = float(choppedUpLine[1])
    
    return trinucMutationRate


# This function takes generates a nucleosome mutation background file from a general mutation background file
# and a file of strongly positioned nucleosome coordinates.
def generateNucleosomeMutationBackground(strongPosNucleosomeFastaFilePath, mutationBackgroundFilePath,
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
    # Initialize the counter for dyad position
    dyadPos = 74

    # Read through the file, adding the mutation rate for every dyad position to the running total in the dictionary
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
                    reverseTrinuc = reverseCompliment(trinuc)
                    
                    # Add the trinuc's mutation rate to the running total in the background dictionaries.
                    plusStrandNucleosomeMutationBackground[dyadPos] += trinucMutationRate[trinuc]
                    minusStrandNucleosomeMutationBackground[dyadPos] += trinucMutationRate[reverseTrinuc]

                    # Increment the dyadPos counter
                    dyadPos += 1

                # Don't forget about the leftovers for this line!
                lineLeftovers = line[-2:]

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
dialog.createMultipleFileSelector("Mutation Background Files:",0,("Text Files",".txt"))
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

# A path to the fasta file of strongly positioned nucleosome coordinates. 
if strongPosNucleosomeFilePath.rsplit('.',1)[0].endswith("_expansion"):
    strongPosNucleosomeFastaFilePath = strongPosNucleosomeFilePath.rsplit("_expansion",1)[0] + ".fa"
else:
    strongPosNucleosomeFastaFilePath = strongPosNucleosomeFilePath.rsplit('.',1)[0] + ".fa"

# Make sure we have a fasta file for strongly positioned nucleosome coordinates
parseStrongPosNucleosomeData(strongPosNucleosomeFilePath,strongPosNucleosomeFastaFilePath)

# Loop through each given mutation background file path, creating a corresponding nucleosome mutation background for each.
for mutationBackgroundFilePath in mutationBackgroundFilePaths:

    mutationBackgroundFileName = mutationBackgroundFilePath.rsplit('/',1)[-1]
    print("Working with file:",mutationBackgroundFileName)

    # Get some information on the file system and generate file paths for convenience.
    workingDirectory = os.path.dirname(mutationBackgroundFilePath) # The working directory for the current data "group"

    # Does it look like we were given a mutation background file?
    if not "_mutation_background" in mutationBackgroundFileName: 
        raise ValueError("Error, expected file with \"_mutation_background\" in the name.")
    # The name of the mutation data set the background pertains to.
    mutationDataGroup = mutationBackgroundFileName.rsplit("_mutation_background",1)[0]

    # A path to the final output file.
    nucleosomeMutationBackgroundFilePath = os.path.join(workingDirectory,mutationDataGroup+"_nucleosome_mutation_background.txt")

    # Generate the nucleosome mutation background file!
    generateNucleosomeMutationBackground(strongPosNucleosomeFastaFilePath,mutationBackgroundFilePath,nucleosomeMutationBackgroundFilePath)