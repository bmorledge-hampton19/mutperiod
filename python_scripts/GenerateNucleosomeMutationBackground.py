# This script, when given a mutation background file, a genome file, and a file of strongly positioned nucleosome coordinates,
# generates a background file with the expected mutations at each dyad position from -73 to 73 (inclusive).

from TkinterDialog import TkinterDialog, Selections
import os, subprocess
from typing import Dict
from UsefulBioinformaticsFunctions import bedToFasta, reverseCompliment, FastaFileIterator

# This function takes a bed file of strongly positioned nucleosomes and expands their coordinates to encompass
# 75 bases on either side of the dyad. (in order to get up to pentanucleotide sequences for positions -73 to 73)
# Also, if a linker offset is requested, the expansion will be even greater.
# Returns the file path to the newly expanded file.
def expandNucleosomeCoordinates(strongPosNucleosomeFilePath, linkerOffset):

    # Generate a file path for the expanded nucleosome coordinate file.
    if includeLinker:
        strongPosNucleosomeExpansionFilePath = strongPosNucleosomeFilePath.rsplit('.',1)[0] + "_expansion+linker.bed"
    else:
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
                choppedUpLine[1] = str(int(choppedUpLine[1]) - 75 - linkerOffset)
                choppedUpLine[2] = str(int(choppedUpLine[2]) + 75 + linkerOffset)

                # Write the results to the expansion file as long as it is not before the start of the chromosome.
                if int(choppedUpLine[1]) != -1: strongPosNucleosomeExpansionFile.write('\t'.join(choppedUpLine) + '\n')
                else: print("Nucleosome at chromosome", choppedUpLine[0], "with expanded start pos", choppedUpLine[1],
                            "extends into invalid positions.  Skipping.")

    # Return the file path to the newly expanded file.
    return strongPosNucleosomeExpansionFilePath


# The user could have given a nucleosome positioning file in one of several reasonable formats... Parse away!
def parseStrongPosNucleosomeData(strongPosNucleosomeFilePath,strongPosNucleosomeFastaFilePath,genomeFilePath,linkerOffset):

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
                print("Given unexpanded nucleosome coordinates.")
                strongPosNucleosomeExpandedFilePath = expandNucleosomeCoordinates(strongPosNucleosomeFilePath,includeLinker)

            # If it is expanded, but not with/without linker DNA as specified, expand as necessary.
            elif int(choppedUpLine[2]) - int(choppedUpLine[1]) == 151:
                print("Given expanded nucleosome coordinates without linker DNA.")
                if linkerOffset:
                    print("Need expanded nucleosome coordinates including linker DNA.  Calling expander...")
                    strongPosNucleosomeFilePath = strongPosNucleosomeFilePath.rsplit("_expansion",1)[0] + ".bed"
                    strongPosNucleosomeExpandedFilePath = expandNucleosomeCoordinates(strongPosNucleosomeFilePath, includeLinker)

            elif int(choppedUpLine[2]) - int(choppedUpLine[1]) == 211:
                print("Given expanded nucleosome cooridantes with linker DNA")
                if not linkerOffset:
                    print("Need expanded nucleosome coordinates not including linker DNA.  Calling expander...")
                    strongPosNucleosomeFilePath = strongPosNucleosomeFilePath.rsplit("_expansion",1)[0] + ".bed"
                    strongPosNucleosomeExpandedFilePath = expandNucleosomeCoordinates(strongPosNucleosomeFilePath, includeLinker)
            
            else:
                raise ValueError("Error: Strongly positioned nucleosome data is not in the expected format.\n" +  
                    "Each coordinate should contain the central base pair in the dyad only, " + 
                    "or, exactly 75 bp on either side,\n" +
                    "or exactly 105 bp on either side (30 bp of linker DNA).")

            # Convert the file to fasta format!
            print("Converting strongly positioned nucleosome coordinates to fasta file...")
            bedToFasta(strongPosNucleosomeExpandedFilePath,genomeFilePath,strongPosNucleosomeFastaFilePath, includeStrand=False)


# Determines the context used to generate the given mutation background
def determineBackgroundContext(mutationBackgroundFilePath):

    with open(mutationBackgroundFilePath, 'r') as mutationBackgroundFile:

        # Skip the header line
        mutationBackgroundFile.readline()

        # Return the length of the context
        choppedUpLine = mutationBackgroundFile.readline().strip().split('\t')
        return len(choppedUpLine[0])


# Returns a dictionary of background mutation rates for the relevant contexts in the associated genome.
def getGenomeBackgroundMutationRates(genomeMutationBackgroundFilePath):

    genomeBackgroundMutationRates = dict()

    # Open the file and read the lines into the dictionary.
    with open(genomeMutationBackgroundFilePath, 'r') as genomeMutationBackgroundFile:
        genomeMutationBackgroundFile.readline() # Skip the line with headers.
        for line in genomeMutationBackgroundFile:
            choppedUpLine = line.strip().split('\t')
            genomeBackgroundMutationRates[choppedUpLine[0]] = float(choppedUpLine[1])
    
    return genomeBackgroundMutationRates


# This function generates a file of context counts in the genome for each dyad position.
def generateDyadPosContextCounts(strongPosNucleosomeFastaFilePath, dyadPosContextCountsFilePath, 
                                 contextNum, linkerOffset):

    # Dictionary of context counts for every dyad position. (Contains a dictionary of either counts for each dyad position)
    plusStrandNucleosomeDyadPosContextCounts = dict()

    observedContexts = dict() # Hash table of observed contexts for lookup.

    # Initialize the dictionary for context counts on the plus strand.
    for dyadPos in range(-73-linkerOffset,74+linkerOffset): 
        plusStrandNucleosomeDyadPosContextCounts[dyadPos] = dict()

    # Read through the file, adding contexts for every dyad position to the running total in the dictionary
    with open(strongPosNucleosomeFastaFilePath, 'r') as strongPosNucleosomeFastaFile:

        trackedPositionNum = 73*2 + 1 + linkerOffset*2 # How many dyad positions we care about.

        for fastaEntry in FastaFileIterator(strongPosNucleosomeFastaFile):

            # Reset dyad position counter
            dyadPos = -73 - linkerOffset
            
            # Determine how much extra information is present in this line at either end for generating contexts.
            extraContextNum = len(fastaEntry.sequence) - trackedPositionNum

            # Make sure we have an even number before dividing by 2 (for both ends)
            if extraContextNum%2 != 0:
                raise ValueError(str(extraContextNum) + " should be even.")
            else: extraContextNum = int(extraContextNum/2)


            # Used to pull out the context of desired length.
            extensionLength = int(contextNum/2)

            # Count all available contexts.
            for i in range(0,trackedPositionNum):

                context = fastaEntry.sequence[i+extraContextNum - extensionLength:i + extraContextNum + extensionLength+1]
                if len(context) != contextNum: raise ValueError("Sequence length does not match expected context length.")
                if context not in observedContexts: observedContexts[context] = None              
                plusStrandNucleosomeDyadPosContextCounts[dyadPos][context] = plusStrandNucleosomeDyadPosContextCounts[dyadPos].setdefault(context, 0) + 1

                # Increment the dyadPos counter
                dyadPos += 1

    # Write the context counts for every dyad position on the plus strand in the output file
    with open(dyadPosContextCountsFilePath, 'w') as dyadPosContextCountsFile:

        # Write the header for the file
        dyadPosContextCountsFile.write("Dyad_Pos\t" + '\t'.join(observedContexts.keys()) + '\n')

        # Write the context counts at each dyad position.
        for dyadPos in plusStrandNucleosomeDyadPosContextCounts.keys():
            
            dyadPosContextCountsFile.write(str(dyadPos))

            for context in observedContexts.keys():
                if context in plusStrandNucleosomeDyadPosContextCounts[dyadPos]:
                    dyadPosContextCountsFile.write('\t' + str(plusStrandNucleosomeDyadPosContextCounts[dyadPos][context]))
                else: dyadPosContextCountsFile.write('\t0')

            dyadPosContextCountsFile.write('\n')
        
    
# This function retrieves the context counts for each dyad position in a genome from a given file.
# The data is returned as a dictionary of dictionaries, with the first key being dyad position and the second
# being a context.
def getDyadPosContextCounts(dyadPosContextCountsFilePath):
    with open(dyadPosContextCountsFilePath, 'r') as dyadPosContextCountsFile:
        
        contexts = list()
        headersHaveBeenRead = False
        dyadPosContextCounts= dict()

        for line in dyadPosContextCountsFile:

            if not headersHaveBeenRead:
                headersHaveBeenRead = True
                contexts = line.strip().split('\t')[1:]
    
            else: 
                choppedUpLine = line.strip().split('\t')
                dyadPos = int(choppedUpLine[0])
                dyadPosContextCounts[dyadPos] = dict()
                for i,context in enumerate(contexts):
                    dyadPosContextCounts[dyadPos][context] = int(choppedUpLine[i+1])

    return dyadPosContextCounts


# This function generates a nucleosome mutation background file from a general mutation background file
# and a file of strongly positioned nucleosome coordinates.
def generateNucleosomeMutationBackgroundFile(dyadPosContextCountsFilePath, mutationBackgroundFilePath, 
                                             nucleosomeMutationBackgroundFilePath, linkerOffset):
    print("Generating nucleosome mutation background file...")

    # Dictionaries of expected mutations for every dyad position included in the analysis, one for each strand.
    plusStrandNucleosomeMutationBackground = dict() 
    minusStrandNucleosomeMutationBackground = dict()

    # Initialize the dictionary
    for dyadPos in range(-73-linkerOffset,74+linkerOffset): 
        plusStrandNucleosomeMutationBackground[dyadPos] = 0
        minusStrandNucleosomeMutationBackground[dyadPos] = 0

    # Get the corresponding mutation background and context counts dictionaries.
    backgroundMutationRate = getGenomeBackgroundMutationRates(mutationBackgroundFilePath)
    dyadPosContextCounts = getDyadPosContextCounts(dyadPosContextCountsFilePath)

    # Calculate the expected mutation rates for each dyad position based on the context counts at that position and that context's mutation rate
    for dyadPos in dyadPosContextCounts:

        for context in dyadPosContextCounts[dyadPos]:

            reverseContext = reverseCompliment(context)
            
            # Add the context's mutation rate to the running total in the background dictionaries.
            plusStrandNucleosomeMutationBackground[dyadPos] += backgroundMutationRate[context] * dyadPosContextCounts[dyadPos][context]
            minusStrandNucleosomeMutationBackground[dyadPos] += backgroundMutationRate[reverseContext] * dyadPosContextCounts[dyadPos][context]

    # Write the results of the dictionary to the nucleosome mutation background file.
    with open(nucleosomeMutationBackgroundFilePath, 'w') as nucleosomeMutationBackgroundFile:

        # Write the headers for the data.
        headers = '\t'.join(("Dyad_Position","Expected_Mutations_Plus_Strand",
        "Expected_Mutations_Minus_Strand","Expected_Mutations_Both_Strands"))

        nucleosomeMutationBackgroundFile.write(headers + '\n')
        
        # Write the data for each dyad position.
        for dyadPos in range(-73-linkerOffset,74+linkerOffset):

            dataRow = '\t'.join((str(dyadPos),str(plusStrandNucleosomeMutationBackground[dyadPos]),
            str(minusStrandNucleosomeMutationBackground[dyadPos]),
            str(plusStrandNucleosomeMutationBackground[dyadPos] + minusStrandNucleosomeMutationBackground[dyadPos])))

            nucleosomeMutationBackgroundFile.write(dataRow + '\n')


def generateNucleosomeMutationBackground(mutationBackgroundFilePaths, genomeFilePath, 
                                         strongPosNucleosomeFilePath, includeLinker):
    
    # Set the linker offset.
    if includeLinker: linkerOffset = 30
    else: linkerOffset = 0

    # Whitespace for readability
    print()

    # Generate a path to the fasta file of strongly positioned nucleosome coordinates (Potentially including linker DNA).
    if (strongPosNucleosomeFilePath.rsplit('.',1)[0].endswith("_expansion") or 
        strongPosNucleosomeFilePath.rsplit('.',1)[0].endswith("_expansion+linker")):
        if includeLinker: strongPosNucleosomeFastaFilePath = strongPosNucleosomeFilePath.rsplit("_expansion",1)[0] + "+linker.fa"
        else: strongPosNucleosomeFastaFilePath = strongPosNucleosomeFilePath.rsplit("_expansion",1)[0] + ".fa"
    else:
        if includeLinker: strongPosNucleosomeFastaFilePath = strongPosNucleosomeFilePath.rsplit('.',1)[0] + "+linker.fa"
        else: strongPosNucleosomeFastaFilePath = strongPosNucleosomeFilePath.rsplit('.',1)[0] + ".fa"

    # Make sure we have a fasta file for strongly positioned nucleosome coordinates
    parseStrongPosNucleosomeData(strongPosNucleosomeFilePath,strongPosNucleosomeFastaFilePath,genomeFilePath,linkerOffset)

    nucleosomeMutationBackgroundFilePaths = list() # A list of paths to the output files generated by the function

    # Loop through each given mutation background file path, creating a corresponding nucleosome mutation background for each.
    for mutationBackgroundFilePath in mutationBackgroundFilePaths:

        mutationBackgroundFileName = os.path.split(mutationBackgroundFilePath)[1]
        print("\nWorking with file:",mutationBackgroundFileName)
        if not "_mutation_background" in mutationBackgroundFileName: 
            raise ValueError("Error, expected file with \"_mutation_background\" in the name.")

        # Determine the context of the mutation background file
        backgroundContextNum = determineBackgroundContext(mutationBackgroundFilePath)

        # Set the name of the type of context being used.
        if backgroundContextNum == 1: contextText = "singlenuc"
        elif backgroundContextNum == 3: contextText = "trinuc"
        elif backgroundContextNum == 5: contextText = "pentanuc"

        print("Given mutation background is in", contextText, "context.")

        # Generate the path to the tsv file of dyad position context counts
        if includeLinker:
            dyadPosContextCountsFilePath = strongPosNucleosomeFastaFilePath.rsplit(".",1)[0] + "_dyad_pos+linker"
        else:
            dyadPosContextCountsFilePath = strongPosNucleosomeFastaFilePath.rsplit(".",1)[0] + "_dyad_pos"

            dyadPosContextCountsFilePath += "_" + contextText + "_counts.tsv"

        # Make sure we have a tsv file with the appropriate context counts at each dyad position.
        if not os.path.exists(dyadPosContextCountsFilePath): 
            print("Dyad position " + contextText + " counts file not found at",dyadPosContextCountsFilePath)
            print("Generating genome wide dyad position " + contextText + " counts file...")
            generateDyadPosContextCounts(strongPosNucleosomeFastaFilePath, dyadPosContextCountsFilePath,
                                        backgroundContextNum, linkerOffset)

        # Get some information on the file system and generate file paths for convenience.
        workingDirectory = os.path.dirname(mutationBackgroundFilePath) # The working directory for the current data "group"

        # The name of the mutation data set the background pertains to.
        mutationDataGroup = mutationBackgroundFileName.split("_mutation_background",1)[0]

        # A path to the final output file.
        nucleosomeMutationBackgroundFilePath = os.path.join(workingDirectory,mutationDataGroup+"_"+contextText+"_nucleosome_mutation_background")
        if includeLinker:
            nucleosomeMutationBackgroundFilePath += "+linker.tsv"
        else: 
            nucleosomeMutationBackgroundFilePath += ".tsv"

        # Generate the nucleosome mutation background file!
        generateNucleosomeMutationBackgroundFile(dyadPosContextCountsFilePath,mutationBackgroundFilePath,
                                                 nucleosomeMutationBackgroundFilePath, linkerOffset)

        nucleosomeMutationBackgroundFilePaths.append(nucleosomeMutationBackgroundFilePath)

    return nucleosomeMutationBackgroundFilePaths


if __name__ == "__main__":
    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(__file__),"..","data"))
    dialog.createMultipleFileSelector("Mutation Background Files:",0,"mutation_background.tsv",("Tab Seperated Values Files",".tsv"))
    dialog.createFileSelector("Genome Fasta File:",1,("Fasta Files",".fa"))
    dialog.createFileSelector("Strongly Positioned Nucleosome File:",2,("Bed or Fasta Files",".bed .fa"))
    dialog.createCheckbox("Include linker DNA",3,0,2)
    dialog.createReturnButton(4,0,2)
    dialog.createQuitButton(4,2,2)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    filePaths = list(selections.getFilePaths())
    mutationBackgroundFilePaths = list(selections.getFilePathGroups())[0] # A list of mutation background file paths
    includeLinker: bool = list(selections.getToggleStates())[0] # Whether or not to include linker DNA on either side of the nucleosomes.
    genomeFilePath = filePaths[0] # The path to the genome fasta file
    strongPosNucleosomeFilePath: str = filePaths[1] # The path to the file containing strongly positioned nucleosomes.
  

    generateNucleosomeMutationBackground(mutationBackgroundFilePaths, genomeFilePath, 
                                         strongPosNucleosomeFilePath, includeLinker)