# This script generates a mutational background file, given a mutation file and 
# a genome fasta file.

import os
from benbiohelpers.CustomErrors import UserInputError, InvalidPathError
from benbiohelpers.DNA_SequenceHandling import reverseCompliment
from benbiohelpers.FileSystemHandling.FastaFileIterator import FastaFileIterator
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (Metadata, DataTypeStr, generateFilePath, getContext,
                                                                  getDataDirectory, getAcceptableChromosomes)
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog, Selections


# This function generates a file containing the frequencies of each sequence in a given context for a given genome
# on the given strand, including N-containing values (but only for "acceptable chromosomes").
def generateGenomeContextFrequencyFile(genomeFilePath, genomeContextFrequencyFilePath, contextNum, contextText, acceptableChromosomes):

    contextCounts=dict() # A dictionary of all the relevant tri/singlenuc contexts and their counts.

    # Get ready to read from the genome file.
    with open(genomeFilePath, 'r') as genomeFile:

        # Used to pull out the context of desired length.
        extensionLength = contextNum/2 - 0.5
        if contextNum % 2 == 0: halfBaseOffset = 0.5
        else: halfBaseOffset = 0

        # Read the genome file line by line, and count context frequencies as we go.
        for fastaEntry in FastaFileIterator(genomeFile, False):

            # Check and make sure that the sequence is one we actually want to count.
            if fastaEntry.sequenceName in acceptableChromosomes:

                print ("Counting context sequences in ",fastaEntry.sequenceName,"...",sep='')

                # Count all available context sequences.
                for i in range(0,len(fastaEntry.sequence)):

                    if i + halfBaseOffset >= extensionLength and i + halfBaseOffset < len(fastaEntry.sequence) - extensionLength: 
                        context = fastaEntry.sequence[int(i+halfBaseOffset-extensionLength):int(i+halfBaseOffset+extensionLength+1)]
                    
                        # Add the context to the dictionary if it's not there, and then increment its count.
                        contextCounts.setdefault(context,0)
                        contextCounts[context] += 1

            else:
                print ("Skipping",fastaEntry.sequenceName)

    totalContextCounts = sum(contextCounts.values())

    # Open the file to write the counts to.
    with open(genomeContextFrequencyFilePath, 'w') as genomeContextFrequencyFile:

        genomeContextFrequencyFile.write("Total " + contextText + "s counted (one strand): " + str(totalContextCounts) + '\n')

        # Output headers for the data.
        genomeContextFrequencyFile.write('\t'.join((contextText,"Occurrences","Frequency")) + '\n')

        # On each line of the file, write the context sequence, how many times it occurred,
        # and the frequency with respect to the total number of contexts in the genome.
        for context in sorted(contextCounts.keys()):
            genomeContextFrequencyFile.write('\t'.join( (context,str(contextCounts[context]),str(int(contextCounts[context])/totalContextCounts)) ))
            genomeContextFrequencyFile.write('\n')


# This function gets the set of genome context counts from a given file path.
# By default, the counts across both strands are given, but this can be changed to return
# one strand or the other.
def getGenomeContextCounts(genomeContextFrequencyFilePath, countPlusStrand = True, countMinusStrand = True):

    # Make sure one strand is actually being counted.
    if not (countMinusStrand or countPlusStrand):
        raise ValueError("Error: At least one strand must be selected.")

    # A function that adds the given counts to its context in a dictionary.
    # Initializes the context if necessary.
    def addFrequency(context, counts):
        contextCounts.setdefault(context,0)
        contextCounts[context] += counts

    contextCounts = dict() # A dictionary to store the context counts.

    # Access the file with context counts.
    with open(genomeContextFrequencyFilePath, 'r') as genomeContextFrequencyFile:

        for lineNum,line in enumerate(genomeContextFrequencyFile):

            # The first two lines are headers, so ignore them.
            if lineNum < 2: continue

            # Get the context sequence, its reverse compliment, and its counts for this line.
            context = line.strip().split('\t')[0]
            reverseContext = reverseCompliment(context)
            counts = int(line.strip().split('\t')[1])
            
            # Add counts to the dictionary based on the parameters set.
            if countPlusStrand:
                addFrequency(context,counts)

            if countMinusStrand:
                addFrequency(reverseContext,counts)

    return contextCounts


# This function generates a file containing the frequencies of each context that appears in a given mutation file.
def generateMutationContextFrequencyFile(mutationFilePath, mutationContextFrequencyFilePath,
                                         contextNum, contextText, acceptableChromosomes):

    contextCounts = dict() # A dictionary of all relevant contexts and their counts.

    # Open the mutation bed file.
    with open(mutationFilePath,'r') as mutationFile:

        # Used to pull out the context of desired length.
        middleIndex = None
        extensionLength = None

        # Read through the lines and count contexts.
        for line in mutationFile:

            choppedUpLine = line.strip().split('\t')
            surroundingBases = choppedUpLine[3]

            # Preform some checks and initialize some helpful variables if it hasn't been done previously
            if middleIndex is None:

                # Make sure the file has sufficient information to generate the requested context
                if len(surroundingBases) < contextNum:
                    raise UserInputError("The given mutation file does not have enough information to produce a " + 
                                      contextText + " context.")

                middleIndex = len(surroundingBases)/2 - 0.5
                extensionLength = contextNum/2 - 0.5

            # Pull out the context of the desired length.
            context = surroundingBases[int(middleIndex-extensionLength):int(middleIndex+extensionLength+1)]

            # Make sure we didn't encounter an invalid chromosome.
            if choppedUpLine[0] not in acceptableChromosomes:
                raise UserInputError("Encountered " + choppedUpLine[0] + " which is not a valid chromosome for this genome.")

            contextCounts.setdefault(context,0)
            contextCounts[context] += 1

    # Get the total number of mutations by summing the context counts
    totalMutations = sum(contextCounts.values())

    # Open the file to write the mutations to.
    with open(mutationContextFrequencyFilePath, 'w') as mutationContextFrequencyFile:

        mutationContextFrequencyFile.write("Total Mutations: " + str(totalMutations) + '\n')
        # Output headers for the data.
        mutationContextFrequencyFile.write('\t'.join((contextText,"Occurrences","Frequency")) + '\n')

        # On each line of the file, write the context sequence, how many times it occurred,
        # and the frequency with respect to the total number of mutations.
        for context in sorted(contextCounts.keys()):
            mutationContextFrequencyFile.write('\t'.join( (context,str(contextCounts[context]),str(int(contextCounts[context])/totalMutations)) ))
            mutationContextFrequencyFile.write('\n')


# This function returns a dictionary with the counts of mutations for each context.
def getMutationContextCounts(mutationContextFrequencyFilePath):

    contextCounts = dict() # A dictionary to store the context frequencies.

    # Access the file with context counts.
    with open(mutationContextFrequencyFilePath, 'r') as mutationContextFrequencyFile:

        for lineNum,line in enumerate(mutationContextFrequencyFile):

            # The first two lines are headers, so ignore them.
            if lineNum < 2: continue

            # Get the context sequence for this line and the number of associated mutations.
            context = line.strip().split('\t')[0]
            count = int(line.strip().split('\t')[1])

            # Add the context and its frequency to the dictionary of counts
            contextCounts.setdefault(context,count)

    return contextCounts
            

# This function takes information on the mutational context frequency and the genome context frequency and uses it
# to generate a list of probabilities that a given mutation will arise in a given context.
def generateMutationBackgroundFile(genomeContextFrequencyFilePath, mutationContextFrequencyFilePath, 
                                   mutationBackgroundFilePath, contextText):
    print("Generating mutation background...")

    # Get data on the necessary data to generate a background mutation rate.
    genomeContextCounts = getGenomeContextCounts(genomeContextFrequencyFilePath)
    mutationContextCounts = getMutationContextCounts(mutationContextFrequencyFilePath)

    backgroundMutationRates = dict() # The rate at which a given context is mutated.
    
    # Initialize the dictionary with every context present in the genome.
    for context in genomeContextCounts:
        backgroundMutationRates.setdefault(context,0)

    # Populate the dictionary using observed mutations in the given context over the total occurences of that context.
    for context in mutationContextCounts:
        backgroundMutationRates[context] = mutationContextCounts[context]/genomeContextCounts[context]

    # Write the results to the mutation background file
    with open(mutationBackgroundFilePath, 'w') as mutationBackgroundFile:

        # Write headers to the file.
        mutationBackgroundFile.write('\t'.join((contextText,"Mutation_Rate")) + '\n')

        # Write the data
        for context in sorted(backgroundMutationRates):
            mutationBackgroundFile.write('\t'.join((context,str(backgroundMutationRates[context]))) + '\n')


def generateMutationBackground(mutationFilePaths, backgroundContextNum):
    
    mutationBackgroundFilePaths = list() # A list of paths to the output files generated by the function

    # A dictionary for converting context numbers to text.
    contextNumToText = {1:"singlenuc", 2:"dinuc", 3:"trinuc", 4:"quadrunuc", 5:"pentanuc", 6:"hexanuc"}

    for mutationFilePath in mutationFilePaths:

        # Retrieve metadata
        metadata = Metadata(mutationFilePath)
        intermediateFilesDirectory = os.path.join(metadata.directory,"intermediate_files")

        # If necessary, adjust the context for files with even-length features.
        if getContext(mutationFilePath, asInt = True) % 2 == 0:
            thisBackgroundContextNum = backgroundContextNum + 1
        else: thisBackgroundContextNum = backgroundContextNum

        # Set the name of the type of context being used.
        assert thisBackgroundContextNum in contextNumToText, "Unexpected background context number: " + str(thisBackgroundContextNum)
        contextText = contextNumToText[thisBackgroundContextNum]

        # Get the list of acceptable chromosomes
        acceptableChromosomes = getAcceptableChromosomes(metadata.genomeFilePath)

        print("\nWorking in:",os.path.split(mutationFilePath)[1])
        if not DataTypeStr.mutations in os.path.split(mutationFilePath)[1]:
            raise InvalidPathError(mutationFilePath, "Given mutation file does not have \"" + DataTypeStr.mutations + 
                                   "\" in the name.",
                                   postPathMessage = "Are you sure you inputted a file from the mutperiod pipeline?")

        # Generate the file path for the genome context frequency file.
        genomeContextFrequencyFilePath = generateFilePath(directory = os.path.dirname(metadata.genomeFilePath),
                                                          dataGroup = metadata.genomeName, context = contextText,
                                                          dataType = "frequency", fileExtension = ".tsv")

        # Generate the file path for the mutation context frequency file.
        mutationContextFrequencyFilePath = generateFilePath(directory = intermediateFilesDirectory,
                                                            dataGroup = metadata.dataGroupName, context = contextText,
                                                            dataType = "mutation_frequencies", fileExtension = ".tsv")

        # Generate the file path for the background mutation rate file.
        mutationBackgroundFilePath = generateFilePath(directory = metadata.directory, dataGroup = metadata.dataGroupName,
                                                      context = contextText, dataType = DataTypeStr.mutBackground,
                                                      fileExtension = ".tsv")

        # If the genome context frequency file doesn't exist, create it.
        if not os.path.exists(genomeContextFrequencyFilePath):
            print("Genome", contextText, "context frequency file not found at path:",genomeContextFrequencyFilePath)
            print("Generating genome " + contextText + " context frequency file...")
            generateGenomeContextFrequencyFile(metadata.genomeFilePath, genomeContextFrequencyFilePath, thisBackgroundContextNum, contextText)

        # Create a directory for intermediate files if it does not already exist...
        if not os.path.exists(intermediateFilesDirectory):
            os.mkdir(intermediateFilesDirectory)

        # Create the mutation context frequency file.
        print("Generating mutation context frequency file...")
        generateMutationContextFrequencyFile(mutationFilePath,mutationContextFrequencyFilePath, thisBackgroundContextNum, 
                                             contextText, acceptableChromosomes)

        # Generate the mutation background file.
        generateMutationBackgroundFile(genomeContextFrequencyFilePath,mutationContextFrequencyFilePath,mutationBackgroundFilePath, contextText)

        mutationBackgroundFilePaths.append(mutationBackgroundFilePath)

    return mutationBackgroundFilePaths


def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Bed Mutation Files:",0,DataTypeStr.mutations + ".bed",("Bed Files",".bed"))
    dialog.createDropdown("Background Context",1,0,("Singlenuc/Dinuc", "Trinuc/Quadrunuc", "Pentanuc/Hexanuc"))

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    mutationFilePaths = list(selections.getFilePathGroups())[0] # A list of paths to the bed mutation files
    backgroundContext: str = list(selections.getDropdownSelections())[0] # What context should be used to generate the background.

    # Convert background context to int
    if backgroundContext == "Singlenuc/Dinuc":
        backgroundContextNum = 1
    elif backgroundContext == "Trinuc/Quadrunuc":
        backgroundContextNum = 3
    elif backgroundContext == "Pentanuc/Hexanuc":
        backgroundContextNum = 5
    else: raise ValueError("Matching strings is hard.")

    generateMutationBackground(mutationFilePaths, backgroundContextNum)

if __name__ == "__main__": main()