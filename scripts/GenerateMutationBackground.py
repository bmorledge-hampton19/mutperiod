# This script generates a mutational background file, given a mutation file and 
# a genome fasta file.

from UsefulBioinformaticsFunctions import reverseCompliment
from TkinterDialog import TkinterDialog, Selections
import os

# This function generates a file containing the frequencies of each trinucleotide context in a genome
# on the given strand, including N-containing values.
def generateGenomeTrinucFrequencyFile(genomeFilePath, genomeTrinucFrequencyFilePath):

    trinucCounts=dict() # A dictionary of all the relevant trinucleotide and their counts.
    totalTrinucleotides = 0 # The number of nucleotides contained in the genome file.        

    # Get ready to read from the genome file.
    with open(genomeFilePath, 'r') as genomeFile:

        lineLeftovers = '' # Any characters left over from the last line to account for trinucleotide sequences spanning 2 lines.
        willCount = False # A boolean value to decide whether or not to count the sequence in the fasta file.

        # Read the genome file line by line, and count trinuc frequencies as we go.
        for line in genomeFile:

            # Go away, whitespace! >:(
            line = line.strip()

            # Check and see if we're starting a new chromosome with this line.
            if line.startswith(">"):
                sequenceTitle = line.split(">")[1]

                # Check and make sure that the sequence is one we actually want to count.
                if not '_' in sequenceTitle and not "chrM" in sequenceTitle and not "pUC19" in sequenceTitle:

                    print ("Counting trinucleotide sequences in ",sequenceTitle,"...",sep='')
                    lineLeftovers = ''
                    willCount = True

                    # Adjust the count for the lack of a trinucleotide context around the last base on either end of the chromosome.
                    totalTrinucleotides -= 2 

                else:
                    print ("Skipping",sequenceTitle)
                    willCount = False
            
            # Otherwise, count the trinuc sequences.
            elif willCount:

                # Increment the total nucleotide count and add any leftovers to the current line.
                totalTrinucleotides += len(line)
                currentSequence = lineLeftovers + line

                # Count all available trinuc sequences.
                for i in range(0,len(currentSequence)-2):
                    trinuc = currentSequence[i:i+3]
                    
                    # Add the trinuc to the dictionary if it's not there, and then increment its count.
                    trinucCounts.setdefault(trinuc,0)
                    trinucCounts[trinuc] += 1

                # Don't forget about the leftovers for this line!
                lineLeftovers = line[-2:]

    # Open the file to write the counts to.
    with open(genomeTrinucFrequencyFilePath, 'w') as genomeTrinucFrequencyFile:

        genomeTrinucFrequencyFile.write("Total trinucleotide sequences (one strand): " + str(totalTrinucleotides) + '\n')

        # Output headers for the data.
        genomeTrinucFrequencyFile.write('\t'.join(("Trinuc","Occurrences","Frequency")) + '\n')

        # On each line of the file, write the trinucleotide sequence, how many times it occurred,
        # and the frequency with respect to the total number of nucleotides in the genome.
        for trinuc in sorted(trinucCounts.keys()):
            genomeTrinucFrequencyFile.write('\t'.join( (trinuc,str(trinucCounts[trinuc]),str(int(trinucCounts[trinuc])/totalTrinucleotides)) ))
            genomeTrinucFrequencyFile.write('\n')


# This function gets the set of trinucleotide counts for a given genome.
# By default, the counts across both strands are given, but this can be changed to return
# one strand or the other.
def getGenomeTrinucCounts(genomeTrinucFrequencyFilePath, countPlusStrand = True, countMinusStrand = True):

    # Make sure one strand is actually being counted.
    if not (countMinusStrand or countPlusStrand):
        raise ValueError("Error: At least one strand must be selected.")

    # A function that adds the given trinucleotide counts to a given trinuc in the dictionary.
    # Initializes the trinuc if necessary.
    def addFrequency(trinuc, counts):
        trinucleotideCounts.setdefault(trinuc,0)
        trinucleotideCounts[trinuc] += counts

    trinucleotideCounts = dict() # A dictionary to store the trinucleotide counts.

    # Access the file with trinucleotide counts.
    with open(genomeTrinucFrequencyFilePath, 'r') as genomeTrinucFrequencyFile:

        for lineNum,line in enumerate(genomeTrinucFrequencyFile):

            # The first two lines are headers, so ignore them.
            if lineNum < 2: continue

            # Get the trinucleotide sequence, its reverse compliment, and its counts for this line.
            trinuc = line.strip().split('\t')[0]
            reverseTrinuc = reverseCompliment(trinuc)
            counts = int(line.strip().split('\t')[1])
            
            # Add counts to the dictionary based on the parameters set.
            if countPlusStrand and countMinusStrand:
                addFrequency(trinuc,counts)
                addFrequency(reverseTrinuc,counts)

            elif countPlusStrand:
                addFrequency(trinuc,counts)

            elif countMinusStrand:
                addFrequency(reverseTrinuc,counts)

    return trinucleotideCounts


# This function generates a file containing the frequencies of each trinucleotide context that appears in a given mutation file.
def generateMutationTrinucFrequencyFile(mutationFilePath, mutationTrinucFrequencyFilePath):

    trinucCounts = dict() # A dictionary of all relevant nucleotides and their counts.

    # Open the mutation bed file.
    with open(mutationFilePath,'r') as mutationFile:

        # Read through the lines and count trinucleotide sequences.
        for line in mutationFile:
            trinuc = line.strip().split('\t')[3]

            # Did we actually get a trinuc sequence?
            if not len(trinuc) == 3:
                raise ValueError("Error, encountered \"" + trinuc + "\" which is not a trinucleotide sequence.")

            trinucCounts.setdefault(trinuc,0)
            trinucCounts[trinuc] += 1

    # Get the total number of mutations by summing the trinucleotide counts
    totalMutations = sum(trinucCounts.values())

    # Open the file to write the mutations to.
    with open(mutationTrinucFrequencyFilePath, 'w') as mutationTrinucFrequencyFile:

        mutationTrinucFrequencyFile.write("Total Mutations: " + str(totalMutations) + '\n')
        # Output headers for the data.
        mutationTrinucFrequencyFile.write('\t'.join(("Trinuc","Occurrences","Frequency")) + '\n')

        # On each line of the file, write the trinucleotide sequence, how many times it occurred,
        # and the frequency with respect to the total number of mutations.
        for trinuc in sorted(trinucCounts.keys()):
            mutationTrinucFrequencyFile.write('\t'.join( (trinuc,str(trinucCounts[trinuc]),str(int(trinucCounts[trinuc])/totalMutations)) ))
            mutationTrinucFrequencyFile.write('\n')


# This function returns a dictionary with the counts of mutations in each trinucleotide context.
def getMutationTrinucCounts(mutationTrinucFrequencyFilePath):

    trinucleotideCounts = dict() # A dictionary to store the trinucleotide frequencies.

    # Access the file with trinucleotide counts.
    with open(mutationTrinucFrequencyFilePath, 'r') as mutationTrinucFrequencyFile:

        for lineNum,line in enumerate(mutationTrinucFrequencyFile):

            # The first two lines are headers, so ignore them.
            if lineNum < 2: continue

            # Get the trinucleotide sequence for this line and the number of associated mutations.
            trinuc = line.strip().split('\t')[0]
            count = int(line.strip().split('\t')[1])

            # Add the trinucleotide and its frequency to the list of 
            trinucleotideCounts.setdefault(trinuc,count)

    return trinucleotideCounts
            

# This function takes information on the mutational trinuc frequency and the genome trinuc frequency and uses it
# to generate a list of probabilities that a given mutation will arise at a given trinucleotide context.
def generateMutationBackgroundFile(genomeTrinucFrequencyFilePath, mutationTrinucFrequencyFilePath, mutationBackgroundFilePath):
    print("Generating mutation background...")

    # Get data on the necessary data to generate a background mutation rate.
    genomeTrinucCounts = getGenomeTrinucCounts(genomeTrinucFrequencyFilePath)
    mutationTrinucCounts = getMutationTrinucCounts(mutationTrinucFrequencyFilePath)

    backgroundMutationRates = dict() # The rate at which a given trinucleotide sequence is mutated.
    
    # Initialize the dictionary with every trinuc present in the genome.
    for trinuc in genomeTrinucCounts:
        backgroundMutationRates.setdefault(trinuc,0)

    # Populate the dictionary using observed mutations in a trinuc context over the total occurences of that trinucleotide.
    for trinuc in mutationTrinucCounts:
        backgroundMutationRates[trinuc] = mutationTrinucCounts[trinuc]/genomeTrinucCounts[trinuc]

    # Write the results to the mutation background file
    with open(mutationBackgroundFilePath, 'w') as mutationBackgroundFile:

        # Write headers to the file.
        mutationBackgroundFile.write('\t'.join(("Trinuc","Mutation_Rate")) + '\n')

        # Write the data
        for trinuc in sorted(backgroundMutationRates):
            mutationBackgroundFile.write('\t'.join((trinuc,str(backgroundMutationRates[trinuc]))) + '\n')


#Create the Tkinter UI
dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(__file__),"..","data"))
dialog.createFileSelector("Bed Mutation File:",0)
dialog.createFileSelector("Genome fasta File:",1)
dialog.createReturnButton(2,0,2)
dialog.createQuitButton(2,2,2)

# Run the UI
dialog.mainloop()

# If no input was received (i.e. the UI was terminated prematurely), then quit!
if dialog.selections is None: quit()

# Get the user's input from the dialog.
selections: Selections = dialog.selections
mutationFilePath: str = list(selections.getFilePaths())[0] # The path to the bed mutation file
genomeFilePath = list(selections.getFilePaths())[1] # The path to the bed mutation file

# Get some information on the file system and generate file paths for convenience.
workingDirectory = os.path.dirname(mutationFilePath) # The working directory for the current data "group"

# The name of the mutation group used to designate data files (between "/" and "_trinuc").
mutationGroupName = mutationFilePath.rsplit("_trinuc",1)[0].rsplit('/',1)[-1]
# Double check that the mutation group name was generated correctly.
if '.' in mutationGroupName: raise ValueError("Error, expected mutation file with \"trinuc\" in the name.")

# Generate the file path for the genome trinuc frequency file.
genomeTrinucFrequencyFilePath = genomeFilePath.rsplit(".")[0] + "_trinuc_frequency.txt"
# Generate the file path for the mutation trinuc frequency file.
mutationTrinucFrequencyFilePath = os.path.join(workingDirectory,"intermediate_files",mutationGroupName+"_trinuc_frequencies.txt")
# Generate the file path for the background mutation rate file.
mutationBackgroundFilePath = os.path.join(workingDirectory,mutationGroupName + "_mutation_background.txt")

# If the genome trinuc frequency file doesn't exist, create it.
if not os.path.exists(genomeTrinucFrequencyFilePath):
    print("Genome trinuc frequency file not found at path:",genomeTrinucFrequencyFilePath)
    print("Generating genome trinuc frequency file...")
    generateGenomeTrinucFrequencyFile(genomeFilePath, genomeTrinucFrequencyFilePath)

# If the mutation trinuc frequency file doesn't exist, create it.
if not os.path.exists(mutationTrinucFrequencyFilePath):
    print("Mutation trinuc frequency file not found at path:",mutationTrinucFrequencyFilePath)

    # Create a directory for intermediate files if it does not already exist...
    if not os.path.exists(os.path.join(workingDirectory,"intermediate_files")):
        os.mkdir(os.path.join(workingDirectory,"intermediate_files"))

    print("Generating mutation trinuc frequency file...")
    generateMutationTrinucFrequencyFile(mutationFilePath,mutationTrinucFrequencyFilePath)

# Generate the mutation background file.
generateMutationBackgroundFile(genomeTrinucFrequencyFilePath,mutationTrinucFrequencyFilePath,mutationBackgroundFilePath)
print("Done!")





# Sanity check :  Do the plus and minus strands have similar trinucleotide frequencies?

# def printFrequencies(trinucFrequencies):
#     for trinuc in sorted(trinucFrequencies):
#         print(trinuc,": ",trinucFrequencies[trinuc])

# print('\n')
# printFrequencies(getGenomeTrinucFrequency(genomeTrinucFrequencyFilePath,countMinusStrand=False))
# print('\n')
# printFrequencies(getGenomeTrinucFrequency(genomeTrinucFrequencyFilePath,countPlusStrand=False))