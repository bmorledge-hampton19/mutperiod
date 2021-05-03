# This script takes a mutation file and the coordinates for DNA bases around nucleosomes and calculates how many
# mutations occured at each dyad position for each (and both) strands.
# NOTE:  Both input files must be sorted for this script to run properly. 
#        (Sorted first by chromosome (string) and then by nucleotide position (numeric))

import os, warnings
from typing import List
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog, Selections
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (Metadata, generateFilePath, getDataDirectory,
                                                                  DataTypeStr, getAcceptableChromosomes)


class MutationData:

    def __init__(self, line, acceptableChromosomes):

        # Read in the next line.
        choppedUpLine = line.strip().split()

        # Assign variables
        self.chromosome = choppedUpLine[0] # The chromosome that houses the mutation.
        self.position = float(choppedUpLine[1]) # The position of the mutation in its chromosome. (0 base)
        self.strand = choppedUpLine[5] # Either '+' or '-' depending on which strand houses the mutation.

        # Make sure the mutation is in a valid chromosome.
        if not self.chromosome in acceptableChromosomes:
            raise ValueError(choppedUpLine[0] + " is not a valid chromosome for the mutation trinuc file.")
        

# Contains data on a single nucleosome position obtained by reading the next available line in a given file.
class NucleosomeData:

    def __init__(self, line):
        
        # Read in the next line.
        choppedUpLine = line.strip().split()

        # Assign variables
        self.chromosome = choppedUpLine[0] # The chromosome that houses the nucleosome.
        self.dyadCenter = int(choppedUpLine[1]) # The position of the base pair at the center of the dyad. (0 base)


# Uses the given nucleosome position file and mutation file to count the number of mutations 
# at a specific radius about each nucleosome (73 bp for single nucleosomes, 1000 for a nucleosome group.).  
# Generates a new file to store these results.
# Optionally, some amount of surrounding linker DNA may be included in the analysis as "linker offset"
# NOTE:  It is VITAL that both files are sorted, first by chromosome number and then by starting coordinate.
#        (Sorted first by chromosome (string) and then by nucleotide position (numeric))
#        This code is pretty slick, but it will crash and burn and give you a heap of garbage as output if the inputs aren't sorted.
class CountsFileGenerator():

    def __init__(self, mutationFilePath, nucPosFilePath, nucleosomeMutationCountsFilePath, 
                 dyadRadius, linkerOffset, acceptableChromosomes):

         # Open the mutation and nucleosome positions files to compare against one another.
        self.mutationFile = open(mutationFilePath, 'r')
        self.nucPosFile = open(nucPosFilePath,'r')

        # Store the other arguments passed to the constructor
        self.acceptableChromosomes = acceptableChromosomes
        self.nucleosomeMutationCountsFilePath = nucleosomeMutationCountsFilePath
        self.dyadRadius = dyadRadius
        self.linkerOffset = linkerOffset

        # Dictionaries holding the number of mutations found at each dyad position from -73 to 73 for each strand (including half positions).
        self.minusStrandMutationCounts = dict() 
        self.plusStrandMutationCounts = dict()
        self.intPositions = list()
        self.halfPositions = list()
        for i in range(-dyadRadius - linkerOffset, dyadRadius + linkerOffset + 1):
            self.minusStrandMutationCounts[i] = 0
            self.plusStrandMutationCounts[i] = 0
            self.intPositions.append(i)
            if i < dyadRadius + linkerOffset:
                self.minusStrandMutationCounts[i+0.5] = 0
                self.plusStrandMutationCounts[i+0.5] = 0
                self.halfPositions.append(i+0.5)

        # Keeps track of mutations that matched to a nucleosome to check for overlap.
        self.mutationsInPotentialOverlap: List[MutationData] = list()

        # The mutation and nucleosome currently being investigated.
        self.currentMutation: MutationData = None
        self.currentNucleosome: NucleosomeData = None


    # Reads in the next mutation from the mutation data into currentMutation
    def readNextMutation(self) -> MutationData:

        # Read in the next line.
        nextLine = self.mutationFile.readline()

        # Check if EOF has been reached.
        if len(nextLine) == 0: 
            self.currentMutation = None
        # Otherwise, read in the next mutation.
        else: self.currentMutation = MutationData(nextLine, self.acceptableChromosomes)

    
    # Reads in the next nucleosome from the nucleosome positioning file into currentNucleosome
    def readNextNucleosome(self) -> NucleosomeData:

        # Read in the next line.
        nextLine = self.nucPosFile.readline()

        # Check if EOF has been reached.
        if len(nextLine) == 0: 
            self.currentNucleosome = None
        # Otherwise, read in the next mutation.
        else: self.currentNucleosome = NucleosomeData(nextLine)

        # Check for mutations in overlapping regions between this nucleosome and the last one.
        if self.currentNucleosome is not None: self.checkMutationsInOverlap()


    # Takes a mutation object and nucleosome object which have unequal chromosomes and read through data until they are equal.
    def reconcileChromosomes(self):
        
        chromosomeChanged = False # A simple flag to determine when to inform the user that a new chromosome is being accessed

        # Until the chromosomes are the same for both mutations and nucleosomes, read through the one with the earlier chromosome.
        while (self.currentMutation is not None and self.currentNucleosome is not None and 
               self.currentMutation.chromosome != self.currentNucleosome.chromosome):
            chromosomeChanged = True
            if self.currentMutation.chromosome < self.currentNucleosome.chromosome: self.readNextMutation()
            else: self.readNextNucleosome()

        if chromosomeChanged and self.currentNucleosome is not None and self.currentMutation is not None: 
            print("Counting in",self.currentNucleosome.chromosome)


    # Determines whether or not the given mutation is past the range of the given radius (+ linker).
    def isMutationPastRadius(self):

        if self.currentMutation is None:
            return True
        elif self.currentMutation.position-self.currentNucleosome.dyadCenter > self.dyadRadius + self.linkerOffset:
            return True
        elif not self.currentMutation.chromosome == self.currentNucleosome.chromosome:
            return True
        else: 
            return False


    # A function which checks to see if a given mutation falls within the radius (+ linker) of the current nucleosome and if it does, 
    # adds it to the list of nucleosome mutations.  
    # (Assumes that the given mutation is not beyond the radius due to isMutationPastRadius.)
    def addMutationIfInRadius(self, mutation: MutationData, firstPass):

        if mutation.position >= self.currentNucleosome.dyadCenter - self.dyadRadius - self.linkerOffset:
            if mutation.strand == '+': 
                self.plusStrandMutationCounts[mutation.position - self.currentNucleosome.dyadCenter] += 1
            elif mutation.strand == '-': 
                self.minusStrandMutationCounts[mutation.position - self.currentNucleosome.dyadCenter] += 1
            else:  raise ValueError("Error:  No strand designation found for mutation")
            
            # Add the mutation to the list of mutations to check for overlap if this is the first time it has been seen.
            if firstPass: self.mutationsInPotentialOverlap.append(mutation)


    # Check to see if any previous mutations that were assigned a dyad position are present in the current nucleosome due to overlap.
    def checkMutationsInOverlap(self):

        # First, get rid of any mutations that fall before the beginning of the range for the current nucleosome.
        self.mutationsInPotentialOverlap = [mutation for mutation in self.mutationsInPotentialOverlap 
                                            if mutation.position >= self.currentNucleosome.dyadCenter - self.dyadRadius - 
                                                                    self.linkerOffset and 
                                               mutation.chromosome == self.currentNucleosome.chromosome]

        # Next, check all remaining mutations to see if they are present in the new nucleosome region.
        for mutation in self.mutationsInPotentialOverlap: self.addMutationIfInRadius(mutation, False)
        

    # Count the mutations at each dyad position based on the given nucleosome range.
    def count(self):
        # Get data on the first mutation and nucleosome and reconcile their chromosomes if necessary to start things off.
        # If either the mutation file or nucleosome positioning file is empty, make sure to bypass the check.
        self.readNextMutation()          
        self.readNextNucleosome()
        if self.currentMutation is None or self.currentNucleosome is None:
            warnings.warn("Empty Mutation or Nucleosome Positions file.  Output will most likely be unhelpful.")
        elif self.currentNucleosome.chromosome == self.currentMutation.chromosome: 
            print("Counting in",self.currentNucleosome.chromosome)
        else: self.reconcileChromosomes()

        # The core loop goes through each nucleosome one at a time and checks mutation positions against it until 
        # one exceeds its rightmost position or is on a different chromosome (or mutations are exhausted).  
        # Then, the next nucleosome is checked, then the next, etc. until no nucleosomes remain.
        while self.currentNucleosome is not None:

            # Read mutations until the mutation is past the range of the current nucleosome.
            while not self.isMutationPastRadius():

                # Check and see if we need to add the mutation to our lists.
                self.addMutationIfInRadius(self.currentMutation, True)
                #Get data on the next mutation.
                self.readNextMutation()

            # Read in a new nucleosome.
            self.readNextNucleosome()

            # Reconcile the mutation data and nucleosome data to be sure
            # that they are looking at the same chromosome for the next iteration
            self.reconcileChromosomes()


    def writeResults(self):

        # Write the results to the output file.
        with open(self.nucleosomeMutationCountsFilePath,'w') as nucleosomeMutationCountsFile:
            
            # Write the headers to the file.
            nucleosomeMutationCountsFile.write('\t'.join(("Dyad_Position","Plus_Strand_Counts",
                                                "Minus_Strand_Counts","Both_Strands_Counts",
                                                "Aligned_Strands_Counts")) + '\n')
            
            # Determine whether or not to include just half positions, just int positions, or all positions in the final output.
            dyadPosRange = self.intPositions
            hasHalfPositions = False

            # Are there half position counts?
            for i in self.halfPositions:
                if self.minusStrandMutationCounts[i] > 0 or self.plusStrandMutationCounts[i] > 0:
                    dyadPosRange = self.halfPositions
                    hasHalfPositions = True
                    break

            # If there were half position counts, are there also int position counts?
            if hasHalfPositions:
                for i in self.intPositions:
                    if self.minusStrandMutationCounts[i] > 0 or self.plusStrandMutationCounts[i] > 0:
                        dyadPosRange = self.minusStrandMutationCounts.keys()
                        break


            # Write the data.
            for i in dyadPosRange:
                nucleosomeMutationCountsFile.write('\t'.join((str(i), str(self.plusStrandMutationCounts[i]), 
                                                              str(self.minusStrandMutationCounts[i]), 
                                                              str(self.plusStrandMutationCounts[i] + self.minusStrandMutationCounts[i]),
                                                              str(self.plusStrandMutationCounts[i] + self.minusStrandMutationCounts[-i]))) 
                                                   + '\n')


def countNucleosomePositionMutations(mutationFilePaths, countSingleNuc, countNucGroup, linkerOffset):

    if not (countSingleNuc or countNucGroup):
        raise ValueError("Must count in either a single nucleosome or group nucleosome radius.")

    nucleosomeMutationCountsFilePaths = list() # A list of paths to the output files generated by the function

    # Loop through each given mutation file path, creating a corresponding nucleosome mutation count file for each.
    for mutationFilePath in mutationFilePaths:

        print("\nWorking with",os.path.split(mutationFilePath)[1])

        # Make sure we have the expected file type.
        if not DataTypeStr.mutations in os.path.basename(mutationFilePath): 
            raise ValueError("Mutation file should have \"" + DataTypeStr.mutations + "\" in the name.")

        # Get metadata and use it to generate a path to the nucleosome positions file.
        metadata = Metadata(mutationFilePath)

        # Get the list of acceptable chromosomes
        acceptableChromosomes = getAcceptableChromosomes(metadata.genomeFilePath)

        # Generate the counts file for a single nucleosome region if requested.
        if countSingleNuc:

            # Generate the output file path
            nucleosomeMutationCountsFilePath = generateFilePath(directory = metadata.directory,
                                                                dataGroup = metadata.dataGroupName, linkerOffset = linkerOffset, 
                                                                fileExtension = ".tsv", dataType = DataTypeStr.rawNucCounts)

            # Ready, set, go!
            print("Counting mutations at each nucleosome position in a 73 bp radius +", str(linkerOffset), "bp linker DNA.")
            counter = CountsFileGenerator(mutationFilePath, metadata.baseNucPosFilePath, 
                                          nucleosomeMutationCountsFilePath, 73, linkerOffset, acceptableChromosomes)
            counter.count()
            counter.writeResults()

            nucleosomeMutationCountsFilePaths.append(nucleosomeMutationCountsFilePath)

        # Generate the counts file for a nucleosome group region if requested.
        if countNucGroup:

            # Generate the output file path
            nucleosomeMutationCountsFilePath = generateFilePath(directory = metadata.directory,
                                                                dataGroup = metadata.dataGroupName, usesNucGroup = True,
                                                                fileExtension = ".tsv", dataType = DataTypeStr.rawNucCounts)

            # Ready, set, go!
            print("Counting mutations at each nucleosome position in a 1000 bp radius.")
            counter = CountsFileGenerator(mutationFilePath, metadata.baseNucPosFilePath, 
                                          nucleosomeMutationCountsFilePath, 1000, 0, acceptableChromosomes)
            counter.count()
            counter.writeResults()

            nucleosomeMutationCountsFilePaths.append(nucleosomeMutationCountsFilePath)

    return nucleosomeMutationCountsFilePaths


def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Mutation Files:",0,DataTypeStr.mutations+".bed",("Bed Files",".bed"))
    
    selectSingleNuc = dialog.createDynamicSelector(1,0)
    selectSingleNuc.initCheckboxController("Count with a single nucleosome radius (73 bp)")
    linkerSelectionDialog = selectSingleNuc.initDisplay(1, "singleNuc")
    linkerSelectionDialog.createCheckbox("Include 30 bp linker DNA on either side of single nucleosome radius.",0,0)
    selectSingleNuc.initDisplay(0)
    selectSingleNuc.initDisplayState()

    dialog.createCheckbox("Count with a nucleosome group radius (1000 bp)", 2, 0)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    mutationFilePaths = selections.getFilePathGroups()[0] # A list of mutation file paths
    if selectSingleNuc.getControllerVar():
        countSingleNuc = True
        includeLinker = selections.getToggleStates("singleNuc")[0]
    else:
        countSingleNuc = False
        includeLinker = False
    countNucGroup = selections.getToggleStates()[0]

    if includeLinker: linkerOffset = 30
    else: linkerOffset = 0

    countNucleosomePositionMutations(mutationFilePaths, countSingleNuc, countNucGroup, linkerOffset)

if __name__ == "__main__": main()