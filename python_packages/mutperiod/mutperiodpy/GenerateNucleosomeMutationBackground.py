# This script, when given a mutation background file, a genome file, and a file of strongly positioned nucleosome coordinates,
# generates a background file with the expected mutations at each dyad position from -73 to 73 (inclusive).

import os
from typing import Dict
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog, Selections
from mutperiodpy.helper_scripts.UsefulBioinformaticsFunctions import bedToFasta, reverseCompliment, FastaFileIterator
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getContext, getIsolatedParentDir, getLinkerOffset, Metadata, generateFilePath, DataTypeStr, getDataDirectory


# This function takes a bed file of strongly positioned nucleosomes and expands their coordinates to encompass
# the given radius plus 2 bases. (in order to get up to hexanucleotide sequences.)
# If a linker offset is requested, the expansion will be even greater to accomodate.
# The expanded bed file is then used to generate a fasta file of nucleosome sequences.
# Returns the file path to the fasta file.
def generateNucleosomeFasta(baseNucPosFilePath, genomeFilePath, dyadRadius, linkerOffset):

    # Generate a path to the fasta file of nucleosome sequences (Potentially including linker DNA).
    if dyadRadius == 73:
        nucPosFastaFilePath = generateFilePath(directory = os.path.dirname(baseNucPosFilePath),
                                               dataGroup = os.path.basename(baseNucPosFilePath).rsplit('.',1)[0],
                                               linkerOffset = linkerOffset, fileExtension = ".fa") 
    elif dyadRadius == 1000:
        nucPosFastaFilePath = generateFilePath(directory = os.path.dirname(baseNucPosFilePath),
                                               dataGroup = os.path.basename(baseNucPosFilePath).rsplit('.',1)[0],
                                               usesNucGroup = True, fileExtension = ".fa") 
    else: raise ValueError("Invalid counting radius: " + str(dyadRadius))

    # Make sure the file doesn't already exist.  If it does, we're done!
    if os.path.exists(nucPosFastaFilePath):
        print("Found relevant nucleosome fasta file:",os.path.basename(nucPosFastaFilePath))
        return nucPosFastaFilePath
    else: print("Nucleosome fasta file not found at: ",nucPosFastaFilePath,"\nGenerating...", sep = '')

    # Generate the (temporary) expanded file path.
    expandedNucPosBedFilePath = generateFilePath(directory = os.path.dirname(baseNucPosFilePath),
                                                 dataGroup = os.path.basename(baseNucPosFilePath).rsplit('.',1)[0],
                                                 dataType = "expanded", fileExtension = ".bed")

    # Expand the bed coordinates.
    print("Expanding nucleosome coordinates...")
    with open(baseNucPosFilePath,'r') as baseNucPosFile:
        with open(expandedNucPosBedFilePath, 'w') as expandedNucPosBedFile:

            # Write the expanded positions to the new file, one line at a time.
            for line in baseNucPosFile:
                choppedUpLine = line.strip().split('\t')
                choppedUpLine[1] = str(int(choppedUpLine[1]) - dyadRadius - linkerOffset - 2)
                choppedUpLine[2] = str(int(choppedUpLine[2]) + dyadRadius + linkerOffset + 2)

                # Write the results to the expansion file as long as it is not before the start of the chromosome.
                if int(choppedUpLine[1]) > -1: expandedNucPosBedFile.write('\t'.join(choppedUpLine) + '\n')
                else: print("Nucleosome at chromosome", choppedUpLine[0], "with expanded start pos", choppedUpLine[1],
                            "extends into invalid positions.  Skipping.")
                                            
    # Convert the expanded bed file to fasta format.
    print("Converting expanded coordinates to fasta file...")
    bedToFasta(expandedNucPosBedFilePath,genomeFilePath,nucPosFastaFilePath, includeStrand=False)

    return nucPosFastaFilePath


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
def generateDyadPosContextCounts(nucPosFastaFilePath, dyadPosContextCountsFilePath, 
                                 contextNum, dyadRadius, linkerOffset):

    # Dictionary of context counts for every dyad position. (Contains a dictionary of either counts for each dyad position)
    plusStrandNucleosomeDyadPosContextCounts = dict()

    observedContexts = dict() # Hash table of observed contexts for lookup.

    # This is a bit weird.  If the context number is even, we need to account for half positions,
    # but if the context number is odd, we need to keep in mind that there's one extra valid position in the dyad range.
    if contextNum % 2 == 0:
        halfBaseOffset = 0.5
        extraDyadPos = 0
    else:
        halfBaseOffset = 0
        extraDyadPos = 1

    # Initialize the dictionary for context counts on the plus strand.
    for dyadPos in range(-dyadRadius-linkerOffset,dyadRadius+linkerOffset+extraDyadPos): 
        plusStrandNucleosomeDyadPosContextCounts[dyadPos + halfBaseOffset] = dict()

    # Read through the file, adding contexts for every dyad position to the running total in the dictionary
    with open(nucPosFastaFilePath, 'r') as nucPosFastaFile:

        trackedPositionNum = dyadRadius*2 + linkerOffset*2 + extraDyadPos # How many dyad positions we care about.

        for fastaEntry in FastaFileIterator(nucPosFastaFile):

            # Reset dyad position counter
            dyadPos = -dyadRadius - linkerOffset + halfBaseOffset
            
            # Determine how much extra information is present in this line at either end for generating contexts.
            extraContextNum = int((len(fastaEntry.sequence) - trackedPositionNum)/2)

            # Used to pull out the context of desired length.
            extensionLength = contextNum/2 - 0.5

            # Count all available contexts.
            for i in range(0,trackedPositionNum):

                context = fastaEntry.sequence[int(i + halfBaseOffset + extraContextNum - extensionLength):
                                              int(i + halfBaseOffset + extraContextNum + extensionLength+1)]
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
                else: dyadPosContextCountsFile.write("\t0")

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
                dyadPos = float(choppedUpLine[0])
                dyadPosContextCounts[dyadPos] = dict()
                for i,context in enumerate(contexts):
                    dyadPosContextCounts[dyadPos][context] = int(choppedUpLine[i+1])

    return dyadPosContextCounts


# This function generates a nucleosome mutation background file from a general mutation background file
# and a file of strongly positioned nucleosome coordinates.
def generateNucleosomeMutationBackgroundFile(dyadPosContextCountsFilePath, mutationBackgroundFilePath, 
                                             nucleosomeMutationBackgroundFilePath, dyadRadius, linkerOffset):

    # Dictionaries of expected mutations for every dyad position included in the analysis, one for each strand.
    plusStrandNucleosomeMutationBackground = dict() 
    minusStrandNucleosomeMutationBackground = dict()

    # This is a bit weird.  If the context number is even, we need to account for half positions,
    # but if the context number is odd, we need to keep in mind that there's one extra valid position in the dyad range.
    if getContext(mutationBackgroundFilePath, asInt = True) % 2 == 0:
        halfBaseOffset = 0.5
        extraDyadPos = 0
    else:
        halfBaseOffset = 0
        extraDyadPos = 1

    # Initialize the dictionary
    for i in range(-dyadRadius - linkerOffset, dyadRadius + linkerOffset + extraDyadPos): 
        dyadPos = i + halfBaseOffset
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
                             "Expected_Mutations_Minus_Strand","Expected_Mutations_Both_Strands",
                             "Expected_Mutations_Aligned_Strands"))

        nucleosomeMutationBackgroundFile.write(headers + '\n')
        
        # Write the data for each dyad position.
        for i in range(-dyadRadius - linkerOffset, dyadRadius + linkerOffset + extraDyadPos):

            dyadPos = i + halfBaseOffset
            dataRow = '\t'.join((str(dyadPos),str(plusStrandNucleosomeMutationBackground[dyadPos]),
            str(minusStrandNucleosomeMutationBackground[dyadPos]),
            str(plusStrandNucleosomeMutationBackground[dyadPos] + minusStrandNucleosomeMutationBackground[dyadPos]),
            str(plusStrandNucleosomeMutationBackground[dyadPos] + minusStrandNucleosomeMutationBackground[-dyadPos])))

            nucleosomeMutationBackgroundFile.write(dataRow + '\n')


def generateNucleosomeMutationBackground(mutationBackgroundFilePaths, nucleosomeMapNames, useSingleNucRadius, 
                                         useNucGroupRadius, linkerOffset):

    if not (useSingleNucRadius or useNucGroupRadius):
        raise ValueError("Must generate background in either a single nucleosome or group nucleosome radius.")

    nucleosomeMutationBackgroundFilePaths = list() # A list of paths to the output files generated by the function

    # Loop through each given mutation background file path, creating the corresponding nucleosome mutation background(s) for each.
    for mutationBackgroundFilePath in mutationBackgroundFilePaths:

        print("\nWorking with",os.path.basename(mutationBackgroundFilePath))
        if not DataTypeStr.mutBackground in os.path.basename(mutationBackgroundFilePath): 
            raise ValueError("Error, expected file with \"" + DataTypeStr.mutBackground + "\" in the name.")

        for nucleosomeMapName in nucleosomeMapNames:

            print("Counting with nucleosome map:",nucleosomeMapName)

            # Get metadata (Assumes that the metadata has already been generated from a call to countNucleosomePositionMutations)
            metadata = Metadata(os.path.join(os.path.dirname(mutationBackgroundFilePath),nucleosomeMapName))

            # Determine the context of the mutation background file
            contextNum = getContext(mutationBackgroundFilePath, asInt=True)
            contextText = getContext(mutationBackgroundFilePath)
            print("Given mutation background is in", contextText, "context.")

            # To avoid copy pasting this code, here is a simple function to change how the background file is generated 
            # based on the desired dyad radius.
            def generateBackgroundBasedOnRadius(usesNucGroup):

                # Set the dyad radius (And linker offset)
                if usesNucGroup: 
                    dyadRadius = 1000
                    currentLinkerOffset = 0
                else: 
                    dyadRadius = 73
                    currentLinkerOffset = linkerOffset

                # Make sure we have a fasta file for strongly positioned nucleosome coordinates
                nucPosFastaFilePath = generateNucleosomeFasta(metadata.baseNucPosFilePath, metadata.genomeFilePath, dyadRadius, currentLinkerOffset)

                # Generate the path to the tsv file of dyad position context counts
                dyadPosContextCountsFilePath = generateFilePath(directory = os.path.dirname(metadata.baseNucPosFilePath),
                                                                dataGroup = metadata.nucPosName,
                                                                context = contextText, linkerOffset = currentLinkerOffset,
                                                                usesNucGroup = usesNucGroup,
                                                                dataType = "dyad_pos_counts", fileExtension = ".tsv")

                # Make sure we have a tsv file with the appropriate context counts at each dyad position.
                if not os.path.exists(dyadPosContextCountsFilePath): 
                    print("Dyad position " + contextText + " counts file not found at",dyadPosContextCountsFilePath)
                    print("Generating genome wide dyad position " + contextText + " counts file...")
                    generateDyadPosContextCounts(nucPosFastaFilePath, dyadPosContextCountsFilePath,
                                                contextNum, dyadRadius, currentLinkerOffset)

                # A path to the final output file.
                nucleosomeMutationBackgroundFilePath = generateFilePath(directory = metadata.directory, dataGroup = metadata.dataGroupName,
                                                                        context = contextText, linkerOffset = currentLinkerOffset,
                                                                        usesNucGroup = usesNucGroup,
                                                                        dataType = DataTypeStr.nucMutBackground, fileExtension = ".tsv")

                # Generate the nucleosome mutation background file!
                generateNucleosomeMutationBackgroundFile(dyadPosContextCountsFilePath,mutationBackgroundFilePath,
                                                        nucleosomeMutationBackgroundFilePath, dyadRadius, currentLinkerOffset)

                nucleosomeMutationBackgroundFilePaths.append(nucleosomeMutationBackgroundFilePath)

            if useSingleNucRadius:
                generateBackgroundBasedOnRadius(False)
            if useNucGroupRadius:
                generateBackgroundBasedOnRadius(True)

    return nucleosomeMutationBackgroundFilePaths


def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Mutation Background Files:",0,DataTypeStr.mutBackground + ".tsv",("Tab Seperated Values Files",".tsv"))
    dialog.createMultipleFileSelector("Nucleosome Map Files:", 1, "nucleosome_map.bed", ("Bed Files", ".bed"))

    selectSingleNuc = dialog.createDynamicSelector(2,0)
    selectSingleNuc.initCheckboxController("Generate background with a single nucleosome radius (73 bp)")
    linkerSelectionDialog = selectSingleNuc.initDisplay(1, "singleNuc")
    linkerSelectionDialog.createCheckbox("Include 30 bp linker DNA on either side of single nucleosome radius.",0,0)
    selectSingleNuc.initDisplay(0)
    selectSingleNuc.initDisplayState()

    dialog.createCheckbox("Generate background with a nucleosome group radius (1000 bp)", 3, 0)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    mutationBackgroundFilePaths = selections.getFilePathGroups()[0] # A list of mutation file paths
    nucleosomeMapNames = [getIsolatedParentDir(nucleosomeMapFile) for nucleosomeMapFile in selections.getFilePathGroups()[1]]
    if selectSingleNuc.getControllerVar():
        useSingleNucRadius = True
        includeLinker = selections.getToggleStates("singleNuc")[0]
    else:
        useSingleNucRadius = False
        includeLinker = False
    useNucGroupRadius = selections.getToggleStates()[0]

    if includeLinker: linkerOffset = 30
    else: linkerOffset = 0

    generateNucleosomeMutationBackground(mutationBackgroundFilePaths, nucleosomeMapNames, useSingleNucRadius, 
                                         useNucGroupRadius, linkerOffset)

if __name__ == "__main__": main()