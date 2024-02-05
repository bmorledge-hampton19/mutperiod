# This script takes a mutation file and the coordinates for DNA bases around nucleosomes and calculates how many
# mutations occured at each dyad position for each (and both) strands.
# NOTE:  Both input files must be sorted for this script to run properly. 
#        (Sorted first by chromosome (string) and then by nucleotide position (numeric))

import os
from benbiohelpers.CustomErrors import UserInputError, InvalidPathError
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog, Selections
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (Metadata, generateFilePath, generateMetadata, getDataDirectory,
                                                                  DataTypeStr, getAcceptableChromosomes, checkDirs, getIsolatedParentDir)

from benbiohelpers.CountThisInThat.Counter import ThisInThatCounter
from benbiohelpers.CountThisInThat.InputDataStructures import EncompassedDataDefaultStrand, EncompassingDataDefaultStrand, EncompassingData
from benbiohelpers.CountThisInThat.CounterOutputDataHandler import AmbiguityHandling, OutputDataWriter


def getCountDerivatives(outputDataWriter: OutputDataWriter, getHeaders):
    if getHeaders: return ["Both_Strands_Counts", "Aligned_Strands_Counts"]
    else:
        thisPlusCounts = outputDataWriter.outputDataStructure[outputDataWriter.previousKeys[0]][True]
        thisMinusCounts = outputDataWriter.outputDataStructure[outputDataWriter.previousKeys[0]][False]
        oppositeMinusCounts = outputDataWriter.outputDataStructure[-outputDataWriter.previousKeys[0]][False]
        return [str(thisPlusCounts+thisMinusCounts),str(thisPlusCounts+oppositeMinusCounts)]


class MutationsInNucleosomesCounter(ThisInThatCounter):

    def setUpOutputDataHandler(self):
        super().setUpOutputDataHandler()
        self.outputDataHandler.addRelativePositionStratifier(self.currentEncompassingFeature, extraRangeRadius = self.encompassingFeatureExtraRadius,
                                                             outputName = "Dyad_Position")
        self.outputDataHandler.addStrandComparisonStratifier(strandAmbiguityHandling = AmbiguityHandling.tolerate)
        self.outputDataHandler.createOutputDataWriter(self.outputFilePath, customStratifyingNames=(None, {True:"Plus_Strand_Counts", False:"Minus_Strand_Counts"}),
                                                      getCountDerivatives = getCountDerivatives)

    def constructEncompassingFeature(self, line) -> EncompassingDataDefaultStrand:
        return EncompassingDataDefaultStrand(line, self.acceptableChromosomes)


class NucleosomesInNucleosomesCounter(MutationsInNucleosomesCounter):

    def constructEncompassedFeature(self, line) -> EncompassedDataDefaultStrand:
        return EncompassedDataDefaultStrand(line, self.acceptableChromosomes)


class MutationsInStrandedNucleosomesCounter(ThisInThatCounter):

    def setUpOutputDataHandler(self):
        super().setUpOutputDataHandler()
        self.outputDataHandler.addRelativePositionStratifier(self.currentEncompassingFeature, extraRangeRadius = self.encompassingFeatureExtraRadius,
                                                             outputName = "Dyad_Position", strandSpecificPos = True)
        self.outputDataHandler.addStrandComparisonStratifier(strandAmbiguityHandling = AmbiguityHandling.tolerate)
        self.outputDataHandler.createOutputDataWriter(self.outputFilePath, customStratifyingNames=(None, {True:"Plus_Strand_Counts", False:"Minus_Strand_Counts"}),
                                                      getCountDerivatives = getCountDerivatives)

    def constructEncompassingFeature(self, line) -> EncompassingData:
        return EncompassingData(line, self.acceptableChromosomes)


def countNucleosomePositionMutations(mutationFilePaths, nucleosomeMapNames, countSingleNuc, countNucGroup, linkerOffset, useNucStrand = False):

    # Check for the special case where a nucleosome map is being counted against itself to determine the nucleosome repeat length.
    if (len(mutationFilePaths) == 1 and len(nucleosomeMapNames) == 1 and 
        os.path.basename(mutationFilePaths[0]).rsplit('.',1)[0] == nucleosomeMapNames[0]):
        
        nucleosomeMapFilePath = mutationFilePaths[0]
        nucleosomeMapName = nucleosomeMapNames[0]

        print("Counting nucleosome map", nucleosomeMapName, "against itself in a 1000 bp radius.")

        countsFilePath = generateFilePath(directory = os.path.dirname(nucleosomeMapFilePath),
                                          dataGroup = nucleosomeMapName, usesNucGroup = True, 
                                          fileExtension = ".tsv", dataType = "self_" + DataTypeStr.rawNucCounts)
        acceptableChromosomes = getAcceptableChromosomes(os.path.dirname(os.path.dirname(nucleosomeMapFilePath)))

        counter = NucleosomesInNucleosomesCounter(nucleosomeMapFilePath, nucleosomeMapFilePath, countsFilePath, 
                                                  encompassingFeatureExtraRadius=1000, acceptableChromosomes=acceptableChromosomes)
        counter.count()

        return [countsFilePath]

    if not (countSingleNuc or countNucGroup):
        raise UserInputError("Must count in either a single nucleosome or group nucleosome radius.")

    nucleosomeMutationCountsFilePaths = list() # A list of paths to the output files generated by the function
    nucleosomeMapSortingChecked = False # Use this to make sure files are checked for sorting only once.

    # Loop through each given mutation file path, creating a corresponding nucleosome mutation count file for each.
    for mutationFilePath in mutationFilePaths:

        print("\nWorking with",os.path.split(mutationFilePath)[1])

        # Make sure we have the expected file type.
        if not DataTypeStr.mutations in os.path.basename(mutationFilePath): 
            raise InvalidPathError(mutationFilePath, "Given mutation file does not have \"" + DataTypeStr.mutations + 
                                   "\" in the name.",
                                   postPathMessage = "Are you sure you inputted a file from the mutperiod pipeline?")

        for nucleosomeMapName in nucleosomeMapNames:

            print("Counting with nucleosome map:",nucleosomeMapName)

            # Generate the path to the nucleosome-map-specific directory.
            nucleosomeMapDataDirectory = os.path.join(os.path.dirname(mutationFilePath),nucleosomeMapName)
            checkDirs(nucleosomeMapDataDirectory)

            # Check to see if the metadata for this directory has been generated before, and if not, set it up!
            if not os.path.exists(os.path.join(nucleosomeMapDataDirectory,".metadata")):

                print("No metadata found.  Generating...")

                parentMetadata = Metadata(mutationFilePath)

                # Check to see if the data name should be altered by this nucleosome map.
                dataGroupName = parentMetadata.dataGroupName

                dataGroupNameSuffixFilePath = os.path.join(os.path.dirname(parentMetadata.genomeFilePath), 
                                                           nucleosomeMapName, "append_to_data_name.txt")
                if os.path.exists(dataGroupNameSuffixFilePath):

                    with open(dataGroupNameSuffixFilePath) as dataGroupNameSuffixFile:
                        dataGroupName += dataGroupNameSuffixFile.readline().strip()

                generateMetadata(dataGroupName, parentMetadata.genomeName, os.path.join("..",parentMetadata.localParentDataPath),
                                 parentMetadata.inputFormat, nucleosomeMapDataDirectory, *parentMetadata.cohorts,
                                 callParamsFilePath = parentMetadata.callParamsFilePath,
                                 associatedNucleosomePositions = nucleosomeMapName)

            # Get metadata and use it to generate a path to the nucleosome positions file.
            metadata = Metadata(nucleosomeMapDataDirectory)

            # Get the list of acceptable chromosomes
            acceptableChromosomes = getAcceptableChromosomes(metadata.genomeFilePath)

            # Determine which counter class to use.
            if useNucStrand: CounterClass = MutationsInStrandedNucleosomesCounter
            else: CounterClass = MutationsInNucleosomesCounter

            # Generate the counts file for a single nucleosome region if requested.
            if countSingleNuc:

                # Generate the output file path
                nucleosomeMutationCountsFilePath = generateFilePath(directory = metadata.directory,
                                                                    dataGroup = metadata.dataGroupName, linkerOffset = linkerOffset, 
                                                                    fileExtension = ".tsv", dataType = DataTypeStr.rawNucCounts)

                # Ready, set, go!
                print("Counting mutations at each nucleosome position in a 73 bp radius +", str(linkerOffset), "bp linker DNA.")
                counter = CounterClass(mutationFilePath, metadata.baseNucPosFilePath, nucleosomeMutationCountsFilePath, 
                                       encompassingFeatureExtraRadius=73 + linkerOffset, acceptableChromosomes=acceptableChromosomes,
                                       checkForSortedFiles = (True, not nucleosomeMapSortingChecked))
                counter.count()

                nucleosomeMutationCountsFilePaths.append(nucleosomeMutationCountsFilePath)

            # Generate the counts file for a nucleosome group region if requested.
            if countNucGroup:

                # Generate the output file path
                nucleosomeMutationCountsFilePath = generateFilePath(directory = metadata.directory,
                                                                    dataGroup = metadata.dataGroupName, usesNucGroup = True,
                                                                    fileExtension = ".tsv", dataType = DataTypeStr.rawNucCounts)

                # Ready, set, go!
                print("Counting mutations at each nucleosome position in a 1000 bp radius.")
                counter = CounterClass(mutationFilePath, metadata.baseNucPosFilePath, nucleosomeMutationCountsFilePath, 
                                       encompassingFeatureExtraRadius=1000, acceptableChromosomes=acceptableChromosomes,
                                       checkForSortedFiles = (True, not nucleosomeMapSortingChecked))
                counter.count()

                nucleosomeMutationCountsFilePaths.append(nucleosomeMutationCountsFilePath)

        nucleosomeMapSortingChecked = True

    return nucleosomeMutationCountsFilePaths


def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory(), title = "Count Nucleosome Position Mutations")
    dialog.createMultipleFileSelector("Mutation Files:",0,DataTypeStr.mutations+".bed",("Bed Files",".bed"))
    dialog.createMultipleFileSelector("Nucleosome Map Files:", 1, "nucleosome_map.bed", ("Bed Files", ".bed"))
    selectSingleNuc = dialog.createDynamicSelector(2,0)
    selectSingleNuc.initCheckboxController("Count with a single nucleosome radius (73 bp)")
    linkerSelectionDialog = selectSingleNuc.initDisplay(1, "singleNuc")
    linkerSelectionDialog.createCheckbox("Include 30 bp linker DNA on either side of single nucleosome radius.",0,0)
    selectSingleNuc.initDisplay(0)
    selectSingleNuc.initDisplayState()
    dialog.createCheckbox("Count with a nucleosome group radius (1000 bp)", 3, 0)
    dialog.createCheckbox("Use strand designation in \"nucleosomes\" file", 4, 0)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    mutationFilePaths = selections.getFilePathGroups()[0] # A list of mutation file paths
    nucleosomeMapNames = [getIsolatedParentDir(nucleosomeMapFile) for nucleosomeMapFile in selections.getFilePathGroups()[1]]
    if selectSingleNuc.getControllerVar():
        countSingleNuc = True
        includeLinker = selections.getToggleStates("singleNuc")[0]
    else:
        countSingleNuc = False
        includeLinker = False
    countNucGroup = selections.getToggleStates()[0]
    useNucStrand = selections.getToggleStates()[1]

    if includeLinker: linkerOffset = 30
    else: linkerOffset = 0

    countNucleosomePositionMutations(mutationFilePaths, nucleosomeMapNames, countSingleNuc, countNucGroup, linkerOffset, useNucStrand)

if __name__ == "__main__": main()