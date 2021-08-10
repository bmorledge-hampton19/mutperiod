# This script takes a nucleosome map and a map of some features that may encompass the nucleosome dyads.
# These two inputs are used to stratify the nucleosome map into only those nucleosomes encompassed by the given features.
# NOTE:  Both input files must be sorted for this script to run properly. 
#        (Sorted first by chromosome (string) and then by nucleotide position (numeric))

import os, sys
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (getDataDirectory, checkDirs)

from benbiohelpers.FileSystemHandling.DirectoryHandling import checkDirs
from benbiohelpers.CountThisInThat.Counter import ThisInThatCounter
from benbiohelpers.CountThisInThat.InputDataStructures import EncompassedDataDefaultStrand, EncompassingDataDefaultStrand, ENCOMPASSED_DATA
from benbiohelpers.CountThisInThat.CounterOutputDataHandler import CounterOutputDataHandler

class NucleosomeStratifier(ThisInThatCounter):

    def setUpOutputDataHandler(self):
        self.outputDataHandler = CounterOutputDataHandler(self.writeIncrementally)
        self.outputDataHandler.addEncompassedFeatureStratifier(outputName = "Nucleosome")
        self.outputDataHandler.addPlaceholderStratifier()
        self.outputDataHandler.createOutputDataWriter(self.outputFilePath, oDSSubs = (None, -1))

    def constructEncompassingFeature(self, line) -> EncompassingDataDefaultStrand:
        return EncompassingDataDefaultStrand(line, self.acceptableChromosomes)

    def constructEncompassedFeature(self, line) -> EncompassedDataDefaultStrand:
        return EncompassedDataDefaultStrand(line, self.acceptableChromosomes)


# Checks position data against a given bed line.  (Assumes strand designations are the same and that the position is a single nucleotide)
def matchesBedEntry(chromosome, startPos, bedEntry: str):

    bedChromosome, bedStartPos = bedEntry.split()[:2]
    return chromosome == bedChromosome and startPos == float(bedStartPos)


# Takes a nucleosome map and files with ranges to check for overlap to stratify by.
# The stratifying features files should each be present within their own nucleosome map directory.
def stratifyNucleosomeMap(nucleosomeMapDir, stratifyingFeaturesMapFilePaths):

    for stratifyingFeaturesMapFilePath in stratifyingFeaturesMapFilePaths:

        assert nucleosomeMapDir != os.path.dirname(stratifyingFeaturesMapFilePath), (
            "Stratifying features filepath is contained within the same directory as the given nucleosome map.")

        print('\n' + "Working in",os.path.basename(stratifyingFeaturesMapFilePath))

        # Get paths to the input and output nucleosome map files.
        originalNucMapFilePath = os.path.join(nucleosomeMapDir, os.path.basename(nucleosomeMapDir) + ".bed")
        stratifiedNucMapDir = os.path.dirname(stratifyingFeaturesMapFilePath)
        stratifiedNucMapFilePath = os.path.join(stratifiedNucMapDir, os.path.basename(stratifiedNucMapDir) + ".bed")
        stratificationConditionsFilePath = os.path.join(stratifiedNucMapDir, "stratification_conditions.txt")

        # Perform the stratification, writing the results as they are determined.
        stratifier = NucleosomeStratifier(originalNucMapFilePath, stratifyingFeaturesMapFilePath, 
                                          stratifiedNucMapFilePath, writeIncrementally = ENCOMPASSED_DATA)
        stratifier.count()

        # Finally, record the conditions of the stratification
        with open(stratificationConditionsFilePath, 'w') as stratificationConditionsFile:
            stratificationConditionsFile.write("Derived from the original nucleosome map: " + os.path.basename(originalNucMapFilePath) + 
                                               " using " + os.path.basename(stratifyingFeaturesMapFilePath) + ".\n")


def parseArgs(args):
    
    # If only the subcommand was given, run the UI.
    if len(sys.argv) == 2: 
        main(); return

    stratifyingFeaturesFilePaths = set()
    for stratifyingFeaturesFilePath in args.stratifyingFeatures:
        assert not os.path.isdir(stratifyingFeaturesFilePath), (
            "Directory was given for stratifying features where file was expected.")
        stratifyingFeaturesFilePaths.add(os.path.abspath(stratifyingFeaturesFilePath))

    if os.path.isfile(args.base_nucleosome_map): 
        baseNucleosomeMap = os.path.dirname(os.path.abspath(args.base_nucleosome_map))
    else: baseNucleosomeMap = os.path.abspath(args.base_nucleosome_map)

    stratifyNucleosomeMap(baseNucleosomeMap, list(stratifyingFeaturesFilePaths))


def main():

    # Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createFileSelector("Original Nucleosome Map Directory:",0, directory = True)    
    dialog.createMultipleFileSelector("Stratifying Feature Ranges:",1,"stratifying_feature_ranges.narrowPeak",
                                      ("Bed Files",".bed"), ("Narrow Peak File", ".narrowPeak"))

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    stratifyNucleosomeMap(dialog.selections.getIndividualFilePaths()[0], dialog.selections.getFilePathGroups()[0])


if __name__ == "__main__": main()