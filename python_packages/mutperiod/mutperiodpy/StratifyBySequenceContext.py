# This script stratifies mutperiod data by sequence context, expanding sequence context as necessary.
import os
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.InputParsing.ParseToIterable import parseToIterable
from benbiohelpers.CustomErrors import UserInputError
from mutperiodpy.ExpandContext import expandContext
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (DataTypeStr, getDataDirectory, checkDirs,
                                                                  getContext,
                                                                  Metadata, generateMetadata, generateFilePath)


# A function to determine whether or not a given sequence fits a given pattern at its center.
# e.g. "TCGA" would fit the pattern "TCNN" or "CG" but not "TC"
def doesSequenceFitPattern(sequence, pattern):

    # Make sure the pattern has the same center (half or full base) as the sequence and that it is not longer.
    assert len(sequence) % 2 == len(pattern) % 2, f"{sequence} and {pattern} do not have the same center (half or full base)."
    assert len(pattern) <= len(sequence), f"Pattern: \"{pattern}\" is longer than sequence: \"{sequence}\"."

    # Get the center of the sequence that is equal in length to the pattern.
    offset = int((len(sequence)-len(pattern))/2)
    sequence = sequence[offset:len(sequence)-offset]

    # Remove N's in the pattern and any corresponding characters in the sequence.
    nonN_Indices = [i for i in range(len(pattern)) if pattern[i] != 'N']
    sequence = ''.join([sequence[i] for i in nonN_Indices])
    pattern = ''.join([pattern[i] for i in nonN_Indices])

    # NOTE: If I wanted to make this fully compatible with IUPAC naming conventions,
    #       I could convert IUPAC letters to an re friendly string (i.e. "Y" > "[CT]")

    # Are the resulting sequences the same?
    return(sequence == pattern)


def stratifyBySequenceContext(mutperiodPositionFilePaths: List[str], sequencesToStratifyBy: List[str]):

    # Create a list to store the paths to newly created files.    
    stratifiedMutperiodPositionFilePaths = list()

    # Make sure that all the given sequence patterns have the same parity and record the maximum sequence length.
    parity = None
    maxSequenceLength = 0
    for sequence in sequencesToStratifyBy:

        if parity is None: parity = len(sequence)%2
        elif parity != len(sequence)%2: raise UserInputError("Not all sequences have the same parity.")

        if len(sequence) > maxSequenceLength: maxSequenceLength = len(sequence)
        if maxSequenceLength > 6:
            raise UserInputError(f"Sequences longer than 6 bases are not supported: {sequence}")

    # Make sure that all the given files match sequence parity and have the sufficient context to 
    # match the pattern, expanding files as necessary.
    validPaths = list()
    for mutperiodPositionFilePath in mutperiodPositionFilePaths:

        print(f"\nChecking sequence context for {os.path.basename(mutperiodPositionFilePath)}...")

        fileContext = getContext(mutperiodPositionFilePath, True)

        if fileContext%2 != parity:
            raise UserInputError("File does not have the same parity as sequences to stratify by.")

        if os.path.sep + "sequence_stratifications" + os.path.sep in mutperiodPositionFilePath:
            print("This file has already been stratified by sequence. Removing.")
            continue

        if fileContext < maxSequenceLength:
            print("Found file with insufficient context.  Expanding...")
            validPaths.append(expandContext([mutperiodPositionFilePath], maxSequenceLength)[0])
        else:
            validPaths.append(mutperiodPositionFilePath)
    mutperiodPositionFilePaths = validPaths

    # Iterate through the given files, stratifying by sequence context for each.
    for mutperiodPositionFilePath in mutperiodPositionFilePaths:

        print(f"\nWorking in {os.path.basename(mutperiodPositionFilePath)}...")

        # Generate output directories, paths, and metadata.
        sequenceStratificationsDir = os.path.join(os.path.dirname(mutperiodPositionFilePath),"sequence_stratifications")
        parentMetadata = Metadata(mutperiodPositionFilePath)
        
        for sequence in sequencesToStratifyBy:

            print(f"Stratifying by sequence pattern: {sequence}")

            sequenceDir = os.path.join(sequenceStratificationsDir, sequence)
            checkDirs(sequenceDir)

            dataGroupName = sequence + '_' + parentMetadata.dataGroupName
            sequenceStratifiedFilePath = generateFilePath(sequenceDir, dataGroupName, getContext(mutperiodPositionFilePath),
                                                          dataType = DataTypeStr.mutations, fileExtension = ".bed")
            
            generateMetadata(dataGroupName, parentMetadata.genomeName, 
                             os.path.join("..","..",os.path.basename(mutperiodPositionFilePath)),
                             parentMetadata.inputFormat, sequenceDir, *parentMetadata.cohorts + [sequence],
                             callParamsFilePath = parentMetadata.callParamsFilePath)

            # Check each line in the input file, and output those lines that match the given sequence pattern.
            with open(mutperiodPositionFilePath, 'r') as mutperiodPositionFile:
                with open(sequenceStratifiedFilePath, 'w') as sequenceStratifiedFile:

                    for line in mutperiodPositionFile:
                        if doesSequenceFitPattern(line.split()[3], sequence): sequenceStratifiedFile.write(line)

            stratifiedMutperiodPositionFilePaths.append(sequenceStratifiedFilePath)

    return stratifiedMutperiodPositionFilePaths


def main():
    # Create the Tkinter dialog.
    with TkinterDialog(workingDirectory= getDataDirectory(), title = "Stratify by Sequence Context") as dialog:
        dialog.createMultipleFileSelector("Bed Feature Files:",0,DataTypeStr.mutations + ".bed", ("Bed Files",".bed"))
        dialog.createTextField("Sequence contexts to filter by: ", 1, 0, defaultText = "TCG, NCG")


    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections = dialog.selections

    stratifyBySequenceContext(selections.getFilePathGroups()[0], 
                              parseToIterable(selections.getTextEntries()[0], castValuesToInt=False))

if __name__ == "__main__": main()