# This script extracts an NRL from a given nucleosome map, without the need to create a full project.
import os, subprocess
from typing import List
from benbiohelpers.FileSystemHandling.DirectoryHandling import getIsolatedParentDir
from mutperiodpy.CountNucleosomePositionMutations import countNucleosomePositionMutations
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import rScriptsDirectory, getDataDirectory
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog


def getNRL(nucleosomeMapFilePaths: List[str], regenerate = False):

    NRLs = list()

    for nucMapFilePath in nucleosomeMapFilePaths:
        nucMapRepeatLengthFilePath = nucMapFilePath.rsplit('.',1)[0] + "_repeat_length.txt"
        if not os.path.exists(nucMapRepeatLengthFilePath) or regenerate:
            print(f"Generating nucleosome repeat length file for {os.path.basename(nucMapFilePath)}...")
            
            nucMapSelfCountsFilePath = countNucleosomePositionMutations((nucMapFilePath,), (getIsolatedParentDir(nucMapFilePath),),
                                                                        None, None, None)[0]
            subprocess.run(("Rscript",os.path.join(rScriptsDirectory,"GetNucleosomeRepeatLength.R"),
                            nucMapSelfCountsFilePath, nucMapRepeatLengthFilePath), check = True)
        with open(nucMapRepeatLengthFilePath, 'r') as nucMapRepeatLengthFile:
            NRLs.append(float(nucMapRepeatLengthFile.readline().strip()))

    return NRLs


def main():

    with TkinterDialog(workingDirectory = getDataDirectory()) as dialog:
        dialog.createMultipleFileSelector("Nucleosome Maps:", 0, "nucleosome_map.bed", ("Bed Files", ".bed"))

    getNRL(dialog.selections.getFilePathGroups()[0])


if __name__ == "__main__": main()