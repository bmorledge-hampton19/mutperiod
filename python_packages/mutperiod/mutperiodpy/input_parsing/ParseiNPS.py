# This script takes output from the iNPS nucleosome calling pipeline and adapts it for mutperiod.
# Input data should be stored within the mutperiod __external_data directory under a directory
# with the following structure: __external_data/[genome]/[nucloeomse_map]
import os
from typing import List
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.CustomErrors import InvalidPathError
from benbiohelpers.FileSystemHandling.DirectoryHandling import getIsolatedParentDir
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory

# Given a list of like_bed files from iNPS, convert each to a bed file of nucleosome dyad centers,
# suitable for the rest of the mutperiod pipeline.
# Input data should be stored within the mutperiod __external_data directory under a directory
# with the following structure: __external_data/[genome_name]/[nucloeomse_map_name]
# Returns a list of the output bed files.
def parseiNPS(likeBedFilePaths: List[str]):

    dyadCenterBedFilePaths = list()

    for likeBedFilePath in likeBedFilePaths:

        print(f"Working with {os.path.basename(likeBedFilePath)}")

        if not os.path.dirname(os.path.dirname(os.path.dirname(likeBedFilePath))).endswith("__external_data"):
            raise InvalidPathError(likeBedFilePath, "Expected file within .../__external_data/[genome_name]/"
                                                    "[nucleosome_map_name] directory.")

        # Create the path to the output file.
        dyadCenterBedFilePath = os.path.join(os.path.dirname(likeBedFilePath), getIsolatedParentDir(likeBedFilePath) + ".bed")

        # Convert between the "like_bed" format from iNPS to a bed format suitable for mutperiod.
        print("Converting to mutperiod-friendly bed format...")
        with open(likeBedFilePath, 'r') as likeBedFile:
            with open(dyadCenterBedFilePath, 'w') as dyadCenterBedFile:

                # Skip the headers in the like_bed file.
                likeBedFile.readline(); likeBedFile.readline()

                for line in likeBedFile:
                    splitLine = line.strip().split('\t')
                    inflectionPointsCenter = round((int(splitLine[1]) + int(splitLine[2])) / 2)

                    dyadCenterBedFile.write('\t'.join((splitLine[0], str(inflectionPointsCenter-1), 
                                                       str(inflectionPointsCenter), '.', '.', '.')) + '\n')

        dyadCenterBedFilePaths.append(dyadCenterBedFilePath)

    return dyadCenterBedFilePaths


def main():

    with TkinterDialog(workingDirectory = getDataDirectory(), title = "Parse iNPS Data") as dialog:
        dialog.createMultipleFileSelector("iNPS like_bed Files:", 0, "Gathering.like_bed", ("\"Like Bed\" Files", ".like_bed"))

    parseiNPS(dialog.selections.getFilePathGroups()[0])


if __name__ == "__main__": main()