# This script cleans up the data directory by deleting any files in "intermediate_files" directories.
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory
import os

def cleanDataDirectory(directory = getDataDirectory()):
    
    itemsRemoved = 0

    # Iterate through the given directory
    for item in os.listdir(directory):
        path = os.path.join(directory,item)

        # When an intermediate_files directory is encountered, delete all the files within.
        if item == "intermediate_files":
            for itemToDelete in os.listdir(path):
                pathToDelete = os.path.join(path,itemToDelete)
                if os.path.isdir(pathToDelete): os.rmdir(pathToDelete)
                else: os.remove(pathToDelete)
                itemsRemoved += 1

        # Recursively search any directories that are not intermediate directories
        elif os.path.isdir(path):
            itemsRemoved += cleanDataDirectory(path)

    return itemsRemoved



def main():
    print("Cleaning data directory...")
    print("Deleted",cleanDataDirectory(),"items within intermediate directories.")

if __name__ == "__main__": main()