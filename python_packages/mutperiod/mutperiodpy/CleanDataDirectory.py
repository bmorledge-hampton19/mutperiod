# This script cleans up the data directory by deleting any files in "intermediate_files" directories.
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
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
    
    with TkinterDialog(workingDirectory = getDataDirectory()) as dialog:
        with dialog.createDynamicSelector(0, 0) as dirDynSel:
            dirDynSel.initCheckboxController("Clean default mutperiod data directory", True)
            dirDynSel.initDisplay(False, "altDir").createFileSelector("Alternative directory:", 0, directory=True)


    if dirDynSel.getControllerVar(): 
        print(f"Cleaning {getDataDirectory()}...")
        itemsRemoved = cleanDataDirectory()
    else: 
        print(f"Cleaning {dialog.selections.getIndividualFilePaths('altDir')[0]}...")
        itemsRemoved = cleanDataDirectory(dialog.selections.getIndividualFilePaths("altDir")[0])
    print(f"Deleted {itemsRemoved} items within intermediate directories.")

if __name__ == "__main__": main()