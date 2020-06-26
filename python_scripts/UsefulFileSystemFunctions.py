# This script contains various functions that I think will often be useful when managing filesystems for projects.
import os

# Recursively searches the given directory for files with the specified ending. Returns a list of the resulting file paths.
def getFilesInDirectory(directory,validEnding, *additionalValidEndings):
    """Recursively searches the given directory(ies) for files of the specified type."""

    filePaths = list() # The files to return.

    # Iterate through the given directory
    for item in os.listdir(directory):
        path = os.path.join(directory,item)

        if os.path.isdir(path):

            filePaths += getFilesInDirectory(path,validEnding, *additionalValidEndings)

        # Send gzipped files to the copyBedData function to be converted
        else:

            if path.endswith(validEnding): filePaths.append(path)

            else:
                for additionalValidEnding in additionalValidEndings:
                    if path.endswith(additionalValidEnding): 
                        filePaths.append(path)
                        break

    return filePaths


# Returns just the name of the first directory above a given path. (e.g. test/file/path.txt would return "file")
def getIsolatedParentDir(filePath: str):
    return filePath.rsplit(os.path.sep,2)[-2]


# Returns the context associated with a given file path as lowercase text (or an int if specified), or none if there is none.
def getContext(filePath: str, asInt = False):

    # A dictionary of contexts which matches them to their respective numbers.
    contexts = {"singlenuc":1, "trinuc":3, "pentanuc":5}

    # Search for each of the contexts in the filename, and return the first (hopefully only) one that is present.
    fileName = os.path.basename(filePath)
    for context in contexts:

        if context in fileName:

            if asInt: return contexts[context]
            else: return context
    
    return None


# Returns the amount of linker DNA associated with the given file path.
def getLinkerDNAAmount(filePath: str):

    # Get the basename of the file path
    fileName = os.path.basename(filePath)

    # If the "linker+" identifier is present, split on it and the preceding underscore to get the amount of linker DNA
    if "linker+" in fileName:
        linkerNum = int(fileName.split("linker+")[0].rsplit('_',1)[-1])
    else:
        linkerNum = 0

    return linkerNum


