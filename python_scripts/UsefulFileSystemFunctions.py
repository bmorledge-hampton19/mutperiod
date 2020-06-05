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