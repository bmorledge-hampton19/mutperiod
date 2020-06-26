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


# Generates a file path in standardized format based on given information about the file.
# Many parts are optional, but the file must contain a directory, group name, and a file extension.
# Additionally, the context parameter can be given as a string or an int representing a context.
def generateFilePath(directory = None, dataGroup = None, context = None, 
                     linkerOffset = 0, dataType = None, fileExtension = None):

    # Check to make sure the minimum necessary components are present.
    if directory is None or dataGroup is None or fileExtension is None:
        raise ValueError("Not all necessary file components are present.  " +
                         "The directory, group name, and file extension are all required.")

    # Bring together the components of the file path from the given parameters.
    filePathPieces = list()
    
    filePathPieces.append(os.path.join(directory,dataGroup))

    if context is not None:

        # A dictionary of numbers which matches them to their respective contexts.
        contexts = {1:"singlenuc", 3:"trinuc", 5:"pentanuc"}

        if isinstance(context, str): 
            if context.lower() not in contexts.values():
                raise ValueError(context + " is not a recognized context.")
            filePathPieces.append(context.lower())
        elif isinstance(context, int): 
            if context not in contexts:
                raise ValueError(str(context) + " is not a valid context value.")
            filePathPieces.append(contexts[context])

    if linkerOffset != 0: filePathPieces.append(str(linkerOffset) + "linker+")

    if dataType is not None: filePathPieces.append(dataType)

    # Stitch together the file path from its components.
    filePath = '_'.join(filePathPieces) + fileExtension
    return filePath


