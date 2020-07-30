# This script contains various functions that I think will often be useful when managing filesystems for projects.

import os

# The directory containing genome and nucleosome data.
dataDirectory = os.path.join(os.path.dirname(__file__),"..","..","..","..","data")
externalDataDirectory = os.path.join(os.path.dirname(__file__),"..","..","..","..","data","__external_data")
RPackagesDirectory = os.path.join(os.path.dirname(__file__),"..","..","..","..","R_packages")


# Stores information about data type identifiers
class DataTypes:

    def __init__(self):

        self.mutations = "context_mutations"
        self.mutBackground = "mutation_background"
        self.nucMutBackground = "nucleosome_mutation_background"
        self.rawNucCounts = "raw_nucleosome_mutation_counts"
        self.normNucCounts = "normalized_nucleosome_mutation_counts"

# The data type identifiers
dataTypes = DataTypes()


# Recursively searches the given directory for files with the specified ending. Returns a list of the resulting file paths.
def getFilesInDirectory(directory,validEnding, *additionalValidEndings):
    """Recursively searches the given directory(ies) for files of the specified type."""

    filePaths = list()

    # Iterate through the given directory
    for item in os.listdir(directory):
        path = os.path.join(directory,item)

        # Recursively search any directories
        if os.path.isdir(path):
            filePaths += getFilesInDirectory(path,validEnding, *additionalValidEndings)

        # Check files for the valid ending(s)
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


# Checks to see if the given directories exist and creates them if they do not.
def checkDirs(*directoryPaths):
    for directoryPath in directoryPaths:
        if not os.path.exists(directoryPath): os.makedirs(directoryPath)


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
def getLinkerOffset(filePath: str):

    # Get the basename of the file path
    fileName = os.path.basename(filePath)

    # If the "linker+" identifier is present, split on it and the preceding underscore to get the amount of linker DNA
    if "linker+" in fileName:
        linkerOffset = int(fileName.split("linker+")[0].rsplit('_',1)[-1])
    else:
        linkerOffset = 0

    return linkerOffset


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


# Generates a .metadata file from the given information.
def generateMetadata(dataGroupName, associatedGenome, associatedNucleosomePositions, 
                     localParentDataPath, metadataDirectory):

    # Open up the metadata file.
    with open(os.path.join(metadataDirectory,".metadata"), 'w') as metadataFile:

        # Write the given data
        metadataFile.write("dataGroupName:\t" + dataGroupName + '\n')

        metadataFile.write("associatedGenome:\t" + associatedGenome + '\n')

        metadataFile.write("associatedNucleosomePositions:\t" + associatedNucleosomePositions + '\n')

        metadataFile.write("localParentDataPath:\t" + localParentDataPath + '\n')

# Keeps track of data about a given data group by accessing the metadata file in the same directory
class Metadata:

    def __init__(self,filePath):

        # Get the path to the metadata file.
        metadataFilePath = os.path.join(os.path.dirname(filePath),".metadata")

        # Read the metadata file and put its contents into and dictionary, key-value pairs in the file.
        self.metadata = dict()
        with open(metadataFilePath, 'r') as metadataFile:
            for line in metadataFile:

                choppedUpLine = line.strip().split('\t')

                if len(choppedUpLine) != 2 or not choppedUpLine[0].endswith(':'):
                    raise ValueError("Malformed metadata line: " + line.strip())

                self.metadata[choppedUpLine[0][:-1]] = choppedUpLine[1]

        # Add the metadata directory to the metadata! (So meta!)
        self.metadata["metadataDirectory"] = os.path.dirname(filePath)

        # Generate quick and easy to use class members to access common metadata!
        self.wrapMetadataInMembers()


    # Search the metadata for a value paired to a given key.
    def getMetadataByKey(self, key):

        if not key in self.metadata:
            raise ValueError("Identifier " + key + " not found in metadata.")
        
        return self.metadata[key]
    

    def wrapMetadataInMembers(self):

        ### Get a variety of common metadata features, quick and easy!

        self.dataGroupName = self.getMetadataByKey("dataGroupName")

        self.genomeName = self.getMetadataByKey("associatedGenome")

        self.nucPosName = self.getMetadataByKey("associatedNucleosomePositions")

        self.localParentDataPath = self.getMetadataByKey("localParentDataPath")

        self.directory = self.getMetadataByKey("metadataDirectory")

        
        ### Get file paths for useful metadata associated files.

        self.genomeFilePath = os.path.join(externalDataDirectory,self.genomeName,self.genomeName+".fa")

        self.baseNucPosFilePath = os.path.join(externalDataDirectory, self.genomeName,
                                               self.nucPosName, self.nucPosName+".bed")

        self.parentDataFilePath = os.path.join(self.directory, self.localParentDataPath)