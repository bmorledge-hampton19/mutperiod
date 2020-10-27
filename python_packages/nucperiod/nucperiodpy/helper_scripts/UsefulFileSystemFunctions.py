# This script contains various functions that I think will often be useful when managing filesystems for projects.

import os, datetime
from enum import Enum

# Get the data directory for nucperiod, creating it from user input if necessary.
def getDataDirectory():

    # Check for the text file which should contain the path to the data directory.
    dataDirectoryTextFilePath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data_dir.txt")

    # If it exists, return the directory path within.
    if os.path.exists(dataDirectoryTextFilePath):
        with open(dataDirectoryTextFilePath, 'r') as dataDirectoryTextFile:
            
            dataDirectory = dataDirectoryTextFile.readline().strip()
            
            # Double check to make sure the data directory is still intact.  
            # If it isn't, inform the user, and progress through the function to recreate it.
            if not os.path.exists(dataDirectory):
                print("Data directory not found at expected location: {}".format(dataDirectory))
                print("Please select a new location to create a data directory.")
            else: return dataDirectory

    else:

        # Create a simple dialog to select a new data directory location.
        from nucperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog, Selections
        dialog = TkinterDialog(workingDirectory=os.path.dirname(os.path.dirname((os.path.abspath(__file__)))))
        dialog.createFileSelector("Location to create new data directory:",0,("Fasta Files",".fa"), directory = True)
        dialog.createExitButtons(1,0)

        # Run the UI
        dialog.mainloop()

        # If no input was received (i.e. the UI was terminated prematurely), then quit!
        if dialog.selections is None: quit()

        selections: Selections = dialog.selections
        dataDirectoryDirectory = selections.getIndividualFilePaths()[0]

        # Make sure a valid directory was given.  Then create the new directory (if it doesn't exist already), 
        # write it to the text file, and return it!
        assert os.path.exists(dataDirectoryDirectory), "Given directory does not exist."

        dataDirectory = os.path.join(dataDirectoryDirectory,"nucperiod_data")
        checkDirs(dataDirectory)
        with open(dataDirectoryTextFilePath, 'w') as dataDirectoryTextFile:
            dataDirectoryTextFile.write(dataDirectory + '\n')
        return dataDirectory


# Get the external data directory, creating it if necessary.
def getExternalDataDirectory(): 
    
    externalDataDirectory = os.path.join(getDataDirectory(), "__external_data")
    checkDirs(externalDataDirectory)
    return externalDataDirectory

# The directory containing the R scripts to call the nucperiodR package functionality.
rScriptsDirectory = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "run_nucperiodR")

# Stores information about data type identifier strings
class DataTypeStr:

    mutations = "context_mutations"
    mutBackground = "mutation_background"
    nucMutBackground = "nucleosome_mutation_background"
    rawNucCounts = "raw_nucleosome_mutation_counts"
    normNucCounts = "normalized_nucleosome_mutation_counts"
    customInput = "custom_input"


# Data input format identifiers and their relevant strings.
class InputFormat(Enum):

    ICGC = "ICGC"
    tXRSeq = "tXR-seq"
    customBed = "customBed"


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
def getIsolatedParentDir(filePath: str, isDir = False):

    if not os.path.isdir(filePath) == isDir:
        if isDir: raise ValueError("Expected directory path, but received the file path: " + filePath)
        else: raise ValueError("Expected file path, but received the directory path: " + filePath)

    if isDir: return filePath.rsplit(os.path.sep,1)[-1]
    else: return filePath.rsplit(os.path.sep,2)[-2]


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


# Returns whether or not the given file path uses the nuc-group radius.
def checkForNucGroup(filePath: str):

    # Get the basename of the file path
    fileName = os.path.basename(filePath)

    # Check for the "nuc-group" identifier.
    return "nuc-group" in fileName


# Generates a file path in standardized format based on given information about the file.
# Many parts are optional, but the file must contain a directory, group name, and a file extension.
# Additionally, the context parameter can be given as a string or an int representing a context.
def generateFilePath(directory = None, dataGroup = None, context = None, linkerOffset = 0, 
                     usesNucGroup = False, dataType = None, fileExtension = None):

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

    if linkerOffset != 0 and usesNucGroup: 
        raise ValueError("No linker offset should be present in a nuc group radius file.")

    if linkerOffset != 0: filePathPieces.append(str(linkerOffset) + "linker+")

    if usesNucGroup: filePathPieces.append("nuc-group")

    if dataType is not None: filePathPieces.append(dataType)

    # Stitch together the file path from its components.
    filePath = '_'.join(filePathPieces) + fileExtension
    return filePath


# Generates a .metadata file from the given information.
def generateMetadata(dataGroupName, associatedGenome, associatedNucleosomePositions, 
                     localParentDataPath, inputFormat, metadataDirectory, callParamsFilePath = None,
                     *cohorts):

    # Open up the metadata file.
    with open(os.path.join(metadataDirectory,".metadata"), 'w') as metadataFile:

        # Write the given data
        metadataFile.write("dataGroupName:\t" + dataGroupName + '\n')

        metadataFile.write("associatedGenome:\t" + associatedGenome + '\n')

        metadataFile.write("associatedNucleosomePositions:\t" + associatedNucleosomePositions + '\n')

        metadataFile.write("localParentDataPath:\t" + localParentDataPath + '\n')

        metadataFile.write("inputFormat: " + inputFormat.value + '\n')

        metadataFile.write("dateTime:\t" + str(datetime.datetime.now()).rsplit(':',1)[0] + '\n')

        if callParamsFilePath is not None:
            metadataFile.write("callParamsFilePath:\t" + callParamsFilePath + '\n')

        if len(cohorts) > 0:
            metadataFile.write("cohorts:\t")
            metadataFile.write(', '.join(cohorts) + '\n')
        else:
            metadataFile.write("cohorts:\tNone\n")


# Keeps track of data about a given data group by accessing the metadata file in the same directory
class Metadata:

    # filePath can be a path to any file in the same directory as the desired metadata or the directory itself.
    def __init__(self,filePath):

        # Get the path to the metadata file.
        if os.path.isdir(filePath):
            self.metadataFilePath = os.path.join(filePath,".metadata")
        else:
            self.metadataFilePath = os.path.join(os.path.dirname(filePath),".metadata")

        # Read the metadata file and put its contents into and dictionary, key-value pairs in the file.
        self.metadata = dict()
        with open(self.metadataFilePath, 'r') as metadataFile:
            for line in metadataFile:

                choppedUpLine = str(line).strip().split(maxsplit = 1)

                if not choppedUpLine[0].endswith(':'):
                    raise ValueError("Malformed metadata line: " + line.strip())

                self.metadata[choppedUpLine[0][:-1]] = choppedUpLine[1]

        # Add the metadata directory to the metadata! (So meta!)
        if os.path.isdir(filePath):
            self.metadata["metadataDirectory"] = filePath
        else:
            self.metadata["metadataDirectory"] = os.path.dirname(filePath)

        # Generate quick and easy to use class members to access common metadata!
        self.wrapMetadataInMembers()


    # Search the metadata for a value paired to a given key.
    def getMetadataByKey(self, key, enforceKey = True):

        if not key in self.metadata:
            if enforceKey: raise ValueError("Identifier " + key + " not found in metadata.")
            else: return None
        
        return self.metadata[key]
    

    def wrapMetadataInMembers(self):

        ### Get a variety of common metadata features, quick and easy!

        self.dataGroupName: str = self.getMetadataByKey("dataGroupName")

        self.genomeName: str = self.getMetadataByKey("associatedGenome")

        self.nucPosName: str = self.getMetadataByKey("associatedNucleosomePositions")

        self.localParentDataPath: str = self.getMetadataByKey("localParentDataPath")

        self.inputFormat = InputFormat(self.getMetadataByKey("inputFormat"))

        self.directory: str = self.getMetadataByKey("metadataDirectory")

        self.dateTime: str = self.getMetadataByKey("dateTime")

        self.callParamsFilePath: str = self.getMetadataByKey("callParamsFilePath", False)

        self.cohorts = list()
        if self.getMetadataByKey("cohorts") != "None":
            self.cohorts += self.getMetadataByKey("cohorts").split('\t')
        
        # Check for addable metadata.
        self.mutationCounts: int = self.getMetadataByKey(self.AddableKeys.mutCounts.value, False)

        ### Get file paths for useful metadata associated files.

        self.genomeFilePath = os.path.join(getExternalDataDirectory(),self.genomeName,self.genomeName+".fa")

        self.baseNucPosFilePath = os.path.join(getExternalDataDirectory(), self.genomeName,
                                               self.nucPosName, self.nucPosName+".bed")

        self.parentDataFilePath = os.path.join(self.directory, self.localParentDataPath)

    class AddableKeys(Enum):

        mutCounts = "mutationCounts"

    # Used to add metadata that cannot be generated when other metadata is initially generated.
    # Note:  The "key" parameter should be a member of the "addableKeys" Enum.
    def addMetadata(self, key: Enum, value):

        # Check the validity of the key, and make sure it doesn't already exist in the metadata.
        assert key in self.AddableKeys, "Given key, \"" + key + "\" is not addable."
        assert self.getMetadataByKey(key.value, False) is None, "Metadata already exists for key: " + key.value

        # Append it to the metadata file.
        with open(self.metadataFilePath, 'a') as metadataFile:
            metadataFile.write(key.value + ':\t' + str(value) + '\n')

        # Re-wrap metadata to include this new addition.
        self.wrapMetadataInMembers()