# This script contains various functions that I think will often be useful when managing filesystems for projects.

import os, datetime, subprocess
from enum import Enum
from benbiohelpers.FileSystemHandling.DirectoryHandling import checkDirs, getIsolatedParentDir
from mutperiodpy.helper_scripts.CustomErrors import InputError, InvalidPathError

# Get the data directory for mutperiod, creating it from user input if necessary.
def getDataDirectory():

    # Check for the text file which should contain the path to the data directory.
    dataDirectoryTextFilePath = os.path.join(os.getenv("HOME"), ".mutperiod", "data_dir.txt")

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
        from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog, Selections
        checkDirs(os.path.dirname(dataDirectoryTextFilePath))
        dialog = TkinterDialog(workingDirectory = os.path.dirname(dataDirectoryTextFilePath))
        dialog.createFileSelector("Location to create new data directory:",0,("Fasta Files",".fa"), directory = True)

        # Run the UI
        dialog.mainloop()

        # If no input was received (i.e. the UI was terminated prematurely), then quit!
        if dialog.selections is None: quit()

        selections: Selections = dialog.selections
        dataDirectoryDirectory = selections.getIndividualFilePaths()[0]

        # Make sure a valid directory was given.  Then create the new directory (if it doesn't exist already), 
        # write it to the text file, and return it!  (Also create the __external_data directory.)
        if not os.path.exists(dataDirectoryDirectory): 
            raise InputError("Given directory: " + dataDirectoryDirectory + " does not exist.")

        dataDirectory = os.path.join(dataDirectoryDirectory,"mutperiod_data")
        checkDirs(dataDirectory)
        checkDirs(os.path.join(dataDirectory,"__external_data"))
        with open(dataDirectoryTextFilePath, 'w') as dataDirectoryTextFile:
            dataDirectoryTextFile.write(dataDirectory + '\n')
        return dataDirectory


# Get the external data directory, creating it if necessary.
def getExternalDataDirectory(): 
    
    externalDataDirectory = os.path.join(getDataDirectory(), "__external_data")
    checkDirs(externalDataDirectory)
    return externalDataDirectory

# The directory containing the R scripts to call the mutperiodR package functionality.
rScriptsDirectory = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "run_mutperiodR")

# Stores information about data type identifier strings
class DataTypeStr:

    mutations = "context_mutations"
    mutBackground = "mutation_background"
    nucMutBackground = "nucleosome_mutation_background"
    customBackgroundInfo = "custom_background_info"
    rawNucCounts = "raw_nucleosome_mutation_counts"
    normNucCounts = "normalized_nucleosome_mutation_counts"
    generalNucCounts = "nucleosome_mutation_counts"
    customInput = "custom_input"


# Data input format identifiers and their relevant strings.
class InputFormat(Enum):

    ICGC = "ICGC"
    tXRSeq_DEPRECATED = "tXR-seq"
    xRSeq = "XR-seq"
    UVDESeq_DEPRECATED = "UVDE-seq"
    standardBed = "standardBed"
    customBed = "customBed"
    prepared = "prepared"


# Given a genome fasta file (or directory containing a genome fasta file), return the chromosomes present in that file.  
def getAcceptableChromosomes(genomeFilePath: str, returnFilePathInstead = False):

    # If given a directory, retrieve the genome file path from it.
    if os.path.isdir(genomeFilePath): genomeFilePath = os.path.join(genomeFilePath,os.path.basename(genomeFilePath) + ".fa")

    # Make sure we were given a reasonable file path.
    if not getIsolatedParentDir(genomeFilePath) in genomeFilePath and genomeFilePath.endswith(".fa"): 
        raise InvalidPathError(
            "Given file path does not appear to be a genome file path internal to mutperiod.  Make sure to "
            "choose a .fa file within the mutperiod_data/__external_data/[genome_name] directory or choose the genome directory itself."
            )

    # Parse the path to the acceptable cohromosomes file from the given path.
    acceptableChromosomesFilePath = genomeFilePath.rsplit(".fa",1)[0] + "_acceptable_chromosomes.txt"

    # If the acceptable chromosomes file has not been generated, do so.
    if not os.path.exists(acceptableChromosomesFilePath):
        print("Acceptable chromosomes file not found at expected location.  Generating from given fasta file...")

        with open(genomeFilePath, 'r') as genomeFile:
            with open(acceptableChromosomesFilePath, 'w') as acceptableChromosomesFile:

                for line in genomeFile:
                    if line.startswith('>'):
                        chromosomeName = line[1:].split(maxsplit = 1)[0]
                        print("Found chromosome:", chromosomeName)
                        acceptableChromosomesFile.write(chromosomeName + '\n')
        
        print("If these chromosome designations seem incorrect, check that the genome fasta file headers are formatted correctly.  "
              "Chromosome names are defined as the string before the first whitespace character and after the '>' in each header line.")

    # Create a list of acceptable chromosome strings from the acceptable chromosomes file and return it.
    # Or, if requested, return the path to the acceptable chromosomes file instead.
    if returnFilePathInstead: return acceptableChromosomesFilePath
    else:
        with open(acceptableChromosomesFilePath, 'r') as acceptableChromosomesFile:
            return [line.strip() for line in acceptableChromosomesFile]


# Returns the context associated with a given file path as lowercase text (or an int if specified), or none if there is none.
def getContext(filePath: str, asInt = False):

    # A dictionary of contexts which matches them to their respective numbers.
    contexts = {"singlenuc":1, "dinuc":2, "trinuc":3, "quadrunuc":4, "pentanuc":5, "hexanuc":6, 
                "polynuc":float('inf'), "custom_context":-1, "mixed_context":0}

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


# Returns the number of mutations that fall within nucleosomes (positions -73 to 73)
# Input can be either a directory or the actual raw nucleosome counts file.
# Only runs on raw nucleosome mutation counts files.
# Returns None if no such file was found. 
# Only counts positions where the absolute value of their dyad position is less than or equal to
# the given dyadPosCutoff.
def getNucMutCounts(rawNucCountsPath, dyadPosCutoff = 60):

    # Set default values
    rawNucCountsFilePath = None
    nucMutCounts = None

    # If given a directory, search for raw counts file path.
    if os.path.isdir(rawNucCountsPath):
        for item in os.listdir(rawNucCountsPath):
            if DataTypeStr.rawNucCounts in item:
                rawNucCountsFilePath = os.path.join(rawNucCountsPath, item)
                
    
    # Otherwise, double check that the given file is a raw counts file.
    elif DataTypeStr.rawNucCounts in os.path.basename(rawNucCountsPath):
        rawNucCountsFilePath = rawNucCountsPath

    # if a raw counts file was found, determine the total number of mutations in the given region from it
    # (Within the given dyadPosCutoff range).
    if rawNucCountsFilePath is not None:
        nucMutCounts = 0
        with open(rawNucCountsFilePath, 'r') as rawNucCountsFile:
            rawNucCountsFile.readline() # Skip header line
            for line in rawNucCountsFile:             
                if abs(int(line.strip().split('\t')[0])) <= dyadPosCutoff:
                    nucMutCounts += int(line.strip().split('\t')[3])

    return nucMutCounts


# Returns whether or not the given file path uses the nuc-group radius.
def checkForNucGroup(filePath: str):

    # Get the basename of the file path
    fileName = os.path.basename(filePath)

    # Check for the "nuc-group" identifier.
    return "nuc-group" in fileName


# Retrieve the expected periods for each of the given counts files, generating them as necessary.
def getExpectedPeriod(nucleosomeMutationCountsFilePath):

    # If this file uses a nuc-group radius (1000bp) the expected period must be determined empirically from the nucleosome map
    if checkForNucGroup(nucleosomeMutationCountsFilePath):

        nucMapFilePath = Metadata(nucleosomeMutationCountsFilePath).baseNucPosFilePath
        nucMapRepeatLengthFilePath = nucMapFilePath.rsplit('.',1)[0] + "_repeat_length.txt"

        # If the file containing the nucleosome repeat length has not been generated, generate it!
        if not os.path.exists(nucMapRepeatLengthFilePath):

            # Import the counting function here.  (Importing at the top of the script creates a circular reference)
            from mutperiodpy.CountNucleosomePositionMutations import countNucleosomePositionMutations

            print("No repeat length file found for nucleosome map ",os.path.basename(nucMapFilePath),".  Generating...", sep = '')
            nucMapSelfCountsFilePath = countNucleosomePositionMutations((nucMapFilePath,), (getIsolatedParentDir(nucMapFilePath),),
                                                                        None, None, None)[0]
            subprocess.run(("Rscript",os.path.join(rScriptsDirectory,"GetNucleosomeRepeatLength.R"),
                            nucMapSelfCountsFilePath, nucMapRepeatLengthFilePath), check = True)

        # Retrieve the repeat length for the nucleosome map.
        with open(nucMapRepeatLengthFilePath, 'r') as nucMapRepeatLengthFile:
            return float(nucMapRepeatLengthFile.readline().strip())
            
    # Otherwise, this file uses a single nucleosome radius, and the expected period is simply 10.2
    else: return 10.2


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
        contexts = {1:"singlenuc", 2:"dinuc", 3:"trinuc", 4:"quadrunuc", 5:"pentanuc", 6:"hexanuc", 
                    float('inf'):"polynuc", -1:"custom_context", 0:"mixed_nucleotide"}

        if isinstance(context, str): 
            if context.lower() not in contexts.values():
                raise ValueError(context + " is not a recognized context.")
            filePathPieces.append(context.lower())
        elif isinstance(context, int) or context == float('inf'): 
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
def generateMetadata(dataGroupName, associatedGenome, localParentDataPath, inputFormat, metadataDirectory, *cohorts,
                     associatedNucleosomePositions = None, callParamsFilePath = None):

    # Open up the metadata file.
    with open(os.path.join(metadataDirectory,".metadata"), 'w') as metadataFile:

        # Write the given data
        metadataFile.write("dataGroupName:\t" + dataGroupName + '\n')

        metadataFile.write("associatedGenome:\t" + associatedGenome + '\n')

        if associatedNucleosomePositions is not None:
            metadataFile.write("associatedNucleosomePositions:\t" + associatedNucleosomePositions + '\n')

        metadataFile.write("localParentDataPath:\t" + localParentDataPath + '\n')

        metadataFile.write("inputFormat: " + inputFormat.value + '\n')

        metadataFile.write("dateTime:\t" + str(datetime.datetime.now()).rsplit(':',1)[0] + '\n')

        if callParamsFilePath is not None:
            metadataFile.write("callParamsFilePath:\t" + callParamsFilePath + '\n')

        if len(cohorts) > 0:
            metadataFile.write("cohorts:\t")
            metadataFile.write(", ".join(cohorts) + '\n')
        else:
            metadataFile.write("cohorts:\tNone\n")


# Keeps track of data about a given data group by accessing the metadata file in the same directory
class Metadata:

    # filePath can be a path to any file in the same directory as the desired metadata or the directory itself.
    def __init__(self,filePath):

        # Make sure the given file path exists. (Otherwise, the isdir call below can be EXTREMELY misleading)
        if not os.path.exists(filePath): raise InvalidPathError("Given path to retrieve metadata from does not exist: " + filePath)

        # Get the path to the metadata file.
        if os.path.isdir(filePath):
            self.metadataFilePath = os.path.join(filePath,".metadata")
        else:
            self.metadataFilePath = os.path.join(os.path.dirname(filePath),".metadata")

        # Make sure the generated metadata file path exists
        if not os.path.exists(filePath): raise InvalidPathError("Metadata not found at expected location: " + self.metadataFilePath)

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

        self.nucPosName: str = self.getMetadataByKey("associatedNucleosomePositions", False)

        self.localParentDataPath: str = self.getMetadataByKey("localParentDataPath")

        self.inputFormat = InputFormat(self.getMetadataByKey("inputFormat"))

        self.directory: str = self.getMetadataByKey("metadataDirectory")

        self.dateTime: str = self.getMetadataByKey("dateTime")

        self.callParamsFilePath: str = self.getMetadataByKey("callParamsFilePath", False)

        self.cohorts = list()
        if self.getMetadataByKey("cohorts") != "None":
            self.cohorts += self.getMetadataByKey("cohorts").split(", ")
        
        # Check for addable metadata.
        self.mutationCounts = self.getMetadataByKey(self.AddableKeys.mutCounts.value, False)
        if self.mutationCounts is not None: self.mutationCounts = int(self.mutationCounts)

        ### Get file paths for useful metadata associated files.

        self.genomeFilePath = os.path.join(getExternalDataDirectory(),self.genomeName,self.genomeName+".fa")

        if self.nucPosName is not None:
            self.baseNucPosFilePath: str = os.path.join(getExternalDataDirectory(), self.genomeName,
                                                   self.nucPosName, self.nucPosName+".bed")
        else: self.baseNucPosFilePath = None

        self.parentDataFilePath = os.path.join(self.directory, self.localParentDataPath)

    class AddableKeys(Enum):

        mutCounts = "mutationCounts"

    # Used to add metadata that cannot be generated when other metadata is initially generated.
    # Note:  The "key" parameter should be a member of the "addableKeys" Enum.
    def addMetadata(self, key: Enum, value):

        # Check the validity of the key, and make sure it doesn't already exist in the metadata.
        if not key in self.AddableKeys: raise ValueError("Given key, \"" + key + "\" is not addable.")
        if self.getMetadataByKey(key.value, False) is not None: raise ValueError("Metadata already exists for key: " + key.value)

        # Append it to the metadata file.
        with open(self.metadataFilePath, 'a') as metadataFile:
            metadataFile.write(key.value + ':\t' + str(value) + '\n')

        # Re-wrap metadata to include this new addition.
        self.wrapMetadataInMembers()