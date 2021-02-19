# This script reads the structure of the mutperiodpy file system to determine the presence/status of 
# data sets and analyses in the project.
# It also provides functionality to delete data sets and analyses.

import os, shutil
from enum import Enum
from mutperiodpy.helper_scripts import UsefulFileSystemFunctions as UFSF
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import dataDirectory, DataTypeStr, getIsolatedParentDir, InputFormat
from mutperiodpy.input_parsing.ParseICGC import parseICGC
from mutperiodpy.input_parsing.ParseCustomBed import parseCustomBed
from mutperiodpy.input_parsing.ParseTXRSeq import parseTXRSeq


# Data stratification identifiers and their relevant strings.
class Stratification(Enum):

    mutSig = "mutation signatures"
    MSAnalysis = "microsatellite analysis"
    individualCohorts = "individual cohorts"


# Data normalization identifiers and their relevant strings.
class Normalization(Enum):

    singlenuc = "single nucleotide context"
    trinuc = "trinucleotide context"
    pentanuc = "pentanucleotide context"


# Generates counting radius identifiers and their relevant strings.
class Radius:

    def __init__(self, nucGroup, linkerOffset = 0):
        
        if nucGroup and linkerOffset != 0:
            raise ValueError("nuc-group should not have any linker offset.")

        self.nucGroup = nucGroup
        self.linkerOffset = linkerOffset

    def value(self):
        
        if nucGroup: return("nucleosome group (1000 bp)")

        elif linkerOffset == 0: return("single nucleosome (73 bp)")
        else: 
            return("single nucleosome +" + str(linkerOffset) + " bp linker (" + str(73 + linkerOffset) + " bp)")

    def __eq__(self, otherRadius): 
        if isinstance(otherRadius, Radius):
            return (self.nucGroup == otherRadius.nucGroup and self.linkerOffset == otherRadius.linkerOffset)
        else: return False


# Keeps track of the ways the data set has been configured (stratifications, flanking nucleotide context, etc.)
class DataSetConfiguration:

    def __init__(self, directory = None):
        
        # Initialize the variables that store information on the configuration.
        self.directory = directory
        self.stratifications = list()
        self.normalizations = list()
        self.radiuses = list()

        # Auto-acquire information on the configuration from the file system if a directory was supplied.
        if directory is not None:
            self.findStratifications()
            self.findNormalizations()
            self.findRadiuses()

    # Adds strings to the stratifications list for every stratification imposed on the data set.
    def findStratifications(self):

        if os.path.isdir(os.path.join(self.directory,"mut_sig_analysis")):
            self.stratifications.append(Stratification.mutSig)

        if os.path.isdir(os.path.join(self.directory,"individual_cohorts")):
            self.stratifications.append(Stratification.individualCohorts)

        if os.path.isdir(os.path.join(self.directory,"microsatellite_analysis")):
            self.stratifications.append(Stratification.MSAnalysis)


    # Adds strings to the stratifications list for every stratification imposed on the data set.
    def findNormalizations(self):

        # Search within normalized counts files for markers of the associated context.
        for item in os.listdir(self.directory):
                
            if not DataTypeStr.mutations in item:

                if "singlenuc" in item and Normalization.singlenuc not in self.normalizations:
                    self.normalizations.append(Normalization.singlenuc)
                elif "trinuc" in item and Normalization.trinuc not in self.normalizations:
                    self.normalizations.append(Normalization.trinuc)
                elif "pentanuc" in item and Normalization.pentanuc not in self.normalizations:
                    self.normalizations.append(Normalization.pentanuc)


    # Adds strings to the stratifications list for every stratification imposed on the data set.
    def findRadiuses(self):

        # Search within normalized counts files for markers of the associated counting radius.
        for item in os.listdir(self.directory):

            if DataTypeStr.normNucCounts in item:

                if "nuc-group" in item and Radius(True) not in self.radiuses:
                    self.radiuses.append(Radius(True))
                else:

                    linkerOffset = UFSF.getLinkerOffset(item) 
                    newRadius = Radius(False, linkerOffset)

                    if newRadius not in self.radiuses: self.radiuses.append(newRadius)

        if len(self.radiuses) == 0: raise ValueError("No counting radius found.")


# Keeps track of the properties for a given data set. (Mostly the configuration, which is itself a separate class)
class DataSet:

    def __init__(self, topLevelDirectory, configuration: DataSetConfiguration = None):

        self.topLevelDirectory = topLevelDirectory

        # Auto-acquire the configuration if none was supplied.
        if configuration is None: self.configuration = DataSetConfiguration(self.topLevelDirectory)
        else: self.configuration: DataSetConfiguration = configuration
        
        self.name = getIsolatedParentDir(self.topLevelDirectory, True)
        self.topLevelMetadata = UFSF.Metadata(self.topLevelDirectory)



    def getStratificationsString(self):
        print()


    # Get strings of relevant metadata output for the user.
    def getMetadataStrings(self):   

        metadataStrings = list()
        metadataStrings.append("Associated Genome: " + self.topLevelMetadata.genomeName)
        metadataStrings.append("Associated Nucleosomes: " + self.topLevelMetadata.nucPosName)
        metadataStrings.append("Date Created: " + self.topLevelMetadata.dateTime)
        return metadataStrings


class FileSystemManager:

    def __init__(self):
        
        # Find data sets in the data directory.
        self.dataSets = list()
        findDataSets(dataDirectory)

        self.illegalNamingSubstrings = ["linker+", "singlenuc", "trinuc", "pentanuc", "nuc-group", DataTypeStr.customInput, DataTypeStr.mutations,
                                        DataTypeStr.mutBackground, DataTypeStr.normNucCounts, DataTypeStr.nucMutBackground, DataTypeStr.rawNucCounts]
        self.illegalNames = [dataSet.name for dataSet in self.dataSets] 
        self.illegalNames += [".metadata", "intermediate_files", "microsatellite_analysis", "mut_sig_analysis"]


    # Creates DataSet objects for every data set in the data directory of the project.
    def findDataSets(self, directory):

        # First, look for a metadata file.  If one is here, create a DataSet object with the current directory.
        if os.path.isfile(os.path.join(directory,".metadata")):
            self.dataSets.append(DataSet(directory))

        # If this directory didn't have a metadata file, search through sub directories of the data directory.  
        else:

            for item in os.listdir(directory):

                path = os.path.join(directory,item)

                # Recursively search any directories
                if os.path.isdir(path):
                    self.findDataSets(path)

    
    # Creates a new data set 
    def createNewDataSet(self, name, inputDataPath, inputFormat, associatedGenome, associatedNucPos, configuration):

        # Make sure the name is legal.
        for illegalNamingSubstring in self.illegalNamingSubstrings:
            if illegalNamingSubstring in name: 
                raise ValueError(illegalNamingSubstring + " is a key substring and cannot be part of a data set name.")
        
        for illegalName in self.illegalNames:
            if name == illegalName:
                raise ValueError(name + " is an illegal name, as it matches another vital file or directory name.")

        # Create the new data set directory and copy the input data into it.
        dataSetDirectory = os.path.join(dataDirectory,name)
        os.mkdir(dataSetDirectory)
        newInputDataPath = os.path.join(dataSetDirectory, os.path.basename(inputDataPath))
        shutil.copy(inputDataPath, newInputDataPath)

            

    
    # Parses the given input data to kick off the pipeline.
    def parseInputData(self, inputDataPath, inputFormat, associatedGenome, associatedNucPos, configuration: DataSetConfiguration):

        # Determine what, if any, stratifications were requested in the given configuration.
        

        # Based on the inputFormat, call the corresponding input_parsing function to initiate the data in the pipeline.
        if inputFormat == InputFormat.ICGC:
            parseICGC(newInputDataPath, associatedGenome, associatedNucPos, )
        elif inputFormat == InputFormat.customBed:
            parseCustomBed(newInputDataPath, associatedGenome, associatedNucPos, )
        elif inputFormat == InputFormat.tXRSeq:
            parseTXRSeq()


    # Runs the rest of the pipeline according to the given configuration parameters.
    def runPipeline(self, dataSet, configuration):
        print()


    # Reconfigures a data set by either adding/regenerating the given configuration or 
    # removing the components of the configuration (if specified in the "remove" parameter)
    def configureDataSet(self, dataSet, configuration, remove = False):
        print()
