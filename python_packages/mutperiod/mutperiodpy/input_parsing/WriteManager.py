# This script contains a class for managing the writing of input data to various stratified data files.
# Once set up, it only needs to be passed data one line at a time.

import os, subprocess
from typing import IO
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import (Metadata, checkDirs, DataTypeStr, generateFilePath,
                                                                  generateMetadata)
from mutperiodpy.input_parsing.IdentifyMSI import MSIIdentifier
from mutperiodpy.input_parsing.IdentifyMutSigs import MutSigIdentifier


class WriteManager:

    def __init__(self, rootDataDir, context):

        # Get the metadata from the given directory
        self.rootDataDir = rootDataDir
        self.rootMetadata = Metadata(self.rootDataDir)

        # Record the data's context.
        self.context = context

        # create and open the output file in the same directory as the root data.
        self.rootOutputFilePath = generateFilePath(directory = self.rootDataDir, dataGroup = self.rootMetadata.dataGroupName,
                                                   context = self.context, dataType = DataTypeStr.mutations, fileExtension = ".bed")
        self.rootOutputFile = open(self.rootOutputFilePath, 'w')
        self.rootMutCounts = 0

        # By default, all other write options are off unless otherwise specified.
        self.stratifyByIndividualCohorts = False
        self.stratifyByMS = False
        self.stratifyByMutSig = False
        self.stratifyBySignature = False


    # Create the necessary functions to use the class with the "with" keyword.
    def __enter__(self): return self
    
    def __exit__(self, type, value, tb):

        # Close all opened files and sort aggregate files.
        self.cleanupAndSort()


    # Prepares the manager to separate data by cohort.
    def setUpForIndividualCohorts(self):

        self.stratifyByIndividualCohorts = True

        self.currentIndividualCohortID = None # The cohort being written to at a point in time.
        self.currentIndividualCohortFile: IO = None # The open file for the current cohort.
        self.completedIndividualCohorts = dict() # A hashtable of cohorts that have been seen before and should NOT be revisited/rewritten.        

        # Create the directory.
        self.rootIndividualCohortsDirectory = os.path.join(self.rootMetadata.directory,"individual_cohorts")
        checkDirs(self.rootIndividualCohortsDirectory)


    # Prepares the manager to separate cohorts by microsatellite stability.
    # Returns an MSIIdentifier object to be "completed" by the function caller.
    def setUpForMSStratification(self) -> MSIIdentifier:

        self.stratifyByMS = True

        self.MSICohorts = dict() # A hashtable of individual cohorts with MSI

        # Create the necessary directories, file paths, and metadata.
        aggregateMSDirectory = os.path.join(self.rootMetadata.directory,"microsatellite_analysis")
        aggregateMSSDirectory = os.path.join(aggregateMSDirectory,"MSS")
        aggregateMSIDirectory = os.path.join(aggregateMSDirectory,"MSI")
        checkDirs(aggregateMSSDirectory, aggregateMSIDirectory)

        generateMetadata("MSS_" + self.rootMetadata.dataGroupName, self.rootMetadata.genomeName, self.rootMetadata.nucPosName,
                         os.path.join('..','..',self.rootMetadata.localParentDataPath), self.rootMetadata.inputFormat, aggregateMSSDirectory, "MSS")
        generateMetadata("MSI_" + self.rootMetadata.dataGroupName, self.rootMetadata.genomeName, self.rootMetadata.nucPosName,
                         os.path.join('..','..',self.rootMetadata.localParentDataPath), self.rootMetadata.inputFormat, aggregateMSIDirectory, "MSI")
        self.aggregateMSSMutCounts = 0
        self.aggregateMSIMutCounts = 0

        self.aggregateMSSFilePath = generateFilePath(directory = aggregateMSSDirectory, dataGroup = "MSS_" + self.rootMetadata.dataGroupName,
                                                     context = self.context, dataType = DataTypeStr.mutations, fileExtension = ".bed")
        self.aggregateMSSFile = open(self.aggregateMSSFilePath, 'w')
        self.aggregateMSIFilePath = generateFilePath(directory = aggregateMSIDirectory, dataGroup = "MSI_" + self.rootMetadata.dataGroupName,
                                                     context = self.context, dataType = DataTypeStr.mutations, fileExtension = ".bed")
        self.aggregateMSIFile = open(self.aggregateMSIFilePath, 'w')

        # Set up the MSIIdentifier to be returned.
        intermediateFilesDir = os.path.join(self.rootDataDir,"intermediate_files")
        checkDirs(intermediateFilesDir)
        MSISeqInputDataFilePath = generateFilePath(directory = intermediateFilesDir, dataGroup = self.rootMetadata.dataGroupName,
                                                   dataType = "MSISeq_data", fileExtension = ".tsv")
        self.MSICohortsFilePath = generateFilePath(directory = aggregateMSDirectory, dataGroup = self.rootMetadata.dataGroupName, 
                                                   dataType = "MSI_cohorts", fileExtension = ".txt")

        self.myMSIIdentifier = MSIIdentifier(MSISeqInputDataFilePath, self.MSICohortsFilePath)
        return(self.myMSIIdentifier)


    # Prepares the manager to separate cohorts by microsatellite stability.
    # Returns a MutSigIdentifier object to be "completed" by the function caller.
    def setUpForMutSigStratification(self) -> MutSigIdentifier:

        self.stratifyByMutSig = True
        self.mutSigDesignations = dict() # A dictionary of the mutation signatures assigned to each cohort.

        mutSigs = ["1A","1B"] + [str(x) for x in list(range(2,22))] + ["R1","R2","R3","U1","U2"]

        # Create the necessary directories, file paths, and metadata.
        parentMutSigDirectory = os.path.join(self.rootMetadata.directory, "mut_sig_analysis")
        self.mutSigFilePaths = dict()
        self.mutSigFiles = dict()
        self.mutSigMutCounts = dict()

        for mutSig in mutSigs:

            thisMutSigDataGroup = "mut_sig_" + mutSig + '_' + self.rootMetadata.dataGroupName

            # Directory
            thisMutSigDirectory = os.path.join(parentMutSigDirectory,mutSig)
            checkDirs(thisMutSigDirectory)

            # Metadata
            generateMetadata(thisMutSigDataGroup, self.rootMetadata.genomeName, self.rootMetadata.nucPosName,
                             os.path.join('..','..',self.rootMetadata.localParentDataPath), 
                             self.rootMetadata.inputFormat, thisMutSigDirectory, "mutSig" + mutSig)

            # Mutation Counter
            self.mutSigMutCounts[mutSig] = 0

            # File path
            self.mutSigFilePaths[mutSig] = generateFilePath(directory = thisMutSigDirectory, dataGroup = thisMutSigDataGroup, 
                                                            context = self.context, dataType = DataTypeStr.mutations, fileExtension = ".bed")
            self.mutSigFiles[mutSig] = open(self.mutSigFilePaths[mutSig], 'w')

        # Set up the MutSigIdentifier object to be returned.
        intermediateFilesDir = os.path.join(self.rootDataDir,"intermediate_files")
        checkDirs(intermediateFilesDir)

        deconstructSigsInputDataFilePath = generateFilePath(directory = intermediateFilesDir, dataGroup = self.rootMetadata.dataGroupName,
                                                            dataType = "deconstructSigs_data", fileExtension = ".tsv")
        self.mutSigDesignationsFilePath = generateFilePath(directory = parentMutSigDirectory, dataGroup = self.rootMetadata.dataGroupName, 
                                                          dataType = "mut_sig_assignments", fileExtension = ".tsv")

        self.mutSigIdentifier = MutSigIdentifier(deconstructSigsInputDataFilePath, self.mutSigDesignationsFilePath)
        return(self.mutSigIdentifier)


    # Open a new individual cohort file for writing (and close the last one, if applicable.)
    def setUpNewIndividualCohort(self, cohortID):

        # If this isn't the first opened cohort file, close and sort the last one.
        if self.currentIndividualCohortFile is not None: 
            self.currentIndividualCohortFile.close()
            subprocess.run(("sort","-k1,1","-k2,2n",self.individualCohortFilePath,"-o",self.individualCohortFilePath), check = True)
            Metadata(self.individualCohortFilePath).addMetadata(Metadata.AddableKeys.mutCounts, self.currentIndividualCohortMutCounts)

        # Make sure this is actually a new cohort.
        if cohortID in self.completedIndividualCohorts:
            raise ValueError("The cohort " + cohortID + " was encountered in more than one distinct block of data.")
        else:
            self.currentIndividualCohortID = cohortID

        individualCohortDirectory = os.path.join(self.rootIndividualCohortsDirectory,self.currentIndividualCohortID)
        individualCohortDataGroup = self.currentIndividualCohortID + "_" + self.rootMetadata.dataGroupName

        checkDirs(individualCohortDirectory)

        # Determine which other set up "umbrella" cohorts this cohort belongs to.
        cohortMembership = [self.currentIndividualCohortID,]
        if self.stratifyByMS:
            if self.currentIndividualCohortID in self.MSICohorts:
                cohortMembership.append("MSI")
            else:
                cohortMembership.append("MSS")
        if self.stratifyByMutSig:
            if self.currentIndividualCohortID in self.mutSigDesignations:
                for mutSig in self.mutSigDesignations[self.currentIndividualCohortID]:
                    cohortMembership.append("mut_sig_" + mutSig)
                

        # Generate the file path and metadata file and open the file for writing.
        self.individualCohortFilePath = generateFilePath(directory = individualCohortDirectory, dataGroup = individualCohortDataGroup, 
                                                         context = self.context, dataType = DataTypeStr.mutations, fileExtension = ".bed")
        self.currentIndividualCohortFile = open(self.individualCohortFilePath, 'w')
        generateMetadata(individualCohortDataGroup, self.rootMetadata.genomeName, self.rootMetadata.nucPosName,
                         os.path.join("..",self.rootMetadata.localParentDataPath),
                         self.rootMetadata.inputFormat, individualCohortDirectory, *cohortMembership)
        self.currentIndividualCohortMutCounts = 0


    # Writes the given data to all the relevant files based on how the manager was set up.
    def writeData(self, chromosome, startPos, endPos, mutFrom, alteration, strand, cohortID = '.'):

        # Make sure we were passed a line that matches the manager's context.
        assert self.context == float("inf") or self.context == 0 or self.context == len(mutFrom), ("mutFrom: \"" + mutFrom + "\" does not " +
            "fit within the manager's context: " + str(self.context))

        # Then, format it for output.
        outputLine = '\t'.join((chromosome, startPos, endPos, mutFrom, alteration, strand)) + '\n'

        # Write data to the root output file.
        self.rootOutputFile.write(outputLine)
        self.rootMutCounts += 1

        # Write to microsatellite designation if it was set up.
        if self.stratifyByMS and cohortID != '.':

            # Make sure the MSICohorts hashtable has actually been propogated.
            if len(self.MSICohorts) == 0:
                if not self.myMSIIdentifier.MSICohortsIdentified:
                    raise ValueError("MSIIdentifier protocol was never completed.")
                else:
                    with open(self.MSICohortsFilePath, 'r') as MSICohortsFile:
                        for line in MSICohortsFile:
                            self.MSICohorts[line.strip()] = None

            if cohortID in self.MSICohorts:
                self.aggregateMSIFile.write(outputLine)
                self.aggregateMSIMutCounts += 1
            else:
                self.aggregateMSSFile.write(outputLine)
                self.aggregateMSSMutCounts += 1

        # Write to signature designations if it was set up.
        if self.stratifyByMutSig and cohortID != '.':

            # Propogate the mutSigDesignations dictionary if it hasn't been.
            if len(self.mutSigDesignations) == 0:
                if not self.mutSigIdentifier.mutSigsIdentified:
                    raise ValueError("MutSigIdentifier protocol was never completed.")
                else:

                    with open(self.mutSigDesignationsFilePath, 'r') as mutSigDesignationsFile:
                        for line in mutSigDesignationsFile:

                            mutSigsCohortID = str(line).strip().split('\t')[0]
                            mutSigs = str(line).strip().split('\t')[1].split(',')

                            if mutSigs[0] == "None": continue
                            else: self.mutSigDesignations[mutSigsCohortID] = mutSigs

            # Write the data to its appropriate mutsig file(s) if its cohort has mutsig designations.
            if cohortID in self.mutSigDesignations:
                
                for mutSig in self.mutSigDesignations[cohortID]:
                    self.mutSigFiles[mutSig].write(outputLine)
                    self.mutSigMutCounts[mutSig] += 1


        # Write to individual cohorts if desired.
        if self.stratifyByIndividualCohorts and cohortID != '.':
            
            # Check to see if we've reached a new cohortID .
            if self.currentIndividualCohortID is None or self.currentIndividualCohortID != cohortID:
                self.completedIndividualCohorts[self.currentIndividualCohortID] = None
                self.setUpNewIndividualCohort(cohortID)

            # Write to the current individual cohort.
            self.currentIndividualCohortFile.write(outputLine)
            self.currentIndividualCohortMutCounts += 1


    # Closes open files to clean up the class after it's done being used.
    def cleanupAndSort(self):

        self.rootOutputFile.close()
        subprocess.run(("sort","-k1,1","-k2,2n",self.rootOutputFilePath,"-o",self.rootOutputFilePath), check = True)
        Metadata(self.rootOutputFilePath).addMetadata(Metadata.AddableKeys.mutCounts, self.rootMutCounts)
        

        if self.stratifyByMS:
            self.aggregateMSIFile.close()
            subprocess.run(("sort","-k1,1","-k2,2n",self.aggregateMSIFilePath,"-o",self.aggregateMSIFilePath), check = True)
            Metadata(self.aggregateMSIFilePath).addMetadata(Metadata.AddableKeys.mutCounts, self.aggregateMSIMutCounts)
            self.aggregateMSSFile.close()
            subprocess.run(("sort","-k1,1","-k2,2n",self.aggregateMSSFilePath,"-o",self.aggregateMSSFilePath), check = True)
            Metadata(self.aggregateMSSFilePath).addMetadata(Metadata.AddableKeys.mutCounts, self.aggregateMSSMutCounts)

        if self.stratifyByMutSig:
            for mutSig in self.mutSigFiles:
                self.mutSigFiles[mutSig].close()
                subprocess.run(("sort","-k1,1","-k2,2n",self.mutSigFilePaths[mutSig], "-o",self.mutSigFilePaths[mutSig]), check = True)
                Metadata(self.mutSigFilePaths[mutSig]).addMetadata(Metadata.AddableKeys.mutCounts, self.mutSigMutCounts[mutSig])

        if self.stratifyByIndividualCohorts and self.currentIndividualCohortFile is not None:
            self.currentIndividualCohortFile.close()
            subprocess.run(("sort","-k1,1","-k2,2n",self.individualCohortFilePath,"-o",self.individualCohortFilePath), check = True)
            Metadata(self.individualCohortFilePath).addMetadata(Metadata.AddableKeys.mutCounts, self.currentIndividualCohortMutCounts)