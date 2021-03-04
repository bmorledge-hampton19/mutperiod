# This script takes stratified mutation data and uses it to generate an input file for the MSIseq R package which
# then generates a list of MSI cohorts.

import subprocess, os
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import rScriptsDirectory


class MSIIdentifier:

    def __init__(self, MSISeqInputDataFilePath, MSICohortsFilePath):
        
        self.MSISeqInputDataFilePath = MSISeqInputDataFilePath
        self.MSISeqInputDataFile = open(self.MSISeqInputDataFilePath, 'w')
        self.MSICohortsIdentified = False

        self.MSICohortsFilePath = MSICohortsFilePath


    # Functions for "with" compatibility.
    def __enter__(self): return self
    
    def __exit__(self, type, value, tb):

        # Close all opened files.
        self.cleanup()


    # Adds an entry to the MSISeq data file.
    # NOTE: startPos and endPos should both be 1-based.
    def addData(self, chromosome, startPos, endPos, mutType, cohortID):

        if mutType not in ("SNP","INS","DEL"):
            raise ValueError(mutType + " is not a valid mutType value. mutType should be SNP, INS, or DEL.")

        if not self.MSICohortsIdentified:
            self.MSISeqInputDataFile.write('\t'.join((chromosome, startPos, endPos, mutType, cohortID)) + '\n')
        else:
            raise ValueError("MSI cohorts already been identified and MSISeq data file already closed!")


    # Generates the MSI cohorts list from the given data.
    def identifyMSICohorts(self, verbose = True):

        if self.MSICohortsIdentified:
            raise ValueError("MSI cohorts have already been identified.")
        else:
            self.MSISeqInputDataFile.close()

        if verbose: print("Calling MSIseq to generate MSI donor list...")
        subprocess.run(("Rscript",os.path.join(rScriptsDirectory,"FindMSIDonors.R"),
                        self.MSISeqInputDataFilePath, self.MSICohortsFilePath), check = True)

        self.MSICohortsIdentified = True
    
    # Close up any open files.
    def cleanup(self):

        if not self.MSICohortsIdentified:       
            self.MSISeqInputDataFile.close()