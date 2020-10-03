# This script takes stratified mutation data and uses it to generate an input file for the deconstructSigs R package which
# then assigns the most prominent mutation signature(s) for each cohort.

import subprocess, os
from nucperiodpy.helper_scripts.UsefulFileSystemFunctions import RPackagesDirectory


class MutSigIdentifier:

    def __init__(self, deconstructSigsInputDataFilePath, deconstructSigsOutputFilePath):
        
        self.deconstructSigsInputDataFilePath = deconstructSigsInputDataFilePath
        self.deconstructSigsInputDataFile = open(self.deconstructSigsInputDataFilePath, 'w')
        self.mutSigsIdentified = False

        self.deconstructSigsOutputFilePath = deconstructSigsOutputFilePath


    # Functions for "with" compatibility.
    def __enter__(self): return self
    
    def __exit__(self, type, value, tb):

        # Close all opened files.
        self.cleanup()


    # Adds an entry to the deconstructSigs input data file.
    def addData(self, cohortID, chromosome, base1Pos, referenceBase, mutantBase):

        if (referenceBase not in ('A','C','G','T') or mutantBase not in ('A','C','G','T')):
            raise ValueError("Given entry should represent an SNP.  \"" + referenceBase + "\" > \"" + 
                             mutantBase + "\" does not represent an SNP.")

        if not self.mutSigsIdentified:
            self.deconstructSigsInputDataFile.write('\t'.join((cohortID, chromosome, base1Pos, referenceBase, mutantBase)) + '\n')
        else:
            raise ValueError("deconstructSigs has already run and the input data file has already been closed!")


    # Generates the mutation signatures output file from the given data.
    def identifyMutSigs(self, verbose = True):

        if self.mutSigsIdentified:
            raise ValueError("Mutation signatures have already been identified.")
        else:
            self.deconstructSigsInputDataFile.close()

        if verbose: print("Calling deconstructSigs to identify mutation signatures")
        subprocess.run(" ".join(("Rscript",os.path.join(RPackagesDirectory,"RunNucPeriod","GetMutSigs.R"),
                                self.deconstructSigsInputDataFilePath, self.deconstructSigsOutputFilePath)), shell = True, check = True)

        self.mutSigsIdentified = True
    
    # Close up any open files.
    def cleanup(self):

        if not self.mutSigsIdentified:       
            self.deconstructSigsInputDataFile.close()