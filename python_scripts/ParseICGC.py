# This script reads one or more "simple somatic mutation" data file(s) from ICGC and 
# writes information on single base substitution mutations to a new bed file or files for further analysis.
import os, gzip, subprocess
from TkinterDialog import TkinterDialog, Selections
from UsefulFileSystemFunctions import dataTypes, generateFilePath, generateMetadata, getIsolatedParentDir
from UsefulBioinformaticsFunctions import reverseCompliment, isPurine, baseChromosomes
from typing import IO, List

# This class represents the mutation data obtained from ICGC in a more precise form.
# It also contains functions to represent that data in the exact format required for bed files or MSIseq, if they are valid.
class ICGCMutation:

    DNABases = ('A','C','T','G') # A list of single DNA nucleotides to check against.

    def __init__(self, choppedUpLine):

        # The relevant data extracted from a given line in the ICGC file.
        self.donorID = choppedUpLine[1]
        self.chromosome = "chr" + choppedUpLine[8]
        if self.chromosome == "chrMT": self.chromosome = "chrM" # Standardize the naming of mitochondrial chromosomes.
        self.startPos = choppedUpLine[9]
        self.endPos = choppedUpLine[10]
        self.mutationType = choppedUpLine[13]
        self.mutatedFrom = choppedUpLine[15]
        self.mutatedTo = choppedUpLine[16]
        
        if choppedUpLine[11] != "1":
            raise ValueError("Error.  Strand field does not contain \"1\" for plus strand.  " + 
                             "Found " + str(choppedUpLine[11]) + " instead.")


    # Is the mutation valid for a bed file of single nucleotide substitutions?
    def isValidForBed(self):
        # Is this a single nucleotide substitution?
        return (self.mutatedFrom.upper() in ICGCMutation.DNABases and self.mutatedTo.upper() in ICGCMutation.DNABases)


    # Return a string that contains the mutation data in tab delimited bed form.
    def formatForBed(self):

        # Did you remember to check to see if it was a single base substitution?
        if not self.isValidForBed(): raise ValueError("Error:  The given line is not valid for bed format.")

        # Piece together the bed formatted line and then return it.

        mutationPos0Base = str(int(self.startPos)-1)

        # If the mutated base is listed as arising from a purine, flip the mutation and the strand.
        if isPurine(self.mutatedFrom):
            mutation = reverseCompliment(self.mutatedFrom) + '>' + reverseCompliment(self.mutatedTo)
            strand = '-'
        else:
            mutation = self.mutatedFrom + '>' + self.mutatedTo
            strand = '+'

        return( '\t'.join((self.chromosome,mutationPos0Base,self.startPos,mutation[0],mutation,strand)) )

    # Return a string that contains the mutation data in tab delimited MSIseq form.
    def formatForMSIseq(self):
        
        # Determine what type of mutation we have.
        if "substitution" in self.mutationType: MSIseqMutationType = "SNP"
        elif "insertion" in self.mutationType: MSIseqMutationType = "INS"
        elif "deletion" in self.mutationType: MSIseqMutationType = "DEL"
        else: raise ValueError("Error: Mutation type, \"" + self.mutationType + "\", does not appear to be a substitution, deletion, or insertion.")

        # Piece together the MSIseq formatted line and then return it.
        return( '\t'.join((self.chromosome,self.startPos,self.endPos,MSIseqMutationType,self.donorID)) )

    # Return a string that contains the data necessary for deisgnating mutational signatures.
    def formatForMutSig(self):
        # Piece together the MutSig formatted line and then return it.
        return( '\t'.join((self.donorID,self.chromosome,self.startPos,self.mutatedFrom,self.mutatedTo)) )


# This class takes an ICGCFile object and, when iterated over, returns exactly once each mutation present in a donor that is
# the result of whole genome sequencing using GRCh37 as the reference genome.
# The mutation is returned as an ICGCMutation object.
class ICGCIterator:

    def __init__(self, ICGCFile: IO):
        
        self.ICGCFile = ICGCFile # The file containing the ICGC data
        self.finishedDonors = list() # A list of donors to make sure we don't encounter one more than once.
        self.currentDonor = '' # The donorID currently being read for mutation data.
        self.currentDonorMutations = dict() # A dictionary of mutations unique to the current donor, to avoid writing duplicate mutations.
        self.previousDonorMutations: List[ICGCMutation] = list() # A list that keeps the last donor's mutations in order to write data on individual donors all at once.

    # Make the class iteratable, returning unique data lines one at a time.
    def __iter__(self):
        return self
    def __next__(self):

        # Read in lines until an acceptable one is found.
        while(True):

            # Read in the next line in the file.
            line = self.ICGCFile.readline()
            # Split the line into its individual data components.  The relevant components (with indices) are:
            #   0: mutation ID
            #   1: donor ID
            #   8: chromosome
            #   9: mutation start position (inclusive)
            #   10: mutation end position (inclusive)
            #   13: mutation type
            #   15: mutated from
            #   16: mutated to
            #   11: strand
            #   12: reference genome version
            #   33: sequencing method
            choppedUpLine = str(line,"utf-8").strip().split('\t')
            if len(choppedUpLine) < 34: raise StopIteration # Check to see if we have reached the end of the file.


            if (choppedUpLine[12] == "GRCh37" and # Is the reference genome hg19?
                choppedUpLine[33] == "WGS" and # Was whole genome sequencing used to generate the data?
                ("chr" + choppedUpLine[8]) in baseChromosomes): # Is the chromosome acceptable?

                # Is this a new donor?
                if choppedUpLine[1] != self.currentDonor:
                    self.finishedDonors.append(self.currentDonor)
                    self.currentDonor = choppedUpLine[1]
                    self.previousDonorMutations = list(self.currentDonorMutations.values())
                    self.currentDonorMutations.clear()
                    if self.currentDonor in self.finishedDonors:
                        raise ValueError("Error:  Donor " + self.currentDonor + " is present in more than one block of data!")
                    print("Reading and writing from donor",self.currentDonor)

                # Have we seen this mutation in this donor?
                if choppedUpLine[0] not in self.currentDonorMutations:
                    newMutation = ICGCMutation(choppedUpLine)
                    self.currentDonorMutations[choppedUpLine[0]] = newMutation
                    return newMutation


# This class handles the (somewhat complex) file systems for the parsing process and makes output file paths easily acessible.
# It takes the directory of the ICGC file as input as well as boolean values to determine how the data is being parsed
# which determines what directories and file paths need to be generated.
# It also takes paths to the associated genome and nucleosome positions files which are used when generating metadata.
class ICGCParserFileManager:

    def __init__(self, ICGCFilePath, convertToBed: bool, convertToMutSig: bool, separatebyMSI: bool, createIndividualDonorFiles: bool,
                 associatedGenomeFilePath: str, associatedNucleosomePositionsFilePath: str):

        self.localRootDirectory = os.path.dirname(ICGCFilePath) # The root directory for this localized file system
        self.mutationGroupName = getIsolatedParentDir(ICGCFilePath) # The name of the mutation group the directory pertains to

        # information used later for generating metadata of sub-groups
        self.associatedGenome = getIsolatedParentDir(associatedGenomeFilePath)
        self.associatedNucleosomePositions = getIsolatedParentDir(associatedNucleosomePositionsFilePath)
        self.ICGCFileBasename = os.path.basename(ICGCFilePath) 

        # Generate the metadata for the main data group.
        generateMetadata(self.mutationGroupName, self.associatedGenome, self.associatedNucleosomePositions, 
                         self.ICGCFileBasename, os.path.dirname(ICGCFilePath))

        if convertToBed: 
            self.prepForBedFileOutput(createIndividualDonorFiles)
            if separatebyMSI: self.prepForMSISeparation(createIndividualDonorFiles)

        if convertToMutSig: self.prepForMutSigFileOutput()


    # Prepare the file system for bed output by generating the paths to the bed ouput files.
    # (Optionally, also prepare the system to output individual donor files)
    def prepForBedFileOutput(self, createIndividualDonorFiles):

        self.bedFilePath = generateFilePath(directory = self.localRootDirectory, dataGroup = self.mutationGroupName,
                                            context = "singlenuc", dataType = dataTypes.mutations, fileExtension = ".bed")

        if createIndividualDonorFiles:

            self.individualDonorsDirectory = os.path.join(self.localRootDirectory,"individual_donors","non-stratified")
            if not os.path.exists(self.individualDonorsDirectory): os.makedirs(self.individualDonorsDirectory)


    # Prepare the file system to separate bed files based on microsatellite stability and generate the paths to output files.
    # (Optionally, also prepare the system to output individual donor files)
    def prepForMSISeparation(self, createIndividualDonorFiles):
        
        self.MSIseqDataFilePath = os.path.join(self.localRootDirectory,self.mutationGroupName+"_MSIseq_data.tsv")
        self.MSIDonorListFilePath = os.path.join(self.localRootDirectory,self.mutationGroupName+"_MSI_donors.txt")

        microsatteliteDirectory = os.path.join(self.localRootDirectory,"microsatellite_analysis")
        if not os.path.exists(microsatteliteDirectory): os.mkdir(microsatteliteDirectory)

        MSIDirectory = os.path.join(microsatteliteDirectory,"MSI")
        if not os.path.exists(MSIDirectory): os.mkdir(MSIDirectory)
        MSSDirectory = os.path.join(microsatteliteDirectory,"MSS")
        if not os.path.exists(MSSDirectory): os.mkdir(MSSDirectory)

        self.MSIMutationGroupName = "MSI_" + self.mutationGroupName
        self.MSSMutationGroupName = "MSS_" + self.mutationGroupName

        # Generate the metadata for each of the new data groups.
        generateMetadata(self.MSIMutationGroupName, self.associatedGenome, self.associatedNucleosomePositions,
                         os.path.join('..','..',self.ICGCFileBasename), MSIDirectory)
        generateMetadata(self.MSSMutationGroupName, self.associatedGenome, self.associatedNucleosomePositions,
                         os.path.join('..','..',self.ICGCFileBasename), MSSDirectory)

        self.MSIBedFilePath = generateFilePath(directory = MSIDirectory, dataGroup = self.MSIMutationGroupName,
                                               context = "singlenuc", dataType = dataTypes.mutations, fileExtension = ".bed")
        self.MSSBedFilePath = generateFilePath(directory = MSSDirectory, dataGroup = self.MSSMutationGroupName,
                                               context = "singlenuc", dataType = dataTypes.mutations, fileExtension = ".bed")

        if createIndividualDonorFiles:

            self.individualMSIDonorsDirectory = os.path.join(os.path.dirname(self.individualDonorsDirectory),"MSI")
            self.individualMSSDonorsDirectory = os.path.join(os.path.dirname(self.individualDonorsDirectory),"MSS")

            if not os.path.exists(self.individualMSIDonorsDirectory): os.makedirs(self.individualMSIDonorsDirectory)
            if not os.path.exists(self.individualMSSDonorsDirectory): os.makedirs(self.individualMSSDonorsDirectory)


    # Prepare the file system for MutSig (deconstructSigs R package) output..
    def prepForMutSigFileOutput(self):
        self.mutSigDataFilePath = os.path.join(self.localRootDirectory,self.mutationGroupName+"_MutSig_data.tsv")


    # Generates a file path for the given donor and returns it.
    # If MSI is NoneType, the file is made without specification of microsatellite stability.
    # Otherwise, the MSI argument dictates whether the donr is classified as microsatellite instable (true) or stable (false).
    def getIndividualDonorFilePath(self, donorID: str, MSI = None):

        if MSI is None: 
            donorDirectory = os.path.join(self.individualDonorsDirectory,donorID)
            donorDataGroup = donorID + "_" + self.mutationGroupName
        elif MSI: 
            donorDirectory = os.path.join(self.individualMSIDonorsDirectory,donorID)
            donorDataGroup = donorID + "_" + self.MSIMutationGroupName
        else: 
            donorDirectory = os.path.join(self.individualMSSDonorsDirectory,donorID)
            donorDataGroup = donorID + "_" + self.MSSMutationGroupName

        if not os.path.exists(donorDirectory): os.mkdir(donorDirectory)

        # Generate the file path and metadata file
        donorFilePath = generateFilePath(directory = donorDirectory, dataGroup = donorDataGroup, context = "singlenuc",
                                         dataType = dataTypes.mutations, fileExtension = ".bed")
        generateMetadata(donorDataGroup, self.associatedGenome, self.associatedNucleosomePositions,
                         os.path.join("..","..","..",self.ICGCFileBasename), donorDirectory)

        return donorFilePath


# Creates an MSI donor list by parsing the ICGC file for MSIseq and running the result through the corresponding R script.
def generateMSIDonorList(ICGCFilePath, fileManager: ICGCParserFileManager):
    
    # Generate the MSIseq input file.
    print("Creating MSIseq input file...")
    with gzip.open(ICGCFilePath, 'r') as ICGCFile:
        with open(fileManager.MSIseqDataFilePath, 'w') as MSIseqDataFile:
            iterator = ICGCIterator(ICGCFile)
            for mutation in iterator:
                MSIseqDataFile.write(mutation.formatForMSIseq() + '\n')

    # Sort the data
    print("Sorting MSIseq data...")
    subprocess.run(" ".join(("sort","-k5,5","-k2,2n",fileManager.MSIseqDataFilePath,"-o",fileManager.MSIseqDataFilePath)), 
                   shell = True, check = True)

    # Pass the path to the newly created MSIseq data file to the R script to generate the MSI Donor list
    print("Calling MSIseq to generate MSI donor list...")
    subprocess.run(" ".join(("Rscript",os.path.join(os.path.dirname(__file__),"..","R_scripts","RunNucPeriod","FindMSIDonors.R"),
                   fileManager.MSIseqDataFilePath)), shell = True, check = True)


# Handles the basic parsing of the script.
def parseICGC(ICGCFilePaths, genomeFilePath, nucPosFilePath, convertToBed, separateByMSI, 
              createIndividualDonorFiles, convertToMutSig):

    if not (convertToBed or convertToMutSig): raise ValueError("Error: No output format selected.")
    if separateByMSI and not convertToBed: raise ValueError("Error: Cannot separate by MSI without converting to bed.")
    if createIndividualDonorFiles and not convertToBed: raise ValueError("Error: Cannot create individual bed files without converting to bed.")

    # Run the parser for each ICGC file given.
    for ICGCFilePath in ICGCFilePaths:

        print("\nWorking in:",os.path.split(ICGCFilePath)[1])

        if not str(ICGCFilePath).endswith(".gz"):
            raise ValueError("Error:  Expected ICGC file to be gzipped (.gz file format).")
        if not "simple_somatic_mutation" in os.path.split(ICGCFilePath)[1]:
            raise ValueError("Error:  Expected ICGC file with \"simple_somatic_mutation\" in the name.\n" +
                            "Note: if a directory was specified to search for ICGC input files, all files ending in .tsv.gz\
                            are considered valid.")

        # Prepare the file system...
        fileManager = ICGCParserFileManager(ICGCFilePath, convertToBed, convertToMutSig, separateByMSI, 
                                            createIndividualDonorFiles, genomeFilePath, nucPosFilePath)

        # Create the MSI donor list if necessary.
        if separateByMSI and not os.path.exists(fileManager.MSIDonorListFilePath):
            print("MSI donor list not found at ", fileManager.MSIDonorListFilePath,".",sep='')
            print("Generating MSI donor list...")
            generateMSIDonorList(ICGCFilePath,fileManager)

        # Read the MSI donor list into a dictionary object.
        MSIDonors = dict()
        if separateByMSI:
            with open(fileManager.MSIDonorListFilePath, 'r') as MSIDonorListFile:
                for donor in MSIDonorListFile: MSIDonors[donor.strip()] = None


        # This function handles the writing of mutation data to individual donor files.
        def writeIndividualDonorFile(donorID, mutations):

            # Get the file paths
            if separateByMSI:
                if donorID in MSIDonors: individualDonorBedFilePath = fileManager.getIndividualDonorFilePath(donorID,MSI = True)
                else: individualDonorBedFilePath = fileManager.getIndividualDonorFilePath(donorID,MSI = False)
            else: individualDonorBedFilePath = fileManager.getIndividualDonorFilePath(donorID)

            # Open
            individualDonorBedFile = open(individualDonorBedFilePath,'w')           

            # Write
            for mutation in mutations:

                if mutation.isValidForBed():
                    individualDonorBedFile.write(mutation.formatForBed() + '\n')

            # Close
            individualDonorBedFile.close()

            # Sort
            subprocess.run(" ".join(("sort","-k1,1","-k2,2n",individualDonorBedFilePath,"-o",individualDonorBedFilePath)),
                           shell = True, check = True)

        # Parse the file!
        with gzip.open(ICGCFilePath, 'r') as ICGCFile:

            # Open the files to write to.
            if convertToBed: 
                bedFile = open(fileManager.bedFilePath, 'w')

                if separateByMSI: 
                    MSIBedFile = open(fileManager.MSIBedFilePath, 'w')
                    MSSBedFile = open(fileManager.MSSBedFilePath, 'w')

            if convertToMutSig: mutSigDataFile = open(fileManager.mutSigDataFilePath, 'w')
            
            iterator = ICGCIterator(ICGCFile) # The object to read through the ICGC data and discard duplicate mutations for each donor.
            currentDonorID = None # Keeps track of when all of a donor's data has been read.       

            # Loop through the ICGC data, writing data to the opened files.
            for mutation in iterator:

                if convertToBed and mutation.isValidForBed(): 
                    bedFile.write(mutation.formatForBed() + '\n')

                    if separateByMSI:
                        if mutation.donorID in MSIDonors: 
                            MSIBedFile.write(mutation.formatForBed() + '\n')
                        else: MSSBedFile.write(mutation.formatForBed() + '\n')

                if convertToMutSig: mutSigDataFile.write(mutation.formatForMutSig() + '\n')

                # Take these extra steps to stratify data by donorID if it was requested by the user and we have read all of a donor's data.
                if currentDonorID is None: currentDonorID = mutation.donorID
                elif createIndividualDonorFiles and mutation.donorID != currentDonorID:
                    writeIndividualDonorFile(currentDonorID, iterator.previousDonorMutations)
                    currentDonorID = mutation.donorID

            # Make sure we don't miss the last donor for individual files!
            if createIndividualDonorFiles: writeIndividualDonorFile(currentDonorID, iterator.currentDonorMutations.values())

            # Close the files that were opened.
            if convertToBed: 
                bedFile.close()

                if separateByMSI: 
                    MSIBedFile.close()
                    MSSBedFile.close()

            if convertToMutSig: mutSigDataFile.close()

        # Sort the mutation data using linux sort, because it's not quite so greedy with memory usage...
        if convertToBed:
            print("Sorting bed data...")
            subprocess.run(" ".join(("sort","-k1,1","-k2,2n",fileManager.bedFilePath,"-o",fileManager.bedFilePath)), shell = True, check = True)

            if separateByMSI:
                subprocess.run(" ".join(("sort","-k1,1","-k2,2n",fileManager.MSIBedFilePath,"-o",fileManager.MSIBedFilePath)), shell = True, check = True)
                subprocess.run(" ".join(("sort","-k1,1","-k2,2n",fileManager.MSSBedFilePath,"-o",fileManager.MSSBedFilePath)), shell = True, check = True)

        if convertToMutSig:
            print("Sorting MutSig data...")
            subprocess.run(" ".join(("sort","-k2,2","-k3,3n",fileManager.mutSigDataFilePath,"-o",fileManager.mutSigDataFilePath)), shell = True, check = True)


if __name__ == "__main__":

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(__file__),"..","data"))
    dialog.createMultipleFileSelector("ICGC Mutation Files:",0,".tsv.gz",("gzip files",".gz"))
    dialog.createFileSelector("Genome Fasta File:",1,("Fasta Files",".fa"))
    dialog.createFileSelector("Strongly Positioned Nucleosome File:",2,("Bed Files",".bed"))
    dialog.createCheckbox("Convert to Bed SNPs", 3, 0)
    dialog.createCheckbox("Separate Bed Mutations by Microsatellite Stability", 3, 1)
    dialog.createCheckbox("Also Create individual bed files for each donor.",4,0)
    dialog.createCheckbox("Convert to MutSig format", 4, 1)
    dialog.createReturnButton(5,0)
    dialog.createQuitButton(5,2)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    ICGCFilePaths = list(selections.getFilePathGroups())[0] # A list of ICGC mutation file paths
    genomeFilePath = list(selections.getIndividualFilePaths())[0]
    nucPosFilePath = list(selections.getIndividualFilePaths())[1]
    convertToBed = list(selections.getToggleStates())[0]
    separateByMSI = list(selections.getToggleStates())[1]
    createIndividualDonorFiles = list(selections.getToggleStates())[2]
    convertToMutSig = list(selections.getToggleStates())[3]

    parseICGC(ICGCFilePaths, genomeFilePath, nucPosFilePath, convertToBed, 
              separateByMSI, createIndividualDonorFiles, convertToMutSig)