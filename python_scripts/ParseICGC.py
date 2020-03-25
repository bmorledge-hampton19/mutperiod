# This script reads one or more "simple somatic mutation" data file(s) from ICGC and 
# writes information on single base substitution mutations to a new bed file or files for further analysis.
import os, gzip, subprocess
from TkinterDialog import TkinterDialog, Selections
from UsefulBioinformaticsFunctions import reverseCompliment, isPurine
from typing import IO


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

        ID = ''.join((self.chromosome,':',mutationPos0Base,'-',self.startPos,'(',strand,')'))

        return( '\t'.join((self.chromosome,mutationPos0Base,self.startPos,ID,mutation,strand)) )

    # Return a string that contains the mutation data in tab delimited MSIseq form.
    def formatForMSIseq(self):
        
        # Determine what type of mutation we have.
        if "substitution" in self.mutationType: MSIseqMutationType = "SNP"
        elif "insertion" in self.mutationType: MSIseqMutationType = "INS"
        elif "deletion" in self.mutationType: MSIseqMutationType = "DEL"
        else: raise ValueError("Error: Mutation type, \"" + self.mutationType + "\", does not appear to be a substitution, deletion, or insertion.")

        # Piece together the MSIseq formatted line and then return it.
        return( '\t'.join((self.chromosome,self.startPos,self.endPos,MSIseqMutationType,self.donorID)) )


# This class takes an ICGCFile object and, when iterated over, returns exactly once each mutation present in a donor that is
# the result of whole genome sequencing using GRCh37 as the reference genome.
# The mutation is returned as an ICGCMutation object.
class ICGCIterator:

    def __init__(self, ICGCFile: IO):
        
        self.ICGCFile = ICGCFile # The file containing the ICGC data
        self.finishedDonors = list() # A list of donors to make sure we don't encounter one more than once.
        self.currentDonor = '' # The donorID currently being read for mutation data.
        self.mutations = dict() # A list of mutations unique to the current donor, to avoid writing duplicate mutations.

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
                choppedUpLine[33] == "WGS"): # Was whole genome sequencing used to generate the data?

                # Is this a new donor?
                if choppedUpLine[1] != self.currentDonor:
                    self.finishedDonors.append(self.currentDonor)
                    self.currentDonor = choppedUpLine[1]
                    self.mutations.clear()
                    if self.currentDonor in self.finishedDonors:
                        raise ValueError("Error:  Donor " + self.currentDonor + " is present in more than one block of data!")
                    print("Reading and writing from donor",self.currentDonor)

                # Have we seen this mutation in this donor?
                if choppedUpLine[0] not in self.mutations:
                    self.mutations[choppedUpLine[0]] = None
                    return ICGCMutation(choppedUpLine)


#Create the Tkinter UI
dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(__file__),"..","data"))
dialog.createMultipleFileSelector("ICGC Mutation Files:",0,("gzip files",".gz"))
dialog.createCheckbox("Convert to Bed SNPs", 1, 0, 2)
dialog.createCheckbox("Convert to MSIseq format", 1, 2, 2)
dialog.createReturnButton(2,0,2)
dialog.createQuitButton(2,2,2)

# Run the UI
dialog.mainloop()

# If no input was received (i.e. the UI was terminated prematurely), then quit!
if dialog.selections is None: quit()

# Get the user's input from the dialog.
selections: Selections = dialog.selections
ICGCFilePaths = list(selections.getFilePathGroups())[0] # A list of ICGC mutation file paths
convertToBed = list(selections.getToggleStates())[0]
convertToMSIseq = list(selections.getToggleStates())[1]
if not convertToBed or convertToMSIseq: raise ValueError("Error: no output format selected.")

# Run the parser for each ICGC file given.
for ICGCFilePath in ICGCFilePaths:

    print("\nWorking in:",os.path.split(ICGCFilePath)[1])
    if not str(ICGCFilePath).endswith(".gz"):
        raise ValueError("Error:  Expected gzipped file.")
    if not "simple_somatic_mutation" in os.path.split(ICGCFilePath)[1]:
        raise ValueError("Error:  Expected file with \"simple_somatic_mutation\" in the name.")
    

    # Create a directory for intermediate files if it does not already exist...
    if not os.path.exists(os.path.join(os.path.dirname(ICGCFilePath),"intermediate_files")):
        os.mkdir(os.path.join(os.path.dirname(ICGCFilePath),"intermediate_files"))

    # Generate paths to output files.
    mutationGroupName = os.path.split(ICGCFilePath)[1].rsplit('.',3)[-3]

    if convertToBed:
        unsortedBedFilePath = os.path.join(os.path.dirname(ICGCFilePath),"intermediate_files",
            mutationGroupName+"_unsorted_singlenuc_context.bed")
        sortedBedFilePath = os.path.join(os.path.dirname(ICGCFilePath),mutationGroupName+"_singlenuc_context.bed")
    
    if convertToMSIseq:
        unsortedMSIseqDataFilePath = os.path.join(os.path.dirname(ICGCFilePath),"intermediate_files",
            mutationGroupName+"_unsorted_MSIseq_data.bed")
        sortedMSIseqDataFilePath = os.path.join(os.path.dirname(ICGCFilePath),mutationGroupName+"_MSIseq_data.bed")

    # Parse the files!
    with gzip.open(ICGCFilePath, 'r') as ICGCFile:

        if convertToBed: unsortedBedFile = open(unsortedBedFilePath, 'w')
        if convertToMSIseq: unsortedMSIseqDataFile = open(unsortedMSIseqDataFilePath, 'w')

        for mutation in ICGCIterator(ICGCFile):
            if convertToBed and mutation.isValidForBed(): unsortedBedFile.write(mutation.formatForBed() + '\n')   
            if convertToMSIseq: unsortedMSIseqDataFile.write(mutation.formatForMSIseq() + '\n')

        if convertToBed: unsortedBedFile.close()
        if convertToMSIseq: unsortedMSIseqDataFile.close()

    # Sort the mutationData using linux sort, because it's not quite so greedy with memory usage...
    if convertToBed:
        print("Sorting bed data...")
        subprocess.run(" ".join(("sort","-k1,1","-k2,2n",unsortedBedFilePath,">",sortedBedFilePath)), shell = True, check = True)
    if convertToMSIseq:
        print("Sorting MSIseq data...")
        subprocess.run(" ".join(("sort","-k5,5","-k2,2n",unsortedMSIseqDataFilePath,">",sortedMSIseqDataFilePath)), shell = True, check = True)