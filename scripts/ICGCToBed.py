# This script reads one or more "simple somatic mutation" data file(s) from ICGC and 
# writes information on single base substitution mutations to a new bed file or files for further analysis.
import os, gzip, subprocess
from TkinterDialog import TkinterDialog, Selections
from UsefulBioinformaticsFunctions import reverseCompliment, isPurine

# This is the function that does the actual parsing from icgc format to bed format.
def ICGCToBed(ICGCFilePath, bedFilePath):

    mutationData = list() # A list of mutation data, of which each entry is a bed formatted list of data.
    DNABases = ('A','C','T','G') # A list of single DNA nucleotides to check against.
    currentDonor = '' # The donorID currently being read for mutation data.
    finishedDonors = list()
    mutations = dict() # A list of mutations unique to the current donor, to avoid writing duplicate mutations.

    # Access the files.
    with gzip.open(ICGCFilePath, 'r') as ICGCFile:
        with open(bedFilePath, 'w') as bedFile:

            # Skip the header line in the icgc file.
            ICGCFile.readline()

            print("Reading information from ICGC file and writing to bed file...")

            # Read through the icgc file one line at a time.
            for line in ICGCFile:

                # Split the line into its individual data components.  The relevant lines are:
                #   1: donor ID
                #   8: chromosome ID
                #   9: mutation start position
                #   15: mutated from
                #   16: mutated to
                #   11: strand
                #   12: reference genome version
                #   33: sequencing method
                choppedUpLine = str(line,"utf-8").strip().split('\t')

                if choppedUpLine[11] != "1":
                    raise ValueError("Error.  Strand field does not contain \"1\" for plus strand.  " +
                                    "Found " + str(choppedUpLine[11]) + " instead.")

                # Make sure the given mutation meets the criteria for incorporation into the bed file.
                if (choppedUpLine[15].upper() in DNABases and choppedUpLine[16].upper() in DNABases and # Is it a single base substitution?
                    choppedUpLine[12] == "GRCh37" and # Is the reference genome hg19?
                    choppedUpLine[33] == "WGS"): # Was whole genome sequencing used to generate the data?

                    # Is this a new donor?
                    if choppedUpLine[1] != currentDonor:
                        finishedDonors.append(currentDonor)
                        currentDonor = choppedUpLine[1]
                        if currentDonor in finishedDonors:
                            raise ValueError("Error:  Donor " + currentDonor + " is present in more than one block of data!")
                        # Write this donor's mutation data to the bed file.
                        for mutationID in mutations: bedFile.write(mutations[mutationID])
                        mutations.clear()
                        print("Reading and writing from donor",currentDonor)

                    # Extract the relevant information for the bed file.

                    # Construct a unique mutation ID.
                    mutationID = ("C"+choppedUpLine[8]+"M"+choppedUpLine[9]+choppedUpLine[15]+choppedUpLine[16])
                    if mutationID in mutations: continue

                    chromosome = "chr" + choppedUpLine[8]
                    mutationPos1Base = choppedUpLine[9]

                    mutationPos0Base = str(int(mutationPos1Base)-1)

                    # If the mutated base is listed as arising from a purine, flip the mutation and the strand.
                    if isPurine(choppedUpLine[15]):
                        mutation = reverseCompliment(choppedUpLine[15]) + '>' + reverseCompliment(choppedUpLine[16])
                        strand = '-'
                    else:
                        mutation = choppedUpLine[15] + '>' + choppedUpLine[16]
                        strand = '+'

                    # Create the ID for the nucleotide coordinate.
                    ID = ''.join((chromosome,':',mutationPos0Base,'-',mutationPos1Base,'(',strand,')'))

                    # Add this mutation to the dictionary for the current donor.
                    mutations[mutationID] = '\t'.join((chromosome,mutationPos0Base,mutationPos1Base,ID,mutation,strand)) + '\n'

            # Don't forget the last donor!
            for mutationID in mutations: bedFile.write(mutations[mutationID])


#Create the Tkinter UI
dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(__file__),"..","data"))
dialog.createMultipleFileSelector("ICGC Mutation Files:",0,("gzip files",".gz"))
dialog.createReturnButton(1,0,2)
dialog.createQuitButton(1,2,2)

# Run the UI
dialog.mainloop()

# If no input was received (i.e. the UI was terminated prematurely), then quit!
if dialog.selections is None: quit()

# Get the user's input from the dialog.
selections: Selections = dialog.selections
ICGCFilePaths = list(selections.getFilePathGroups())[0] # A list of ICGC mutation file paths

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

    # Generate a path to the output file.
    mutationGroupName = os.path.split(ICGCFilePath)[1].rsplit('.',3)[-3]
    unsortedBedFilePath = os.path.join(os.path.dirname(ICGCFilePath),"intermediate_files",
        mutationGroupName+"_unsorted_singlenuc_context.bed")
    sortedBedFilePath = os.path.join(os.path.dirname(ICGCFilePath),mutationGroupName+"_singlenuc_context.bed")

    # Parse the file!
    ICGCToBed(ICGCFilePath, unsortedBedFilePath)

    # Sort the mutationData using linux sort, because it doesn't absolutely destroy memory usage...
    print("Sorting bed data...")
    subprocess.run(" ".join(("sort","-k1,1","-k2,2n",unsortedBedFilePath,">",sortedBedFilePath)), shell = True, check = True)