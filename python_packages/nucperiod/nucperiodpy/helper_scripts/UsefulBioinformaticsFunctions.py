# Various functions useful when working with (usually DNA related) bioinformatics data.

import subprocess
from typing import IO, List


# A dictionary that converts from one base to its reverse compliment
reverser = {'A':'T','T':'A','G':'C','C':'G','N':'N',
            'a':'t','t':'a','g':'c','c':'g','n':'n'}

def reverseCompliment(DNA):

    #Reverse the string using a nifty map!
    reverseComplimentMap = map(lambda base: reverser[base], DNA[::-1])
    reverseComplimentList = list(reverseComplimentMap)
    reverseCompliment = ''.join(reverseComplimentList)
    return reverseCompliment

#Determines whether or not a given base is a purine.
def isPurine(nucleotide: str) -> bool:
    return (nucleotide.upper() == "G" or nucleotide.upper() == "A")


#Returns a list of sequence identifiers from a fasta name resulting from the bedtools "getfasta" command.  
#The returned list contains chromosome, start pos, end pos, and strand, if present, in that order)
def parseFastaDescription(fastaSequenceName: str):
    """
    Returns a list of sequence identifiers from a fasta name resulting from the bedtools "getfasta" command.  
    The returned list contains chromosome, start pos, end pos, and strand, if present, in that order)
    """

    # Remove the leading '>' if necessary.
    if fastaSequenceName.startswith('>'): fastaSequenceName = fastaSequenceName[1:]

    # Get the chromosome. 
    # (The split result is taken from the rear in case there is a leading bedtools formatted name, which is separated by '::')
    splitSequence = fastaSequenceName.split(':')
    chromosome = splitSequence[-2]
    theRest = splitSequence[-1]

    # Get the strand. 'theRest' should look something like: 123-456(+) or 123-456() or 123-456
    if not '(' in theRest or theRest[-2] == '(': strand = None
    else: 
        strand = theRest[-2]
        if not strand in ('+','-','.'):
            raise ValueError("Error.  Unexpected symbol \'" + strand + "\' found where strand designation was expected.")
    theRest = theRest.split('(')[0]

    # Get the nucleotide positions. ('theRest' should look something like: 123-456)
    splitSequence = theRest.split('-')
    startPos = splitSequence[0]
    endPos = splitSequence[1]

    return (chromosome, startPos, endPos, strand)


# Parses fasta files one entry at a time.
# Designed to work with output from the bedtools getfasta function.
class FastaFileIterator:
    """
    Parses fasta files one entry at a time.
    Designed to work with output from the bedtools getfasta function.
    """


    # The object to hold information on each entry.
    class FastaEntry:
        def __init__(self, sequenceLocation: List[str], sequence: str, sequenceName: str):
            self.sequenceLocation = sequenceLocation # The list containing all the relevant information to locate 
                                                     # The fasta sequence in a given genome.
            self.sequenceName = sequenceName # The full fasta sequence header, without the angle bracket.

            # The information relevant to finding the sequence in the genome.
            self.chromosome = sequenceLocation[0]
            self.startPos = sequenceLocation[1]
            self.endPos = sequenceLocation[2]
            self.strand = sequenceLocation[3]

            self.sequence = sequence # The DNA sequence itself


    # Initialize the FastaFileIterator with an open fasta file object.
    def __init__(self, fastaFile: IO, containsLocationInformation = True):
        self.fastaFile = fastaFile # The file that will be parsed and read through.
        self.containsLocationInformation = containsLocationInformation # Whether or not the fasta file contains full
                                                                       # sequence location information
        self.eof = False # A flag for when the end of the file has been reached.

        # Initialize the iterator with the first entry.
        self.nextSequenceName = self.fastaFile.readline().strip()[1:]
        # If the file is empty (and thus, the first line is blank) set the eof flag to true.
        if not self.nextSequenceName: self.eof = True


    # Reads in the next fasta entry and returns it.
    def readEntry(self):

        # Return None if we've reached the end of the file.
        if self.eof: return None

        # Get the name for the upcoming sequence.
        sequenceName = self.nextSequenceName

        # Get the sequence location for the current entry (if available).
        if self.containsLocationInformation: 
            sequenceLocation = parseFastaDescription(sequenceName)
        else: sequenceLocation = (None,None,None,None)

        # Read through lines until we get to the next entry, adding to the sequence as we go.
        line = self.fastaFile.readline().strip()
        sequence = ''

        while not line.startswith('>'):

            sequence += line.upper()

            line = self.fastaFile.readline().strip()

            if not line:
                self.eof = True
                break

        # Prep for the next entry if we aren't at eof.
        if not self.eof:
            self.nextSequenceName = line[1:]

        # Return the current fasta entry.
        return self.FastaEntry(sequenceLocation, sequence, sequenceName)


    # Make the class iteratable, returning each fasta entry one at a time.
    def __iter__(self):
        return self
    def __next__(self):

        # Make sure we haven't reached the end of the file.
        if self.eof: raise StopIteration

        return self.readEntry()



def bedToFasta(bedFilePath, genomeFilePath, fastaOutputFilePath, 
               incorporateBedName = False, includeStrand = True, verbose = False):
    "Uses bedtools to convert a bed file to fasta format."

    if verbose:
        print("Calling shell subprocess to use bedtools to generate a fasta file from the given bed file...")

    optionalParameters = list()
    if includeStrand: optionalParameters.append("-s")
    if incorporateBedName: optionalParameters.append("-name")

    subprocess.run(" ".join(("bedtools","getfasta") + tuple(optionalParameters) + ("-fi",genomeFilePath,
        "-bed",bedFilePath,">",fastaOutputFilePath)), shell = True, check = True)