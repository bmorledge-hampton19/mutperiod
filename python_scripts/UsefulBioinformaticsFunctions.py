import subprocess

#Create a dictionary that converts from one base to its reverse compliment
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

def bedToFasta(bedFilePath, genomeFilePath, fastaOutputFilePath):
    "Uses bedtools to convert a bed file to fasta format."

    print("Calling shell subprocess to use bedtools to generate a fasta file from the given bed file...")

    if incorporateBedName:
        subprocess.run(" ".join(["bedtools","getfasta","-s","-name","-fi",genomeFilePath,
            "-bed",bedFilePath,">",fastaOutputFilePath]), shell = True, check = True)
    else:
        subprocess.run(" ".join(["bedtools","getfasta","-s","-fi",genomeFilePath,
            "-bed",bedFilePath,">",fastaOutputFilePath]), shell = True, check = True)