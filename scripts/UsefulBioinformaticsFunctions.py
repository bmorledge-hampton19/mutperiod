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
    return (nucleotide.lower() == "g" or nucleotide.lower() == "a")