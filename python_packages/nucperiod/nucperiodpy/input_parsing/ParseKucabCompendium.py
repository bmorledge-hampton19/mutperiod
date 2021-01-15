# This script takes data from the Kucab et al. mutation compendium paper and converts it to
# a trinucleotide context bed file.

import os, subprocess
from nucperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog, Selections
from nucperiodpy.helper_scripts.UsefulBioinformaticsFunctions import reverseCompliment, isPurine
from nucperiodpy.helper_scripts.UsefulFileSystemFunctions import (getIsolatedParentDir, generateFilePath, dataDirectory,
                                                                  dataTypes, generateMetadata, getAcceptableChromosomes)
                                                                  

def parseKucabCompendium(kucabSubstitutionsFilePaths, genomeFilePath, nucPosFilePath, includeAllPAHs):

    for kucabSubstitutionsFilePath in kucabSubstitutionsFilePaths:

        print("\nWorking in:",os.path.basename(kucabSubstitutionsFilePath))

        if not kucabSubstitutionsFilePath.endswith("final.txt"):
            raise ValueError("Error:  Expected input file from Kucab data which should end in \"final.txt\".")

        # Prepare the output file path.
        localRootDirectory = os.path.dirname(kucabSubstitutionsFilePath)
        dataGroupName = getIsolatedParentDir(kucabSubstitutionsFilePath)
        if includeAllPAHs:
            outputDirectory = os.path.join(localRootDirectory,"all_PAHs")
            dataGroupName += "_all_PAHs"
        else: 
            dataGroupName += "_smoker_lung"
            outputDirectory = os.path.join(localRootDirectory,"smoker_lung")
        
        # Make sure the data directory exists.
        if not os.path.exists(outputDirectory): os.mkdir(outputDirectory)

        # Generate the output file path and metadata
        outputTrinucBedFilePath = generateFilePath(directory = outputDirectory, dataGroup = dataGroupName,
                                                   context = "trinuc", dataType = dataTypes.mutations, fileExtension = ".bed")
        generateMetadata(dataGroupName, getIsolatedParentDir(genomeFilePath), getIsolatedParentDir(nucPosFilePath),
                         os.path.join("..",os.path.basename(kucabSubstitutionsFilePath)), outputDirectory)

        # Get the list of acceptable chromosomes
        acceptableChromosomes = getAcceptableChromosomes(genomeFilePath)

        # These are the designations for PAH mutation signatures, the ones related to tobacco smoke that we want to study.
        PAHDesignations = ("MSM0.54","MSM0.26","MSM0.92","MSM0.2","MSM0.42","MSM0.74","MSM0.103"
                           "MSM0.14","MSM0.82","MSM0.130","MSM0.12","MSM0.132","MSM0.13","MSM0.96")
        # These designations specifically mimic the indel signature in smokers' lung cancer tumors.
        LungCancerSpecificDesignations = ("MSM0.26","MSM0.92","MSM0.2","MSM0.103","MSM0.14")

        # Set the designations that will be used to collect data based on the input to the function.
        if includeAllPAHs:
            relevantDesignations = PAHDesignations
        else: relevantDesignations = LungCancerSpecificDesignations

        print("Reading data and writing to trinuc bed file...")
        with open(kucabSubstitutionsFilePath, 'r') as kucabSubstitutionsFile:
            with open(outputTrinucBedFilePath, 'w') as outputTrinucBedFile:

                firstLineFlag = True
                for line in kucabSubstitutionsFile:
                    
                    # Skip the first line with headers.
                    if firstLineFlag:
                        firstLineFlag = False
                        continue

                    # The lines are separated by tabs.  The relevant data have the following indices in a tab-separated list:
                    # 15: mutagen designation
                    # 4: Chromosome
                    # 5: Start Pos (1 base)
                    # 6: Reference base
                    # 7: Mutated base
                    # 13: pre-base context
                    # 14: post-base context
                    choppedUpLine = line.strip().split('\t')

                    # Skip the mutation if it does not belong to the relevant group.
                    if not choppedUpLine[15] in relevantDesignations: continue

                    # Compile the necessary information for the bed file.
                    chromosome = "chr" + choppedUpLine[4]

                    # Handle the weird chromsome formatting and then check for invalid chromosomes.
                    if chromosome == "chr23": chromosome = "chrX"
                    if chromosome == "chr24": chromosome = "chrY"
                    if not chromosome in acceptableChromosomes: continue
                    startPos1Base = choppedUpLine[5]
                    startPos0Base = str(int(startPos1Base)-1)

                    mutatedFrom = choppedUpLine[6]
                    mutatedTo = choppedUpLine[7]
                    trinucContext = ''.join((choppedUpLine[13],mutatedFrom,choppedUpLine[14]))

                    # If the mutated base is listed as arising from a purine, flip the mutation and the strand.
                    if isPurine(mutatedFrom):
                        mutation = reverseCompliment(mutatedFrom) + '>' + reverseCompliment(mutatedTo)
                        strand = '-'
                        trinucContext = reverseCompliment(trinucContext)
                    else:
                        mutation = mutatedFrom + '>' + mutatedTo
                        strand = '+'

                    # Write the information to the trinuc bed file.
                    outputTrinucBedFile.write('\t'.join((chromosome, startPos0Base, startPos1Base,
                                                         trinucContext, mutation, strand)) + '\n')

        # Sort the output file.
        print("Sorting output file...")
        subprocess.run(" ".join(("sort","-k1,1","-k2,2n",outputTrinucBedFilePath,"-o",outputTrinucBedFilePath)),
                           shell = True, check = True)


if __name__ == "__main__":

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=dataDirectory)
    dialog.createMultipleFileSelector("Kucab Substitutions File Paths:",0,"final.txt",("text files",".txt")) #NOTE: Weird file ending?
    dialog.createFileSelector("Genome Fasta File:",1,("Fasta Files",".fa"))
    dialog.createFileSelector("Strongly Positioned Nucleosome File:",2,("Bed Files",".bed"))
    dialog.createCheckbox("Include all PAH Designations",3,0)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    kucabSubstitutionsFilePaths = list(selections.getFilePathGroups())[0]
    genomeFilePath = list(selections.getIndividualFilePaths())[0]
    nucPosFilePath = list(selections.getIndividualFilePaths())[1]
    includeAllPAHs = list(selections.getToggleStates())[0]

    parseKucabCompendium(kucabSubstitutionsFilePaths, genomeFilePath, nucPosFilePath, includeAllPAHs)