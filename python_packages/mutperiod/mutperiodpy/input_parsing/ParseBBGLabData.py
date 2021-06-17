# Parses data from the bbglab's format.  More of a proof of concept and mostly deprecated, tbh.

import os, gzip
from benbiohelpers.DNA_SequenceHandling import reverseCompliment, isPurine


# This function takes the input from the bbglab's mutation data and converts it to bed format.
# Inputs are (respectively) the directory containing the bbglab data and the directory the bed data will be created in.
def convertDirectoryToBedFiles(preBedDirectory, bedDirectory):
    
    # Iterate through the given directory
    for fileName in os.listdir(preBedDirectory):

        # Copy the directory structure from the bbglab data to the bed data.
        if os.path.isdir(preBedDirectory + "/" + fileName):

            print ("Found directory: " + fileName)

            if not os.path.exists(bedDirectory + "/" + fileName):
                 os.mkdir(bedDirectory + "/" + fileName)
                
            convertDirectoryToBedFiles(preBedDirectory + "/" + fileName,
            bedDirectory + "/" + fileName)

        # Send gzipped files to the copyBedData function to be converted
        elif fileName.endswith("tsv.gz"): 
            copyBedData(preBedDirectory, bedDirectory, fileName)
        
        # Error case if extra files are included that shouldn't be.
        else: 
            print ("Error:  Unexpected file type")

# This function converts a single mutation file from bbglab format to bed format
def copyBedData(preBedDirectory, bedDirectory, fileName):

    print ("Reading from " + fileName + ".")

    # A list of all of the individual mutations in bed format
    bedData = list()

    # Use gzip to read the file contents and generate a bed formatted output.
    with gzip.open(preBedDirectory+"/"+fileName, "r") as preBedData:
        
        # Each line contains one mutation to be converted and added to bedData.
        for mutationData in preBedData.readlines():
            
            # The start of the bed formatted data.
            bedMutation = ""

            # Convert the bbglab data into a list of data entries.
            preBedDataCols = list()
            for data in mutationData.strip().split():
                preBedDataCols.append(str(data,"utf-8"))

            # Construct the bed data format from the bbglab data.
            bedMutation += preBedDataCols[0] + "\t" # Chromosome number
            bedMutation += str(int(preBedDataCols[1])-1) + "\t"  # base-0 start (hence the "-1")
            bedMutation += preBedDataCols[1] + "\t" # base-1 end
            
            # Based on the nature of the mutation, asign it to either the + or - strand and output the mutation accordingly.
            # Mutations are assumed to have arisen in pyrimidines.
            if isPurine(preBedDataCols[2]):
                bedMutation += reverseCompliment(preBedDataCols[2]) + '\t'
                bedMutation += reverseCompliment(preBedDataCols[2]) + ">" + reverseCompliment(preBedDataCols[3]) + "\t"
                bedMutation += "-\n"
            else:
                bedMutation += preBedDataCols[2] + '\t'
                bedMutation += preBedDataCols[2] + ">" + preBedDataCols[3] + "\t"
                bedMutation += "+\n"

            # Add the mutation entry to the list of bed data.
            bedData.append(bedMutation)

    # Write (gzipped) the bed formatted mutation data to the bed directory.
    print ("Writing bed data.")
    with gzip.open(bedDirectory + "/" + fileName, "w") as bedFile:
        for data in bedData:
            bedFile.write(data.encode())

#Get the current working directory.
workingDirectory = os.path.dirname(os.path.realpath(__file__))

#Designate the directories with the data to be converted to bed format, and the location for the bed data to be written.
preBedDirectory = os.path.join(workingDirectory,'..') + "/bbglab_data"
bedDirectory = os.path.join(workingDirectory,'..') + "/bed_data"

# Go! Go! Go!
convertDirectoryToBedFiles(preBedDirectory,bedDirectory)