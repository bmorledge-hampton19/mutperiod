# This script takes normalized nucleosome mutation counts files and passes them to an R script
# which outputs relevant data about them such as periodicity snr, assymetry, and differences between MSI and MSS data.
import os, subprocess
from UsefulFileSystemFunctions import dataTypes
from TkinterDialog import TkinterDialog, Selections
from typing import List

def runNucleosomeMutationAnalysis(normalizedNucleosomeMutationCountsFilePaths: List[str], outputFilePath: str, 
                                  MSIFilePaths: List[str], MSSFilePaths: List[str]):

    # Check for potential value errors.
    if len(MSIFilePaths) * len(MSSFilePaths) == 0 and len(MSIFilePaths) + len(MSSFilePaths) != 0:
        raise ValueError("Either MSI or MSS files were given, but the other group was not.\n" +
                         "Both or neither should be present.")

    if len(normalizedNucleosomeMutationCountsFilePaths) == 0:
        raise ValueError("No normalized counts files given.")

    if not outputFilePath.endswith(".rda"): raise ValueError("Output file should end with \".rda\".")

    # Write the inputs to a temporary file to be read by the R script
    rScriptDirectory = os.path.join(os.path.dirname(__file__),"..","R_scripts","RunNucPeriod")
    inputsFilePath = os.path.join(rScriptDirectory,"inputs.txt")

    with open(inputsFilePath, 'w') as inputsFile:
        if (len(MSIFilePaths) == 0 and len(MSSFilePaths) == 0):
            print("Generating inputs to run analysis without microsatellite designation...")
            inputsFile.write('\n'.join(('$'.join(normalizedNucleosomeMutationCountsFilePaths), outputFilePath)) + '\n')
            
        else:
            print("Generating inputs to run analysis with microsatellite designation...")
            inputsFile.write('\n'.join(('$'.join(normalizedNucleosomeMutationCountsFilePaths), outputFilePath,
                                        '$'.join(MSIFilePaths), '$'.join(MSSFilePaths))) + '\n')

    # Call the R script
    print("Calling R script...")
    subprocess.run(" ".join(("Rscript",os.path.join(rScriptDirectory,"RunNucleosomeMutationAnalysis.R"),inputsFilePath)),
                   shell = True, check = True)

    print("Results can be found at",outputFilePath)


if __name__ == "__main__":

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=os.path.join(os.path.dirname(__file__),"..","data"))
    dialog.createMultipleFileSelector("Normalized Nucleosome Mutation Counts files:",0,
                                      dataTypes.normNucCounts + ".tsv",("Tab Seperated Values Files",".tsv"))
    dialog.createFileSelector("Output File", 1, ("R Data File", ".rda"), newFile = True)
    dialog.createMultipleFileSelector("MSI files (optional):",2,
                                      dataTypes.normNucCounts + ".tsv",("Tab Seperated Values Files",".tsv"))
    dialog.createMultipleFileSelector("MSS files (optional):",3,
                                      dataTypes.normNucCounts + ".tsv",("Tab Seperated Values Files",".tsv"))                                  
    dialog.createReturnButton(4,0,2)
    dialog.createQuitButton(4,2,2)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    normalizedNucleosomeMutationCountsFilePaths = list(selections.getFilePathGroups())[0]
    outputFilePath = list(selections.getIndividualFilePaths())[0]
    MSIFilePaths = list(selections.getFilePathGroups())[1]
    MSSFilePaths = list(selections.getFilePathGroups())[2]

    runNucleosomeMutationAnalysis(normalizedNucleosomeMutationCountsFilePaths, outputFilePath,
                                  MSIFilePaths, MSSFilePaths)