# This script takes paths to files containing nucleosome counts and exports them to an R script which generates
# some nice plots for the files and exports them to a given location.

import os, subprocess
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import dataDirectory, RPackagesDirectory


def generateFigures(tsvFilePaths, rdaFilePaths, exportPath, omitOutliers):

    # Determine whether the export path is a directory or file and set variables accordingly.
    if os.path.isdir(exportPath):
        exportDir = exportPath
        exportFileName = ''
    else:
        exportDir = os.path.dirname(exportPath)
        exportFileName = os.path.basename(exportPath)

    # Create the temporary inputs file to pass to the R script
    rScriptDirectory = os.path.join(RPackagesDirectory,"RunNucPeriod")
    inputsFilePath = os.path.join(rScriptDirectory,"inputs.txt")

    # Write the inputs
    with open(inputsFilePath, 'w') as inputsFile:
        inputsFile.write('$'.join(tsvFilePaths) + '\n')
        inputsFile.write('$'.join(rdaFilePaths) + '\n')
        inputsFile.write(exportDir + '\n')
        inputsFile.write(exportFileName + '\n')

    # Call the R script to generate the figures.
    print("Calling R script...")
    subprocess.run(" ".join(("Rscript",os.path.join(rScriptDirectory,"GenerateFigures.R"),inputsFilePath, str(omitOutliers))),
                   shell = True, check = True)


if __name__ == "__main__":

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=dataDirectory)
    dialog.createMultipleFileSelector("Nucleosome Counts Files:",0,
                                      "nucleosome_mutation_counts.tsv",("tsv files",".tsv"))
    dialog.createMultipleFileSelector("R Nucleosome Mutation Analysis Files:",1,
                                      ".rda",("rda files",".rda"))

    fileNumSelector = dialog.createDynamicSelector(2,0)
    fileNumSelector.initCheckboxController("Export to one file?")
    oneFileDialog = fileNumSelector.initDisplay(1, "oneFile")
    oneFileDialog.createFileSelector("Export File", 0, ("pdf file",".pdf"), newFile = True)
    manyFilesDialog = fileNumSelector.initDisplay(0, "manyFiles")
    manyFilesDialog.createFileSelector("Export Directory", 0, directory = True)
    fileNumSelector.initDisplayState()

    dialog.createCheckbox("Omit Outliers?", 3, 0)
    dialog.createExitButtons(4,0)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections = dialog.selections
    tsvFilePaths = selections.getFilePathGroups()[0]
    rdaFilePaths = selections.getFilePathGroups()[1]
    if fileNumSelector.getControllerVar(): exportPath = selections.getIndividualFilePaths("oneFile")[0]
    else: exportPath = selections.getIndividualFilePaths("manyFiles")[0]
    omitOutliers = bool(selections.getToggleStates()[0])

    generateFigures(tsvFilePaths, rdaFilePaths, exportPath, omitOutliers)
