# This script takes paths to files containing nucleosome counts and exports them to an R script which generates
# some nice plots for the files and exports them to a given location.

import os, subprocess, sys
from typing import List
from benbiohelpers.CustomErrors import UserInputError, InvalidPathError, checkIfPathExists
from benbiohelpers.TkWrappers.TkinterDialog import TkinterDialog
from benbiohelpers.FileSystemHandling.DirectoryHandling import getFilesInDirectory
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory, getExpectedPeriod, rScriptsDirectory, DataTypeStr


def generateFigures(tsvFilePaths: List[str], rdaFilePaths: List[str], exportPath, omitOutliers, smoothNucGroup, strandAlign):

    # Check for invalid arguments.
    if len(rdaFilePaths) + len(tsvFilePaths) == 0:
        raise UserInputError("No input files were found to generate graphs from.")

    if not (os.path.isdir(exportPath) or exportPath.endswith(".pdf")): 
        raise InvalidPathError("The given export path is neither a directory nor a pdf file.")

    for tsvFilePath in tsvFilePaths:
        if not tsvFilePath.endswith(DataTypeStr.generalNucCounts + ".tsv"):
            raise InvalidPathError(tsvFilePath, "The given nucleosome counts tsv file does not end with \"" + 
                                   DataTypeStr.generalNucCounts + ".tsv\" as expected.")

    for rdaFilePath in rdaFilePaths:
        if not rdaFilePath.endswith(".rda"):
            raise InvalidPathError(rdaFilePath, "The given nucleosome counts rda file does not end with \".rda\" as expected.")

    # Retrieve the expected periods for each of the given tsv counts files.
    tsvExpectedPeriods = [str(getExpectedPeriod(tsvFilePath)) for tsvFilePath in tsvFilePaths]

    # Determine whether the export path is a directory or file and set variables accordingly.
    if os.path.isdir(exportPath):
        exportDir = exportPath
        exportFileName = ''
    else:
        exportDir = os.path.dirname(exportPath)
        exportFileName = os.path.basename(exportPath)

    # Create the temporary inputs file to pass to the R script
    inputsFilePath = os.path.join(os.getenv("HOME"), ".mutperiod","R_inputs.txt")

    # Write the inputs
    with open(inputsFilePath, 'w') as inputsFile:
        inputsFile.write('$'.join(tsvFilePaths) + '\n')
        inputsFile.write('$'.join(tsvExpectedPeriods) + '\n')
        inputsFile.write('$'.join(rdaFilePaths) + '\n')
        inputsFile.write(exportDir + '\n')
        inputsFile.write(exportFileName + '\n')

    # Call the R script to generate the figures.
    print("Calling R script...")
    subprocess.run(("Rscript",os.path.join(rScriptsDirectory,"GenerateFigures.R"),inputsFilePath, str(omitOutliers),
                    str(smoothNucGroup), str(strandAlign)), check = True)


def parseArgs(args):
    
    # If only the subcommand was given, run the UI.
    if len(sys.argv) == 2: 
        main(); return

    # Get valid tsv paths from the given input.
    tsvFilePaths = list()
    if args.tsv_paths is not None:
        for tsvFilePath in args.tsv_paths:
            checkIfPathExists(tsvFilePath)
            if os.path.isdir(tsvFilePath):
                tsvFilePaths += [os.path.abspath(filePath) for filePath in getFilesInDirectory(tsvFilePath, DataTypeStr.generalNucCounts + ".tsv")]
            else: tsvFilePaths.append(os.path.abspath(tsvFilePath))
    
    # Get valid rda paths from the given input.
    rdaFilePaths = list()
    if args.rda_paths is not None:
        for rdaFilePath in args.rda_paths:
            checkIfPathExists(rdaFilePath)
            if os.path.isdir(rdaFilePath):
                rdaFilePaths += [os.path.abspath(filePath) for filePath in getFilesInDirectory(rdaFilePath, ".rda")]
            else: rdaFilePaths.append(os.path.abspath(rdaFilePath))

    if args.output_directory is not None: exportPath = os.path.abspath(args.output_directory)
    elif args.output_file is not None: exportPath = os.path.abspath(args.output_file)
    else: raise UserInputError("No output path given.")

    # Pass the given commands to the generateFigures function
    generateFigures(list(set(tsvFilePaths)), list(set(rdaFilePaths)), exportPath, args.remove_outliers, args.smooth_nuc_group, args.align_strands)


def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Nucleosome Counts Files:",0,
                                      DataTypeStr.generalNucCounts + ".tsv",("tsv files",".tsv"))
    dialog.createMultipleFileSelector("R Nucleosome Mutation Analysis Files:",1,
                                      ".rda",("rda files",".rda"))

    fileNumSelector = dialog.createDynamicSelector(2,0)
    fileNumSelector.initCheckboxController("Export to one file (as opposed to one file for each graph)")
    oneFileDialog = fileNumSelector.initDisplay(1, "oneFile")
    oneFileDialog.createFileSelector("Export File", 0, ("pdf file",".pdf"), newFile = True)
    manyFilesDialog = fileNumSelector.initDisplay(0, "manyFiles")
    manyFilesDialog.createFileSelector("Export Directory", 0, directory = True)
    fileNumSelector.initDisplayState()

    dialog.createCheckbox("Omit Outliers", 3, 0)
    dialog.createCheckbox("Smooth Nuc Group results", 3, 1)
    dialog.createCheckbox("Strand align results", 4, 0)
    #dialog.createCheckbox("Use normalized values from rda input", 5, 0)
    #dialog.createCheckbox("Use raw values from rda input", 5, 1)

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
    smoothNucGroup = bool(selections.getToggleStates()[1])
    strandAlign = bool(selections.getToggleStates()[2])

    generateFigures(tsvFilePaths, rdaFilePaths, exportPath, omitOutliers, smoothNucGroup, 
                    strandAlign)

if __name__ == "__main__": main()