# This script takes paths to files containing nucleosome counts and exports them to an R script which generates
# some nice plots for the files and exports them to a given location.

import os, subprocess, sys
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory, rScriptsDirectory, checkDirs


def generateFigures(tsvFilePaths, rdaFilePaths, exportPath, omitOutliers, smoothNucGroup, 
                    includeNorm, includeRaw, strandAlign):

    # Check for invalid arguments.
    assert len(rdaFilePaths) == 0 or (includeNorm or includeRaw), ("at least one rda file was included, but " +
                                                                    "neither normalized nor raw counts are chosen to be used.")

    assert os.path.isdir(exportPath) or exportPath.endswith(".pdf"), ("The export path: " + exportPath + " is neither a directory nor a pdf file.")

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
        inputsFile.write('$'.join(rdaFilePaths) + '\n')
        inputsFile.write(exportDir + '\n')
        inputsFile.write(exportFileName + '\n')

    # Call the R script to generate the figures.
    print("Calling R script...")
    subprocess.run(("Rscript",os.path.join(rScriptsDirectory,"GenerateFigures.R"),inputsFilePath, str(omitOutliers),
                    str(smoothNucGroup), str(includeNorm), str(includeRaw), str(strandAlign)), check = True)


def parseArgs(args):
    
    # If only the subcommand was given, run the UI.
    if len(sys.argv) == 2: 
        main(); return

    # Format the arguments.
    if args.tsv_paths is None: args.tsv_paths = list()
    if args.rda_paths is None: args.rda_paths = list()

    if args.output_directory is not None: exportPath = args.output_directory
    elif args.output_file is not None: exportPath = args.output_file
    else: raise ValueError("No output path given.")

    # Pass the given commands to the generateFigures function
    generateFigures(args.tsv_paths, args.rda_paths, exportPath, args.omit_outliers, args.smooth_nuc_group,
                    args.include_normalized, args.include_raw, args.align_strands)


def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createMultipleFileSelector("Nucleosome Counts Files:",0,
                                      "nucleosome_mutation_counts.tsv",("tsv files",".tsv"))
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
    dialog.createCheckbox("Use normalized values from rda input", 4, 0)
    dialog.createCheckbox("Use raw values from rda input", 4, 1)
    dialog.createCheckbox("Strand align results", 5, 0)

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
    includeNorm = bool(selections.getToggleStates()[2])
    includeRaw = bool(selections.getToggleStates()[3])
    strandAlign = bool(selections.getToggleStates()[4])

    generateFigures(tsvFilePaths, rdaFilePaths, exportPath, omitOutliers, smoothNucGroup, 
                    includeNorm, includeRaw, strandAlign)

if __name__ == "__main__": main()