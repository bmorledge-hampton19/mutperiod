# This script takes a given cohort group (e.g. microsatellite satellite stability) and groups it by a potential secondary cohort feature (e.g. mut sigs)
import os
from typing import Dict
from mutperiodpy.Tkinter_scripts.TkinterDialog import TkinterDialog, Selections
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import getDataDirectory


# VERY TEMPORARY EXECUTION.  Not very adaptable at all.
def groupCohorts(dataSetDirectory):

    dataSetName = os.path.basename(dataSetDirectory)

    # Assume we are grouping by MS stability and mut sigs.
    cohortMSStatus = dict()
    mSICohortsFilePath = os.path.join(dataSetDirectory,"microsatellite_analysis",
                                      dataSetName+"_MSI_cohorts.txt")
    
    cohortMutSig = dict()
    mutSigDesignationsFilePath = os.path.join(dataSetDirectory,"mut_sig_analysis",
                                              dataSetName+"_mut_sig_assignments.tsv")

    mutSigByMSStatus: Dict[str,Dict[str,int]] = dict()

    # Get cohort mut sig designations
    with open(mutSigDesignationsFilePath,'r') as mutSigDesignationsFile:

        # Ignore headers (Should probably remove these during generation?)
        mutSigDesignationsFile.readline()

        for line in mutSigDesignationsFile:
            choppedUpLine = line.strip().split('\t')
            cohortMutSig[choppedUpLine[0]] = choppedUpLine[1:]

    # Get microsatellite designations
    mSICohorts = list()
    with open(mSICohortsFilePath, 'r') as mSICohortsFile:

        for line in mSICohortsFile: mSICohorts.append(line.strip())

    for cohort in cohortMutSig:

        if cohort in mSICohorts: cohortMSStatus[cohort] = "MSI"
        else: cohortMSStatus[cohort] = "MSS"

    # Now, record the mutation signatures seen for each MS designation.
    mutSigByMSStatus["MSS"] = dict()
    mutSigByMSStatus["MSI"] = dict()

    for cohort in cohortMSStatus:
        mSStatus = cohortMSStatus[cohort]
        for mutSig in cohortMutSig[cohort]:
            mutSigByMSStatus[mSStatus][mutSig] = mutSigByMSStatus[mSStatus].setdefault(mutSig,0) + 1

    for mSStatus in mutSigByMSStatus:
        print(mSStatus)
        mSStatusCount = sum(mutSigByMSStatus[mSStatus][mutSig] for mutSig in mutSigByMSStatus[mSStatus])
        for mutSig in mutSigByMSStatus[mSStatus]:
            print('\t', mutSig, ": ", mutSigByMSStatus[mSStatus][mutSig], "/", mSStatusCount, sep = '')


def main():

    #Create the Tkinter UI
    dialog = TkinterDialog(workingDirectory=getDataDirectory())
    dialog.createFileSelector("Data Set Directory:",0,directory=True)
    
    # TODO: Maybe allow for different cohort selections here.  Right now, I'm assuming that we group MS data by mut sigs.
    #       This will probably work best if it communicates with the project manager to some capacity (what cohort data is available?)

    # Run the UI
    dialog.mainloop()

    # If no input was received (i.e. the UI was terminated prematurely), then quit!
    if dialog.selections is None: quit()

    # Get the user's input from the dialog.
    selections: Selections = dialog.selections
    dataSetDirectory = selections.getFilePaths()[0]

    groupCohorts(dataSetDirectory)


if __name__ == "__main__": main()