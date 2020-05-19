# The class to hold the data resulting from the periodicity and asymmetry analyses.
NucPeriodData = setClass("NucPeriodData", slots = list(periodicityResults = "data.table",
                                                       asymmetryResults = "data.table",
                                                       MSIInputs = "char", MSSInputs = "char", wilcoxinResult = "list"))

generateNucPeriodData = function(rawCountsFilePaths, MSIFilePaths = NA, MSSFilePaths = NA) {



}
