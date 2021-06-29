# This script will be called from the command line to execute other scripts.
from argparse import ArgumentParser
from mutperiodpy import (RunNucleosomeMutationAnalysis, RunAnalysisSuite, GenerateFigures)
from mutperiodpy.input_parsing import (ParseCustomBed, ParseICGC)
from mutperiodpy.helper_scripts.UsefulFileSystemFunctions import DataTypeStr
import argparse, importlib.util
if importlib.util.find_spec("shtab") is not None: 
        import shtab
        fileCompletion = shtab.FILE
else:
        fileCompletion = None

        


def formatParseICGCParser(parseICGCParser: ArgumentParser):

    parseICGCParser.set_defaults(func = ParseICGC.parseArgs)
    parseICGCParser.add_argument("ICGCFilePaths", nargs = '*',
                                 help = "One or more paths to ICGC files to parse.  Should be gzipped.  If given a directory, "
                                        "the directory will be recursively searched for files ending in \".tsv.gz\".").complete = fileCompletion

    parseICGCParser.add_argument("-g", "--genome-file", help = "The associated genome fasta file").complete = fileCompletion

    parseICGCParser.add_argument("-d", "-c", "--stratify-by-donors", action = "store_true", 
                                 help = "Stratify results by individual donors")
    parseICGCParser.add_argument("-m", "--stratify-by-Microsatellite", action = "store_true", 
                                 help = "Stratify results by microsatellite stability")
    parseICGCParser.add_argument("-s", "--stratify-by-Mut-Sigs", action = "store_true", 
                                 help = "Stratify results by mutation signature")


def formatParseBedParser(parseBedParser: ArgumentParser):

    parseBedParser.set_defaults(func = ParseCustomBed.parseArgs)
    parseBedParser.add_argument("bedFilePaths", nargs = '*',
                                help = "One or more paths to bed files to parse.  If given a directory, "
                                       "the directory will be recursively searched for files ending in \"custom_input.bed\".").complete = fileCompletion

    parseBedParser.add_argument("-g", "--genome-file", help = "The associated genome fasta file").complete = fileCompletion

    parseBedParser.add_argument("-c", "--stratify-by-cohorts", action = "store_true", 
                                help = "Stratify results by individual donors")
    parseBedParser.add_argument("-m", "--stratify-by-Microsatellite", action = "store_true", 
                                help = "Stratify results by microsatellite stability")
    parseBedParser.add_argument("-s", "--stratify-by-Mut-Sigs", action = "store_true", 
                                help = "Stratify results by mutation signature")

    parseBedParser.add_argument("-1", "--only-sbs", action = "store_true",
                                help = "Discard all entries that are not single base substitutions")
    parseBedParser.add_argument("-i", "--include-indels", action = "store_true",
                                help = "Include insertion and deletion entries in the output file")


def formatMainPipelineParser(mainPipelineParser: ArgumentParser):

    mainPipelineParser.set_defaults(func = RunAnalysisSuite.parseArgs)
    mainPipelineParser.add_argument("--mutation-file-paths", nargs = '*',
                                    help = "One or more bed mutation files to run through the pipeline.  These files should be "
                                           "output from either parseICGC or parseBed.  If given a directory, the directory "
                                           "will be recursively searched for files ending in "
                                           "\"" + DataTypeStr.mutations + ".bed\".").complete = fileCompletion
    mainPipelineParser.add_argument("--nucleosome-maps", nargs = '*',
                                    help = "One or more nucleosome map files or directories containing them.  "
                                           "Each directory should contain files of or related to exactly ONE nucleosome map, and "
                                           "the base nucleosome map file name should match the directory name (without the file extension).")

    normGroup = mainPipelineParser.add_mutually_exclusive_group()
    normGroup.add_argument("-c", "--context-normalization", type = int, choices = (1,2,3,4,5,6),
                           help = "Choose to normalize by the surrounding context of each dyad position in a nucleosome.  "
                                  "Choice should be (1 or 2) for singlenuc/dinuc, (3 or 4) for trinuc/quadrunuc, "
                                  "or (5 or 6) for pentanuc/hexanuc.")
    normGroup.add_argument("-b", "--background", 
                           help = "The base mutation counts bed file for the data set to be used as "
                                  "the background data set for normalization.  Alternatively, "
                                  "the directory containing the data set can be supplied.").complete = fileCompletion
    normGroup.add_argument("-n", "--no-normalization", action = "store_true", 
                           help = "Choose to skip any normalization and simply produce raw nucleosome mutation counts.  "
                                  "(This is the default option)")
    mainPipelineParser.add_argument("-i", "--generate-background-immediately", action = "store_true",
                                    help = "Generate (or regenerate) nucleosome mutation counts for the background file before the "
                                           "rest of the analysis is run.  If background counts don't already exist, "
                                           "this option should definitely be selected.")

    mainPipelineParser.add_argument("-s", "--singlenuc-radius", action = "store_true",
                                    help = "Generate output files where mutations are counted within a 73 base pair radius "
                                           "of each dyad center to roughly encompass the width of a single nucleosome.")
    mainPipelineParser.add_argument("-l", "--add-linker", action = "store_true",
                                    help = "Increase the single nucleosome radius by 30 base pairs to include linker DNA "
                                           "just outside the nucleosome.")
    mainPipelineParser.add_argument("-g", "--nuc-group-radius", action = "store_true",
                                    help = "Generate output files where mutations are counted within a 1000 base pair radius "
                                           "of each dyad center to cover a group of several nucleosomes.")


def formatPeriodicityAnalysisParser(periodicityAnalysisParser: ArgumentParser):

    periodicityAnalysisParser.set_defaults(func = RunNucleosomeMutationAnalysis.parseArgs)

    periodicityAnalysisParser.add_argument("nucleosomeMutationFilePaths", nargs = '*',
                                           help = "One or more nucleosome mutation counts file paths.  These files should be "
                                                  "output from the main mutperiod pipeline.  If given a directory, the directory "
                                                  "will be recursively searched for files ending in "
                                                  "\"" + DataTypeStr.generalNucCounts + ".tsv\"").complete = fileCompletion
    periodicityAnalysisParser.add_argument("-o", "--output-file-path", 
                                           help = "The output file path for the periodicity analysis results.  "
                                                  "Should be a .rda file for complete output or a .tsv file to only "
                                                  "output the periodicity results.  The analysis process will attempt to create "
                                                  "the file if it does not already exist.").complete = fileCompletion

    periodicityAnalysisParser.add_argument("-e", "--use-expected-periodicity", action = "store_true",
                                           help = "Use the expected periodicity instead of the peak periodicity when calculating "
                                                  "statistics such as the power of the periodicity or signal-to-noise ratio.")   
    periodicityAnalysisParser.add_argument("-a", "--align-strands", action = "store_true",
                                           help = "Align the strands in the counts data so that they both run 5' to 3'.  "
                                                  "By default, the strand data is interpreted in the context of the double helix, "
                                                  "with the two strands running antiparallel to one another.")

    groupComparison = periodicityAnalysisParser.add_argument_group("Periodicity Comparison")
    groupComparison.add_argument("--group-1", nargs = '*',
                                 help = "One or more nucleosome file paths, similar to the nucleosomeMutationFilePaths argument.  "
                                        "The periodicities of this group are compared to the periodicities in the --group-2 "
                                        "argument using a Wilcoxon rank sum test.  They are also added to the regular analysis, "
                                        "so it is not necessary to include them in the nucleosomeMutationFilePaths argument as well.").complete = fileCompletion
    groupComparison.add_argument("--group-2", nargs = '*',
                                       help = "The counterpart to --group-1").complete = fileCompletion


def formatGenerateFiguresParser(generateFiguresParser: ArgumentParser):

    generateFiguresParser.set_defaults(func = GenerateFigures.parseArgs)

    generateFiguresParser.add_argument("--rda-paths", nargs = '*',
                                        help = "One or more paths to .rda files resulting from the periodicity analysis "
                                               "or directories to recursively search for those files.").complete = fileCompletion
    generateFiguresParser.add_argument("--tsv-paths", nargs = '*',
                                       help = "One or more paths to .tsv nucleosome counts files or directories to "
                                              "recursively search for those files.").complete = fileCompletion

    outputGroup = generateFiguresParser.add_mutually_exclusive_group()

    outputGroup.add_argument("-f", "--output-file",
                                       help = "A single output file to send all generated figures to.  "
                                              "Should be a single pdf file.").complete = fileCompletion

    outputGroup.add_argument("-d", "--output-directory",
                                       help = "A directory to to send each of the generated figures to.  "
                                              "Each figure is sent to a separate pdf file.").complete = fileCompletion

    generateFiguresParser.add_argument("-r", "--remove-outliers", action = "store_true",
                                       help = "Remove any statistical outliers from the data before figure generation.")

    generateFiguresParser.add_argument("-s", "--smooth-nuc-group", action = "store_true",
                                       help = "Smooth all nuc-group (translational periodicity) data in "
                                              "a 5 base pair radius to suppress rotational periodicity.")

    generateFiguresParser.add_argument("-a", "--align-strands", action = "store_true",
                                       help = "Invert counts on the minus strand so that counts for each strand "
                                              "run 5' to 3'.")
    

def getMainParser():

    # Initialize the argument parser.
    parser = argparse.ArgumentParser(description = "Run mutperiod to analyze nucleosome mutational periodicity.",
                                     prog = "mutperiod")
    subparsers = parser.add_subparsers(title='mutperiod command', required = True, dest = "mutperiod command")

    ### Create the subparsers for each relevant script.

    parseICGCParser = subparsers.add_parser("parseICGC", description = "Pass in one or more ICGC simple somatic mutation files "
                                                                       "to be parsed for use in the main pipeline.  If no additional "
                                                                       "arguments are given, the UI opens instead.")               
    formatParseICGCParser(parseICGCParser)

    
    parseBedParser = subparsers.add_parser("parseBed", description = "Pass in one or more custom bed files "
                                                                     "to be parsed for use in the main pipeline.  If no additional "
                                                                     "arguments are given, the UI opens instead.")
    formatParseBedParser(parseBedParser)                                                                  
    

    # For RunAnalysisSuite...
    mainPipelineParser = subparsers.add_parser("mainPipeline", description = "Pass in one or more mutation files formatted for mutperiod.  "
                                                                             "The mutation files are run through the primary pipeline "
                                                                             "in order to derive nucleosome mutation counts for "
                                                                             "periodicity analysis.")
    formatMainPipelineParser(mainPipelineParser)
    

    # For PeriodicityAnalysis...
    periodicityAnalysisParser = subparsers.add_parser("periodicityAnalysis", description = "Pass in one or more nucleosome mutation counts "
                                                                                           "files for periodicity analysis using a lomb-"
                                                                                           "scargle periodogram.")
    formatPeriodicityAnalysisParser(periodicityAnalysisParser)
    

    # For GenerateFigures...
    generateFiguresParser = subparsers.add_parser("generateFigures", description = "Generates figures from nucleosome counts data or "
                                                                                   "the output from the periodicity analysis.")
    formatGenerateFiguresParser(generateFiguresParser)                                                                                


    return parser


def main():

    # Run the relevant function for the subparser given.
    args = getMainParser().parse_args()
    args.func(args)


if __name__ == "__main__": main()