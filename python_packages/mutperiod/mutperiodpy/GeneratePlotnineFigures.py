import pandas
from typing import Tuple
from plotnine import ggplot, aes, geom_path, scale_color_identity, coord_cartesian, xlab, ylab, ggtitle, theme
from benbiohelpers.Plotting.PlotnineHelpers import defaultTextScaling, blankBackground


# Defines individual minor-in and minor-out positions for coloring.
# From Cui and Zhurkin, 2010 with 1 added to each side and half-base positions included.
minorInPositions = set()
for start, stop in zip((4, 14, 25, 36, 46, 56, 66), (9, 19, 30, 41, 51, 61, 71)):
    minorInPositions.update(range(start-74,stop-74+1), [i - 0.5 for i in range(start+1-74,stop-74+1)])
    minorInPositions.update([abs(p) for p in minorInPositions])

minorOutPositions = set()
for start, stop in zip((9, 20, 31, 41, 51, 61), (14, 24, 35, 46, 56, 66)):
    minorOutPositions.update(range(start-74,stop-74+1), [i - 0.5 for i in range(start+1-74,stop-74+1)])
    minorOutPositions.update([abs(p) for p in minorOutPositions])

# Set persistent column names
DYAD_POS_COL = "Dyad_Position"
PERIODIC_POS_COLOR_COL = "Periodic_Position_Color"


class NucleosomeColorPalette:

    def __init__(self, minorIn = "#1E8449", minorOut = "#7D3C98", intermediate = "black",
                 nucleosomal = "#B03A2E", linker = "#2874A6", plus = "red", minus = "black"):
        self.minorIn = minorIn
        self.minorOut = minorOut
        self.intermediate = intermediate
        self.nucleosomal = nucleosomal
        self.linker = linker
        self.plus = plus
        self.minus = minus

def _getPeriodicityTypes(nucleosomeCountsData: pandas.DataFrame) -> Tuple[bool, bool, bool]:
    "Determines whether the data is rotational, rotational+linker, or translational"
    rotational = False
    rotationalPlus = False
    translational = False
    minDyadPos = min(nucleosomeCountsData[DYAD_POS_COL])
    if minDyadPos >= -73:
        rotational = True
    elif minDyadPos > -999:
        rotational = True
        rotationalPlus = True
    else: translational = True

    return (rotational, rotationalPlus, translational)


def parseNucleosomeCountsDataForPlotting(nucleosomeCountsData: pandas.DataFrame, rotationalOnlyCutoff = 60, dataCol = "Normalized_Both_Strands",
                                         smoothTranslational = True, nucRepLen = None, nucleosomeColorPalette = NucleosomeColorPalette()):

    # Determine whether the data is rotational, rotational+linker, or translational.
    rotational, rotationalPlus, translational = _getPeriodicityTypes(nucleosomeCountsData)

    # Create a copy of the original data frame so that it is not altered by this function.
    nucleosomeCountsData = nucleosomeCountsData.copy()

    # If only rotational, trim to the cutoff value
    if rotational and not rotationalPlus:
        nucleosomeCountsData = nucleosomeCountsData.loc[
            (nucleosomeCountsData[DYAD_POS_COL] >= -rotationalOnlyCutoff) & (nucleosomeCountsData[DYAD_POS_COL] <= rotationalOnlyCutoff)
        ].copy()

    # Smooth if translational and requested
    if translational and smoothTranslational:
        nucleosomeCountsData[dataCol] = nucleosomeCountsData[dataCol].rolling(11, center = True, min_periods=1).mean()

    # Color positions
    if rotational:
        # Color rotational positioning
        nucleosomeCountsData.loc[(pos in minorInPositions for pos in nucleosomeCountsData[DYAD_POS_COL]),
                                 PERIODIC_POS_COLOR_COL] = nucleosomeColorPalette.minorIn
        nucleosomeCountsData.loc[(pos in minorOutPositions for pos in nucleosomeCountsData[DYAD_POS_COL]),
                                 PERIODIC_POS_COLOR_COL] = nucleosomeColorPalette.minorOut
        nucleosomeCountsData.loc[nucleosomeCountsData[PERIODIC_POS_COLOR_COL] == "nan", PERIODIC_POS_COLOR_COL] = nucleosomeColorPalette.intermediate

    if rotationalPlus:
        # Color linker DNA in linker+ plots.
        nucleosomeCountsData.loc[(nucleosomeCountsData[DYAD_POS_COL] <= -73) | (nucleosomeCountsData[DYAD_POS_COL] >= 73), PERIODIC_POS_COLOR_COL] = nucleosomeColorPalette.linker

    if translational:

        # Derive linker and nucleosome positions from the expected period of the data.
        if nucRepLen is None: raise ValueError("Translational data found, but no nucleosome repeat length given")

        nucRepLen = round(nucRepLen)

        nucleosomePositions = set()
        nucleosomePositions.update(range(74))
        nucleosomePositions.update([i + 0.5 for i in range(0,73)])
        for i in range(1,11):
            nucleosomePositions.update(range(-73+i*nucRepLen, 73+i*nucRepLen+1))
            nucleosomePositions.update([j + 0.5 for j in range(-73+i*nucRepLen, 72+i*nucRepLen+1)])

        linkerPositions = set()
        for i in range(8):
            linkerPositions.update(range(74+i*nucRepLen,-74+(i+1)*nucRepLen+1))
            linkerPositions.update([j + 0.5 for j in range(73+i*nucRepLen,-74+(i+1)*nucRepLen+1)])

        # Color translational positioning
        nucleosomeCountsData.loc[(pos in nucleosomePositions for pos in nucleosomeCountsData[DYAD_POS_COL].abs()),
                                 PERIODIC_POS_COLOR_COL] = nucleosomeColorPalette.nucleosomal
        nucleosomeCountsData.loc[(pos in linkerPositions for pos in nucleosomeCountsData[DYAD_POS_COL].abs()),
                                 PERIODIC_POS_COLOR_COL] = nucleosomeColorPalette.linker

    return(nucleosomeCountsData)


# Plot the given periodicity data.
def plotPeriodicity(nucleosomeCountsData, dataCol = "Normalized_Both_Strands", title = "Periodicity", ylim = None,
                    xAxisLabel = "Position Relative to Dyad (bp)", yAxisLabel = "Repair/Damage", nucleosomeColorPalette = NucleosomeColorPalette()):

    periodicityPlot = (
        ggplot(nucleosomeCountsData, aes(DYAD_POS_COL, dataCol, color = PERIODIC_POS_COLOR_COL)) +
        geom_path(aes(group=1), size = 1.25, lineend = "round")
    )
    
    # Determine whether the data is rotational, rotational+linker, or translational.
    rotational, rotationalPlus, translational = _getPeriodicityTypes(nucleosomeCountsData)

    if rotationalPlus:
        periodicityPlot = (
            periodicityPlot +
            scale_color_identity(name = ' ', guide = "legend",
                                 breaks = (nucleosomeColorPalette.minorIn, nucleosomeColorPalette.minorOut, nucleosomeColorPalette.linker),
                                 labels = ("Minor-in", "Minor-out", "Linker"))
        )
    elif rotational:
        periodicityPlot = (
            periodicityPlot +
            scale_color_identity(name = ' ', guide = "legend",
                                 breaks = (nucleosomeColorPalette.minorIn, nucleosomeColorPalette.minorOut),
                                 labels = ("Minor-in", "Minor-out"))
        )
    elif translational:
        periodicityPlot = (
            periodicityPlot +
            scale_color_identity(name = ' ', guide = "legend",
                                 breaks = (nucleosomeColorPalette.nucleosomal, nucleosomeColorPalette.linker),
                                 labels = ("Nucleosomal","Linker"))
        )

    periodicityPlot = (
        periodicityPlot +
        coord_cartesian(ylim = ylim) +
        xlab(xAxisLabel) + ylab(yAxisLabel) + ggtitle(title) +
        defaultTextScaling + blankBackground + theme(figure_size = (12,6))
    )

    return periodicityPlot


# Combines the parsePeriodicityData and plotPeriodicity functions into one convenience function
def parseAndPlotPeriodicity(nucleosomeCountsData, rotationalOnlyCutoff = 60, dataCol = "Normalized_Both_Strands",
                            smoothTranslational = True, nucRepLen = None, title = "Periodicity", ylim = None,
                            xAxisLabel = "Position Relative to Dyad (bp)", yAxisLabel = "Counts"):
    parsedData = parseNucleosomeCountsDataForPlotting(nucleosomeCountsData, rotationalOnlyCutoff, dataCol, smoothTranslational, nucRepLen)
    return plotPeriodicity(parsedData, dataCol, title, ylim, xAxisLabel, yAxisLabel)


# Plot the plus and minus strands (aligned) as two lines on the same graph.
def plotPlusAndMinus(nucleosomeCountsData: pandas.DataFrame, title = "Stranded Periodicity", ylim = None,
                     xAxisLabel = "Position Relative to Dyad (bp)", yAxisLabel = "Repair/Damage", nucleosomeColorPalette = NucleosomeColorPalette()):

    for columnName in nucleosomeCountsData.columns:
        if "Plus_Strand" in columnName: plusStrandColumnName = columnName
        if "Minus_Strand" in columnName: minusStrandColumnName = columnName

    plusAndMinusPlot = (
        ggplot(nucleosomeCountsData, aes(x = DYAD_POS_COL)) +
        geom_path(aes(y = plusStrandColumnName, color = nucleosomeColorPalette.plus), size = 1.25, lineend = "round") +
        geom_path(aes(y = minusStrandColumnName, color = nucleosomeColorPalette.minus), size = 1.25, lineend = "round") +
        scale_color_identity(name = '', guide = "legend",
                                breaks = (nucleosomeColorPalette.plus, nucleosomeColorPalette.minus),
                                labels = ("Plus Strand", "Minus Strand")) +
        coord_cartesian(ylim = ylim) + 
        xlab(xAxisLabel) + ylab(yAxisLabel) + ggtitle(title) +
        defaultTextScaling + blankBackground + theme(figure_size = (12,6))
    )

    return plusAndMinusPlot