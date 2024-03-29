import pandas
from typing import Tuple
from plotnine import ggplot, aes, geom_path, scale_color_identity, coord_cartesian, xlab, ylab, ggtitle, theme
from benbiohelpers.Plotting.PlotnineHelpers import defaultTextScaling, blankBackground


# Defines individual minor-in and minor-out positions for coloring.
# From Cui and Zhurkin, 2010 with 1 added to each side and half-base positions included.
# NOTE: These ranges are not quite suitable for coloring geom_path, since they color from the given position to the next.
#       I.e., the coloring will extend one position beyond the given range. There is also overlapping of regions which may lead to undesirable behavior.
#       However, I don't want to change the actual ranges, since they would misrepresent the actual features on the DNA.
#       I recommend that future me creates a class to wrap these data with special access methods for coloring ranges.
#       (Present me just needs to ship the analysis for the sake of my dissertation :P)
minorInPositions = set()
for start, stop in zip((4, 14, 25, 36, 46, 56, 66), (9, 19, 30, 41, 51, 61, 71)):
    minorInPositions.update(range(start-74,stop-74+1), [i - 0.5 for i in range(start+1-74,stop-74+1)])
    minorInPositions.update([abs(p) for p in minorInPositions])

minorOutPositions = set()
for start, stop in zip((9, 20, 31, 41, 51, 61), (14, 24, 35, 46, 56, 66)):
    minorOutPositions.update(range(start-74,stop-74+1), [i - 0.5 for i in range(start+1-74,stop-74+1)])
    minorOutPositions.update([abs(p) for p in minorOutPositions])

# Add the central region to minor out positions.
minorOutPositions.update([-2,-1.5,-1,-0.5,0,0.5,1,1.5,2])

# Determine transitioning positions between minor-in and minor-out regions
minorInToMinorOutPositions = set()
lastMinorIn = min(minorInPositions)
for minorIn in sorted(minorInPositions)[1:]:

    if minorIn - lastMinorIn != 0.5:

        # Determine what kind of a transition we have. Overlapping or adjacent?
        if lastMinorIn in minorOutPositions: start = lastMinorIn - 2
        else: start = lastMinorIn - 1

        # Add the relevant positions.
        for pos in range(start, lastMinorIn+2):
            minorInToMinorOutPositions.add(pos)
            minorInToMinorOutPositions.add(pos + 0.5)
        minorInToMinorOutPositions.add(lastMinorIn+2)

    lastMinorIn = minorIn

minorOutToMinorInPositions = set()
lastMinorOut = min(minorOutPositions)
for minorOut in sorted(minorOutPositions)[1:]:

    if minorOut - lastMinorOut != 0.5:

        # Determine what kind of a transition we have. Overlapping or adjacent?
        if lastMinorOut in minorInPositions: start = lastMinorOut - 2
        else: start = lastMinorOut - 1

        # Add the relevant positions.
        for pos in range(start, lastMinorOut+2):
            minorOutToMinorInPositions.add(pos)
            minorOutToMinorInPositions.add(pos + 0.5)
        minorOutToMinorInPositions.add(lastMinorOut+2)

    lastMinorOut = minorOut

# Don't forget to add the positions at the end!
minorOutToMinorInPositions.update([63,63.5,64,64.5,65,65.5,66,66.5,67])

# Set persistent column names
DYAD_POS_COL = "Dyad_Position"
PERIODIC_POS_COLOR_COL = "Periodic_Position_Color"
PERIODIC_POS_UNDERLAID_COLOR_COL = "Periodic_Position_Underlaid_Color"


class NucleosomeColorPalette:

    def __init__(self, minorIn = "#1E8449", minorOut = "#7D3C98", minorInToMinorOut = "#DC0AB4", minorOutToMinorIn = "#FFA300", intermediate = "black",
                 nucleosomal = "#B01200", linker = "#0063A6", nucleosomalUnderlaid = "#B08884", linkerUnderlaid = "#7C95A6",
                 plus = "#CC0000", minus = "#000000", plusUnderlaid = "#CC9999", minusUnderlaid = "#BFBFBF"):
        self.minorIn = minorIn
        self.minorOut = minorOut
        self.minorInToMinorOut = minorInToMinorOut
        self.minorOutToMinorIn = minorOutToMinorIn
        self.intermediate = intermediate
        self.nucleosomal = nucleosomal
        self.linker = linker
        self.nucleosomalUnderlaid = nucleosomalUnderlaid
        self.linkerUnderlaid = linkerUnderlaid
        self.plus = plus
        self.minus = minus
        self.plusUnderlaid = plusUnderlaid
        self.minusUnderlaid = minusUnderlaid


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
                                         smoothTranslational = True, nucRepLen = None, nucleosomeColorPalette = NucleosomeColorPalette(),
                                         colorRotationalTransitions = False):

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
        nucleosomeCountsData[dataCol+"_Smoothed"] = nucleosomeCountsData[dataCol].rolling(11, center = True, min_periods=1).mean()

    # Color positions
    if rotational:

        if not colorRotationalTransitions:
            # Color rotational positioning minor-in or minor-out
            nucleosomeCountsData.loc[(pos in minorInPositions for pos in nucleosomeCountsData[DYAD_POS_COL]),
                                    PERIODIC_POS_COLOR_COL] = nucleosomeColorPalette.minorIn
            nucleosomeCountsData.loc[(pos in minorOutPositions for pos in nucleosomeCountsData[DYAD_POS_COL]),
                                    PERIODIC_POS_COLOR_COL] = nucleosomeColorPalette.minorOut
            nucleosomeCountsData.loc[nucleosomeCountsData[PERIODIC_POS_COLOR_COL] == "nan", PERIODIC_POS_COLOR_COL] = nucleosomeColorPalette.intermediate
        else:
            # Color transitions between minor-in and minor-out
            nucleosomeCountsData.loc[(pos in minorInToMinorOutPositions for pos in nucleosomeCountsData[DYAD_POS_COL]),
                                    PERIODIC_POS_COLOR_COL] = nucleosomeColorPalette.minorInToMinorOut
            nucleosomeCountsData.loc[(pos in minorOutToMinorInPositions for pos in nucleosomeCountsData[DYAD_POS_COL]),
                                    PERIODIC_POS_COLOR_COL] = nucleosomeColorPalette.minorOutToMinorIn
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
        nucleosomeCountsData.loc[(pos in nucleosomePositions for pos in nucleosomeCountsData[DYAD_POS_COL].abs()),
                                 PERIODIC_POS_UNDERLAID_COLOR_COL] = nucleosomeColorPalette.nucleosomalUnderlaid
        nucleosomeCountsData.loc[(pos in linkerPositions for pos in nucleosomeCountsData[DYAD_POS_COL].abs()),
                                 PERIODIC_POS_UNDERLAID_COLOR_COL] = nucleosomeColorPalette.linkerUnderlaid

    return(nucleosomeCountsData)


# Plot the given periodicity data.
def plotPeriodicity(nucleosomeCountsData, dataCol = "Normalized_Both_Strands", title = "Periodicity", ylim = None,
                    xAxisLabel = "Position Relative to Dyad (bp)", yAxisLabel = "Repair/Damage", nucleosomeColorPalette = NucleosomeColorPalette(),
                    overlaySmoothedAndNormal = False, colorRotationalTransitions = False):

    # If overlaid smoothing is requested, determine if the relevant columns are actually present.
    if overlaySmoothedAndNormal:
        if dataCol.endswith("_Smoothed"):
            smoothedDataCol = dataCol
            dataCol = dataCol.rsplit("_Smoothed",1)[0]
        else: smoothedDataCol = dataCol + "_Smoothed"

    # Ensure that the expected data columns are actually present.
    if dataCol not in nucleosomeCountsData:
        raise ValueError(f"Expected data column named {dataCol} but it is not present in the given data frame.")
    if overlaySmoothedAndNormal and smoothedDataCol not in nucleosomeCountsData:
        raise ValueError(f"Expected smoothed data column named {smoothedDataCol} but it is not present in the given data frame.")

    periodicityPlot = ggplot(nucleosomeCountsData, aes(DYAD_POS_COL, dataCol, color = PERIODIC_POS_COLOR_COL))

    if overlaySmoothedAndNormal:
        periodicityPlot = (periodicityPlot +
            geom_path(aes(color = PERIODIC_POS_UNDERLAID_COLOR_COL, group = 1), size = 1, lineend = "round") +
            geom_path(aes(y = smoothedDataCol, group = 2), size = 1.25, lineend = "round")
        )
    else: periodicityPlot = periodicityPlot + geom_path(aes(group=1), size = 1.25, lineend = "round")

    # Determine whether the data is rotational, rotational+linker, or translational.
    rotational, rotationalPlus, translational = _getPeriodicityTypes(nucleosomeCountsData)

    if rotationalPlus:
        if not colorRotationalTransitions:
            periodicityPlot = (
                periodicityPlot +
                scale_color_identity(name = ' ', guide = "legend",
                                    breaks = (nucleosomeColorPalette.minorIn, nucleosomeColorPalette.minorOut, nucleosomeColorPalette.linker),
                                    labels = ("Minor-in", "Minor-out", "Linker"))
            )
        else:
            periodicityPlot = (
                periodicityPlot +
                scale_color_identity(name = ' ', guide = "legend",
                                    breaks = (nucleosomeColorPalette.minorInToMinorOut, nucleosomeColorPalette.minorOutToMinorIn, nucleosomeColorPalette.linker),
                                    labels = ("Minor-in to -out", "Minor-out to -in", "Linker"))
            )
    elif rotational:
        if not colorRotationalTransitions:
            periodicityPlot = (
                periodicityPlot +
                scale_color_identity(name = ' ', guide = "legend",
                                    breaks = (nucleosomeColorPalette.minorIn, nucleosomeColorPalette.minorOut),
                                    labels = ("Minor-in", "Minor-out"))
            )
        else:
            periodicityPlot = (
                periodicityPlot +
                scale_color_identity(name = ' ', guide = "legend",
                                    breaks = (nucleosomeColorPalette.minorInToMinorOut, nucleosomeColorPalette.minorOutToMinorIn),
                                    labels = ("Minor-in to -out", "Minor-out to -in"))
            )
    elif translational:
        periodicityPlot = (
            periodicityPlot +
            scale_color_identity(name = ' ', guide = "legend",
                                 breaks = (nucleosomeColorPalette.nucleosomal, nucleosomeColorPalette.linker),
                                 labels = ("Nucleosomal", "Linker"))
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
                            xAxisLabel = "Position Relative to Dyad (bp)", yAxisLabel = "Counts", overlaySmoothedAndNormal = False,
                            nucleosomeColorPalette = NucleosomeColorPalette(), colorRotationalTransitions = False):
    parsedData = parseNucleosomeCountsDataForPlotting(nucleosomeCountsData, rotationalOnlyCutoff, dataCol, smoothTranslational, nucRepLen,
                                                      nucleosomeColorPalette, colorRotationalTransitions)
    if smoothTranslational and dataCol+"_Smoothed" in parsedData:
        dataCol += "_Smoothed"
        smoothedDataPresent = True
    else: smoothedDataPresent = False
    if not smoothedDataPresent and overlaySmoothedAndNormal:
        raise ValueError("Smoothed overlay was requested, but no smoothed data was created. Is this a rotational data set?")
    return plotPeriodicity(parsedData, dataCol, title, ylim, xAxisLabel, yAxisLabel, nucleosomeColorPalette, overlaySmoothedAndNormal, colorRotationalTransitions)


# Plot the plus and minus strands (aligned) as two lines on the same graph.
def plotPlusAndMinus(nucleosomeCountsData: pandas.DataFrame, title = "Stranded Periodicity", ylim = None,
                     xAxisLabel = "Position Relative to Dyad (bp)", yAxisLabel = "Repair/Damage", nucleosomeColorPalette = NucleosomeColorPalette(),
                     smoothData = False, overlaySmoothedAndNormal = False):

    if overlaySmoothedAndNormal and not smoothData:
        raise ValueError("Cannot overlay smoothed and normal data if data is not smoothed in the first place.")

    for columnName in nucleosomeCountsData.columns:
        if "Plus_Strand" in columnName: plusStrandColumnName = columnName
        if "Minus_Strand" in columnName: minusStrandColumnName = columnName

    if smoothData:
        smoothedPlusStrandColumnName = plusStrandColumnName + "_Smoothed"
        smoothedMinusStrandColumnName = minusStrandColumnName + "_Smoothed"
        nucleosomeCountsData[smoothedPlusStrandColumnName] = nucleosomeCountsData[plusStrandColumnName].rolling(11, center = True, min_periods=1).mean()
        nucleosomeCountsData[smoothedMinusStrandColumnName] = nucleosomeCountsData[minusStrandColumnName].rolling(11, center = True, min_periods=1).mean()

    plusAndMinusPlot = ggplot(nucleosomeCountsData, aes(x = DYAD_POS_COL))

    if overlaySmoothedAndNormal:
        plusAndMinusPlot = (plusAndMinusPlot +
            geom_path(aes(y = plusStrandColumnName, color = 'nucleosomeColorPalette.plusUnderlaid'), size = 1, lineend = "round") +
            geom_path(aes(y = minusStrandColumnName, color = 'nucleosomeColorPalette.minusUnderlaid'), size = 1, lineend = "round") +
            geom_path(aes(y = smoothedPlusStrandColumnName, color = 'nucleosomeColorPalette.plus'), size = 1, lineend = "round") +
            geom_path(aes(y = smoothedMinusStrandColumnName, color = 'nucleosomeColorPalette.minus'), size = 1, lineend = "round")
        )
    elif smoothData:
        plusAndMinusPlot = (plusAndMinusPlot +
            geom_path(aes(y = smoothedPlusStrandColumnName, color = 'nucleosomeColorPalette.plus'), size = 1, lineend = "round") +
            geom_path(aes(y = smoothedMinusStrandColumnName, color = 'nucleosomeColorPalette.minus'), size = 1, lineend = "round")
        )
    else:
        plusAndMinusPlot = (plusAndMinusPlot +
            geom_path(aes(y = plusStrandColumnName, color = 'nucleosomeColorPalette.plus'), size = 1, lineend = "round") +
            geom_path(aes(y = minusStrandColumnName, color = 'nucleosomeColorPalette.minus'), size = 1, lineend = "round")
        )


    plusAndMinusPlot = (plusAndMinusPlot +
        scale_color_identity(name = ' ', guide = "legend",
                                breaks = (nucleosomeColorPalette.minus, nucleosomeColorPalette.plus),
                                labels = ("Minus Strand", "Plus Strand")) +
        coord_cartesian(ylim = ylim) + 
        xlab(xAxisLabel) + ylab(yAxisLabel) + ggtitle(title) +
        defaultTextScaling + blankBackground + theme(figure_size = (12,6))
    )

    return plusAndMinusPlot