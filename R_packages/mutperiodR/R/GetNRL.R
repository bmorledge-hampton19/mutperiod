#' @export
getNRL = function(nucleosomeCountsFilePath, nucleosomeExclusionBoundary = 147, plot = FALSE) {

  # This function uses lomb scargle to derive a period from "nucleosome-in-nucleosome" data,
  # taking into account where nucleosomes are not recorded due to proximity to the encompassing nucleosome.
  countsTable = data.table::fread(nucleosomeCountsFilePath)

  if (is.na(nucleosomeExclusionBoundary)) {
    counts = countsTable$Both_Strands_Counts
  } else {
    counts = countsTable[Dyad_Position > nucleosomeExclusionBoundary |
                           Dyad_Position < -nucleosomeExclusionBoundary, Both_Strands_Counts]
  }

  if (is.na(nucleosomeExclusionBoundary)) {
    times = countsTable$Dyad_Position
  } else {
    times = countsTable[Dyad_Position > nucleosomeExclusionBoundary |
                          Dyad_Position < -nucleosomeExclusionBoundary, Dyad_Position]
  }

  lombResult = lomb::lsp(counts, times, type = "period",
                         from = 50, to = 250, ofac = 100, plot = plot)

  return(lombResult$peak.at[1])

}
