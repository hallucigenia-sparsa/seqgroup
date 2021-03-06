% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compareGroups.R
\name{makeColorsGivenMetadata}
\alias{makeColorsGivenMetadata}
\title{Assign colors such that they encode a metadata item}
\usage{
makeColorsGivenMetadata(
  metadataItem = c(),
  groups = c(),
  binning = "equalfreq",
  numBins = 0,
  binLabels = c(),
  groupColorStrategy = "maxfreq",
  returnBinned = FALSE
)
}
\arguments{
\item{metadataItem}{a vector with binary, categoric or numeric values}

\item{groups}{an optional group vector with as many entries as samples in the target data set}

\item{binning}{binning strategy, only required for numeric metadata, either a string (equalfreq: the same number of samples in each group or equaldist: the same range of values) or a vector with thresholds}

\item{numBins}{the number of bins into which to bin a numeric metadata item, if zero: as many bins as samples (only required for numeric metadata and equalfreq/equaldist binning)}

\item{binLabels}{optional names for bins}

\item{groupColorStrategy}{maxfreq (assign the most frequent bin color to the group), binnum (color encodes number of bins per group)}

\item{returnBinned}{return the binned metadataItem}
}
\value{
a list with color vector and a color map, if returnBinned is true, in addition bins
}
\description{
For a binary or categoric metadata item, as many colors as categories exist are assigned.
For a numeric metadata item, the item is binned according to the equalfreq or equalwidth
strategy in as many bins as requested. Alternatively, the user can also provide thresholds to define bins.
Each bin is assigned its color. If the groups vector is empty, a color vector is returned with as many entries as samples.
If groups are given, a group color vector is returned with as many entries as groups. The color map is returned along with
the color vector.
}
\details{
Notation: a round bracket means that the end point is excluded, a square bracket that the end point is included.
For groupColorStrategy maxfreq: if there are several equally frequent bins in the group, the first one is selected.
Missing values receive color "gray".
}
