% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exportToCoNet.R
\name{exportToCoNet}
\alias{exportToCoNet}
\title{Export taxon abundances in a format suitable for CoNet}
\usage{
exportToCoNet(
  abundances,
  metadata = NULL,
  lineages = NULL,
  groups = c(),
  export.folder = "",
  root.name = "conet",
  metadata.to.skip = c(),
  omitNeg = TRUE,
  yes = "Y",
  no = "N",
  date.items = c(),
  date.format = "\%d/\%m/\%y",
  taxa.are.rownames = FALSE
)
}
\arguments{
\item{abundances}{a matrix with taxa as rows and samples as columns}

\item{metadata}{an optional data frame with metadata items as columns where sample names match sample names in abundances; data types (factor vs numeric) are supposed to be correct}

\item{lineages}{an optional lineage matrix where row names match row names of taxa; taxa must be provided in the same order as in abundances; first column gives kingdom, following columns give phylum, class, order, family, genus and species (not all are required)}

\item{groups}{an optional group membership vector with as many entries as abundances has samples; if provided, data are exported separately for each group}

\item{export.folder}{the export folder into which files should be saved}

\item{root.name}{a name included in all exported files}

\item{metadata.to.skip}{names of metadata items that should not be exported}

\item{omitNeg}{if true, metadata with negative values are removed, since CoNet run with Bray Curtis or Kullback Leibler dissimilarity cannot handle them}

\item{yes}{the string encoding 1 in a binary metadata item}

\item{no}{the string encoding 0 in a binary metadata item}

\item{date.items}{names of metadata items that are dates (they will be converted into days since 1/1/1900)}

\item{date.format}{the date format (an example date fitting the default format is 26/1/80)}

\item{taxa.are.rownames}{if true, row names instead of the last column of the lineage matrix provide the lowest-level taxa (e.g. for OTUs or sequencing variants)}
}
\description{
If provided, metadata are exported into a CoNet feature file. Prior to the export, categoric
metadata are binarised, such that the resulting feature file only contains binary (as 0/1) and
numeric features, since CoNet cannot handle categoric features. Constant metadata are removed.
If binary metadata are provided as strings, the strings encoding 0 (no) and 1 (yes) have to be indicated.
If lineages are provided, they are exported into a CoNet metadata file, which can be used for all group-specific
data sets, in case groups are provided. Either the last column or the row names of the lineage matrix have to match the
row names of the abundance matrix.
}
