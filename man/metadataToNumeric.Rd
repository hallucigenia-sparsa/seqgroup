% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metadataToNumeric.R
\name{metadataToNumeric}
\alias{metadataToNumeric}
\title{Convert metadata into numeric form}
\usage{
metadataToNumeric(
  metadata,
  yes = "Y",
  no = "N",
  na.threshold = 100,
  to.skip = c(),
  binarize = TRUE,
  date.items = c(),
  format = "\%d/\%m/\%y",
  referenceDate = "1/1/1900",
  remove.neg = TRUE
)
}
\arguments{
\item{metadata}{a dataframe}

\item{yes}{the symbol used for the first value in a binary metadata item (e.g. "Y")}

\item{no}{the symbol used for the second value in a binary metadata item (e.g. "N")}

\item{na.threshold}{remove metadata with more than the maximum allowed percentage of missing values}

\item{to.skip}{names of metadata to skip from binarization}

\item{binarize}{convert categoric metadata items into as many binary metadata items as there are categories (if false, metadata with more than 2 categories are removed)}

\item{date.items}{names of metadata items that represent dates}

\item{format}{the format used for the date items (an example date fitting the default format is 26/1/80)}

\item{referenceDate}{reference date used for conversion of dates into numbers (the reference date format is always d/m/Y)}

\item{remove.neg}{remove metadata items with negative values}
}
\value{
a purely numeric dataframe
}
\description{
Binary metadata items are converted into binary numeric metadata items (with 0/1).
Categoric metadata with more than two categories can be binarized, such that
each category is represented by a separate binary metadata item. Numeric metadata are kept as is.
Dates, when specified, are converted into days since the reference date.
Note that constant metadata are removed.
}
