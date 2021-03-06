% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pca.R
\name{get_top_components}
\alias{get_top_components}
\title{Get a subset of PCA components}
\usage{
get_top_components(p)
}
\arguments{
\item{p}{pca object}
}
\value{
dataset with reduced dimensionality
}
\description{
Keeps only the top n principal components using the elbow method
}
\examples{
df <- data.frame(icgc_donor_id = c('DO10324', 'DO10486', 'DO10558', 'DO10631'),
                 A1BG = rnorm(4),
                 A2ML1 = rnorm(4),
                 A4GALT = rnorm(4),
                 A4GNT = rnorm(4))
p <- do_pca(df)
get_top_components(p)
}
