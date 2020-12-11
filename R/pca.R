#' Perform PCA
#'
#' Perform PCA to reduce dimension of the dataset
#' @param data a dataframe with subjects/observations on the rows and variables on the columns
#' @param transform boolean indicating whether the first column needs to be converted to rownames, default is TRUE
#'
#' @return pca object
#' @export
#'
#' @examples
#' df <- data.frame(icgc_donor_id = c('DO10324', 'DO10486', 'DO10558', 'DO10631'),
#'                  A1BG = rnorm(4),
#'                  A2ML1 = rnorm(4),
#'                  A4GALT = rnorm(4),
#'                  A4GNT = rnorm(4))
#' do_pca(df)
do_pca <- function(data, transform=TRUE){

  matrix <- data
  if(transform){
    matrix <- textshape::column_to_rownames(data)
  }

  return(PCAtools::pca(mat = matrix, removeVar = 0.1, transposed = TRUE))
}


#' Get a subset of PCA components
#'
#' Keeps only the top n principal components using the elbow method
#' @param p pca object
#'
#' @return dataset with reduced dimensionality
#' @export
#'
#' @examples
#' df <- data.frame(icgc_donor_id = c('DO10324', 'DO10486', 'DO10558', 'DO10631'),
#'                  A1BG = rnorm(4),
#'                  A2ML1 = rnorm(4),
#'                  A4GALT = rnorm(4),
#'                  A4GNT = rnorm(4))
#' p <- do_pca(df)
#' get_top_components(p)
get_top_components <- function(p){
  num_components <- PCAtools::findElbowPoint(p$variance)
  return(p$rotated[,1:num_components])
}




