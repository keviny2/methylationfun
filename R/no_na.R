#' Checks for NA values in a vector
#'
#' If there are any NA values in the vector, the function will return FALSE. Otherwise, it will return TRUE.
#' If there are any errors, the function will return FALSE
#' @param vec vector
#'
#' @return boolean
#' @export
#'
#' @examples
#' no_na(c(1,2,3,NA,4,5)) # Returns FALSE
#' no_na(c(1,2,3,4,5,6)) # Returns TRUE
no_na <- function(vec){
  res <- tryCatch(!any(is.na(vec)), error = function(x){
    print("There was an error!")
    return(FALSE)
  })
  return(res)
}