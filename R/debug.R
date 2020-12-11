#' Debugger
#'
#' Used to debug dplyr code by printing out custom messages
#' @param df current dataframe in the dplyr pipeline
#' @param message custom message used for debugging purposes
#'
#' @return the same dataframe to pipe into next dplyr function
#' @export
#'
#' @examples
#' df <- data.frame(x = seq(1:10),
#'                  y = rnorm(10))
#' debug(df, "debugging!")
debug <- function(df, message){
  print(message)
  return(df)
}