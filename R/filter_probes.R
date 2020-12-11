#' Filter probes based on mean absolute deviation
#'
#' Processes dataset and returns a subset based on most variable probes
#' @param input data frame with icgc_donor_id, probe_id and methylation_value columns, or path to file
#' @param sep delimiter for data_file, default = ','
#' @param average_duplicate average methylation value for duplicate (probe_id, methylation_value) pairs, default=FALSE
#' @param array_platform array platform to filter on
#' @param num_probes number of probes to subset
#' @param cols_to_read cols_only object indicating which columns to read if input is a file path,
#' default = cols_only(icgc_donor_id ='c', probe_id = 'c', methylation_value = 'd')
#'
#' @return long data frame subsetted on probes
#' @import dplyr
#' @import data.table
#' @export
#'
#' @examples
#' input <- data.frame(icgc_donor_id = c('DO1000','DO2000'),
#'                     probe_id = c('cg10861017','cg10861017'),
#'                     methylation_value = c(0.74, 0.92),
#'                     array_platform = rep('HumanMethylation450_after_2011_08_02', 2))
#' filter_probes(input)
filter_probes <- function(input,
                          sep=',',
                          average_duplicate=TRUE,
                          array_platform='HumanMethylation450_after_2011_08_02',
                          num_probes=5000,
                          cols_to_read=readr::cols_only(icgc_donor_id ='c',
                                                        probe_id = 'c',
                                                        methylation_value = 'd',
                                                        array_platform='c')){

  data <- load_data(input, sep, cols_to_read) %>%
    dplyr::filter(array_platform == array_platform) %>%
    dplyr::select(-c('array_platform'))

  if(average_duplicate){
    print('averaging duplicate probe_ids!')
    data <- setDT(data)[, .(methylation_value=mean(methylation_value)),
                        by='icgc_donor_id,probe_id']
  }

  data <- data %>%
    tidyr::pivot_wider(names_from = probe_id,
                       values_from = methylation_value,
                       values_fill = NA) %>%
    dplyr::select(where(no_na))


  data <- as.data.frame(data)
  rownames(data) <- data$icgc_donor_id
  data <- data %>% select(-c('icgc_donor_id'))

  # convert data to a matrix in order to use colMads
  matrix <- as.matrix(data)

  # compute mean absolute deviation of columns
  top_probes <- Rfast::colMads(matrix, method="mean")

  # decreasing order
  top_probes_arg <- order(top_probes, decreasing=TRUE)

  if(num_probes > ncol(data)){
    num_probes <- ncol(data)
  }

  subset <- data[, top_probes_arg[1:num_probes], drop=FALSE] %>%
    dplyr::mutate(icgc_donor_id = rownames(data)) %>%
    dplyr::select(icgc_donor_id, everything()) %>%
    tidyr::pivot_longer(data=.,
                        cols=names(.)[-1],
                        names_to = 'probe_id',
                        values_to = 'methylation_value')

  return(subset)
}



