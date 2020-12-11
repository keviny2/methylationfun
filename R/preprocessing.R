#' Preprocesses raw methylation data
#'
#' Process methylation data (in long form) from icgc portal through pipeline and outputs a list (methylation_matrix, donor_ids, gene_ids)
#' methylation_matrix - donors on rows, genes on columns, methyaltion values as elements of matrix
#' donor_ids - vector of donor ids for rows of methylation matrix
#' gene_ids - vector of genes for columns of methylation matrix
#' @param input data frame with icgc_donor_id, probe_id and methylation_value columns, or path to file
#' @param sep delimiter for data_file, default = ','
#' @param normalization_method method to normalize methylation values, default = 'alfers_dinges'
#' @param chunk_size size of chunks when converting probe ids to genes
#' @param average_duplicate average methylation value for duplicate (probe_id, methylation_value) pairs, default=FALSE
#' @param cols_to_read cols_only object indicating which columns to read if input is a file path,
#' default = cols_only(icgc_donor_id ='c', probe_id = 'c', methylation_value = 'd')
#'
#' @return list object: data.frame matrix whose rows are donors and columns are genes, gene dataframe, donor id dataframe
#' @import dplyr
#' @import data.table
#' @export
#'
#' @examples
#' data <- data.frame(icgc_donor_id = c('DO1000'),
#'                    probe_id = 'cg10861017',
#'                    methylation_value = 0.74)
#' process_raw_methylation_data(data, 'alfers_dinges')
process_raw_methylation_data <- function(input,
                                         sep=',',
                                         normalization_method='alfers_dinges',
                                         chunk_size=10000,
                                         average_duplicate=FALSE,
                                         cols_to_read=cols_only(icgc_donor_id ='c',
                                                                probe_id = 'c',
                                                                methylation_value = 'd')) {

  library(FDb.InfiniumMethylation.hg19)

  data <- load_data(input, sep, cols_to_read)

  if(average_duplicate){
    print('averaging duplicate probe_ids!')
    data <- setDT(data)[, .(methylation_value=mean(methylation_value)),
                        by='icgc_donor_id,probe_id']
  }

  # get normalization method
  normal_approx <- get_normalization_method(normalization_method)

  # load a GRanges object for use later in order to speed up runtime
  hm450 <- get450k()


  methylation_matrix <- data %>%

    debug('first pivot_wider!') %>%

    # going from long to wide dataframe
    tidyr::pivot_wider(names_from = probe_id, values_from = methylation_value, values_fill = NA) %>%

    debug('keeping probe_ids shared across all donors!') %>%

    select(where(no_na)) %>%

    debug('first pivot_longer!') %>%

    # revert back to long dataframe
    tidyr::pivot_longer(data=.,
                        cols=names(.)[-1],
                        names_to = 'probe_id',
                        values_to = 'methylation_value') %>%

    debug('converting probes to genes!') %>%

    group_by(icgc_donor_id) %>%
    do(data.frame(
      icgc_donor_id = .$icgc_donor_id,
      nearest_gene_symbol = convert_probe_to_gene(.$probe_id, unique(.$icgc_donor_id), hm450, chunk_size),
      methylation_value_norm_approx = normal_approx(.$methylation_value)
    )) %>%
    ungroup() %>%


    debug('filtering ERRORS!') %>%
    # remove rows that have 'ERROR'
    # (since donors are processed by groups, it should remove all data for the patient if there was an error)
    filter(nearest_gene_symbol != 'ERROR') %>%

    debug('averaging duplicate genes!') %>%

    # average duplicate genes (perhaps two probes map to the same gene)
    group_by(icgc_donor_id, nearest_gene_symbol) %>%
    summarize(methylation_value_norm_approx = mean(methylation_value_norm_approx)) %>%
    ungroup() %>%

    debug('second pivot_wider!') %>%

    tidyr::pivot_wider(names_from = nearest_gene_symbol, values_from = methylation_value_norm_approx, values_fill = NA) %>%

    debug('only keep genes shared across all donors!') %>%

    # only keep genes that are shared across all donors
    select(where(no_na))

  donor_ids <- methylation_matrix %>% select(c('icgc_donor_id'))
  methylation_matrix <- methylation_matrix %>% select(-c('icgc_donor_id'))
  gene_ids <- as.data.frame(colnames(methylation_matrix))
  # tissue_locations <- get_tissue_locations()  TODO: add this to the function

  return(list(methylation_matrix=methylation_matrix,
              donor_ids=donor_ids,
              gene_ids=gene_ids))
}


#' Bind multiple (wide) data frames together
#'
#' @param ... list of data frames
#' @param remove_na remove columns with NA values? defualt = TRUE
#'
#' @return binded dataset
#' @import dplyr
#' @export
#'
#' @examples
#' df1 <- data.frame(icgc_donor_id = 'DO1111',
#'                   A1BG = 0.4,
#'                   AAAS = 0.5)
#' df2 <- data.frame(icgc_donor_id = 'DO2222',
#'                   A1BG = 0.4,
#'                   AAAS = 0.5)
#' bind_data_frames(df1,df2)
bind_data_frames <- function(..., remove_na = TRUE){
  binded <- bind_rows(...)

  if(remove_na){
    binded <- dplyr::select(binded, where(no_na))

  }

  return(binded)
}


#' User friendly data loader
#'
#' @param input can be a file name or a data frame
#' @param sep delimiter if input is a file name
#' @param cols_to_read cols_only object indicating which columns to read if input is a file path, default = NULL
#'
#' @return loaded data
#' @export
#'
#' @examples
#' load_data(data.frame(x = c(1,2,3),
#'                      y = c(4,5,6)))
load_data <- function(input, sep, cols_to_read=NULL){
  if(is.character(input) & length(input) == 1){
    data <- readr::read_delim(input,
                              sep,
                              col_types=cols_to_read)
  } else{
    data <- input
  }

  return(data)
}


#' Map donor ids to tissue locations using project code
#'
#' @param input donor dataset that maps icgc donor ids to project codes or file name
#' (only requires icgc_donor_id and project_code columns)
#' @param sep character, delimiter if input is file name
#' @param icgc_donor_ids dataframe with the relevant donor ids
#'
#' @return dataframe with the mapped tissues
#' @export
#'
#' @examples
#' input <- data.frame(icgc_donor_id = c('DO1111', 'DO2222'),
#'                     project_code = c('BRCA-US', 'OV-AU'))
#' icgc_donor_ids <- data.frame(icgc_donor_id = c('DO1111', 'DO2222'))
#' get_tissue_locations(input=input, icgc_donor_ids=icgc_donor_ids)
get_tissue_locations <- function(input, sep=NA, icgc_donor_ids){
  donor <- load_data(input, sep, readr::cols_only(icgc_donor_id='c', project_code='c'))
  donor <- donor[match(icgc_donor_ids$icgc_donor_id, donor$icgc_donor_id),]
  mapped <- dplyr::inner_join(donor, dict, by = c('project_code'))
  mapped <- mapped[!duplicated(mapped$icgc_donor_id),]

  return(data.frame(tissue = mapped$tissue))
}

#' Entire preprocessing pipeline. First filters probes then does probe to gene conversion
#'
#' @param input data frame with icgc_donor_id, probe_id and methylation_value columns, or path to file
#' @param filename filename for writing matrix to .csv file
#' @param sep delimiter for input, default = '\t'
#' @param average_duplicate average methylation value for duplicate (probe_id, methylation_value) pairs, default=FALSE
#' @param array_platform array platform to filter on
#' @param num_probes number of probes to subset
#' @param normalization_method method to normalize methylation values, default = 'alfers_dinges'
#' @param chunk_size size of chunks when converting probe ids to genes
#' @param cols_to_read_probe cols_only object indicating which columns to read if input is a file path,
#' default = cols_only(icgc_donor_id ='c', probe_id = 'c', methylation_value = 'd')
#' @param cols_to_read_meth cols_only object indicating which columns to read if input is a file path,
#' default = cols_only(icgc_donor_id ='c', probe_id = 'c', methylation_value = 'd')
#' @param tissue_input data frame with icgc_donor_id, probe_id and methylation_value columns, or path to file
#' @param tissue_sep delimiter if tissue_input is a filename
#'
#' @return
#' @export
#'
#' @examples
#' input <- data.frame(icgc_donor_id = c('DO1000','DO2000'),
#'                     probe_id = c('cg10861017','cg10861017'),
#'                     methylation_value = c(0.74, 0.92),
#'                     array_platform = rep('HumanMethylation450_after_2011_08_02', 2))
#'
#' tissue_input <- data.frame(icgc_donor_id = c('DO1111', 'DO2222'),
#'                            project_code = c('BRCA-US', 'OV-AU'))
#'
#' methylation_pipeline(input=input,
#'                      tissue_input=tissue_input)
methylation_pipeline <- function(input,
                                 filename = NA,
                                 sep='\t',
                                 average_duplicate=TRUE,
                                 array_platform='HumanMethylation450_after_2011_08_02',
                                 num_probes=5000,
                                 normalization_method='alfers_dinges',
                                 chunk_size=10000,
                                 cols_to_read_probe=readr::cols_only(icgc_donor_id ='c',
                                                                     probe_id = 'c',
                                                                     methylation_value = 'd',
                                                                     array_platform ='c'),
                                 cols_to_read_meth=readr::cols_only(icgc_donor_id ='c',
                                                                    probe_id = 'c',
                                                                    methylation_value = 'd'),
                                 tissue_input='donor.tsv',
                                 tissue_sep='\t'){

  res <- filter_probes(input = input,
                       sep = sep,
                       average_duplicate = average_duplicate,
                       array_platform = array_platform,
                       num_probes = num_probes,
                       cols_to_read = cols_to_read_probe) %>%
           process_raw_methylation_data(chunk_size = num_probes,
                                 cols_to_read = cols_to_read_meth)

  res[['tissue_locations']] <- get_tissue_locations(tissue_input, tissue_sep, res$donor_ids)

  if(!is.na(filename)) {
    readr::write_csv(res$methylation_matrix, filename)
    readr::write_csv(res$donor_ids, 'icgc_donor_id.csv')
    readr::write_csv(res$gene_ids, 'gene_id.csv')
    readr::write_csv(res$tissue_locations, 'tissue_location_data.csv')
  }

  return(res)
}



