#' Create partitions for a vector
#'
#' Given a vector length and the size of each partition, create a vector which
#' represents those partitions. There may be one more partition of size
#' less than the equally divided parts. For example, if the vector length is
#' 14 and the desired partitions are size 3, then there will be 4 partitions
#' of length 3 and one of length 2: 1 1 1 2 2 2 3 3 3 4 4 4 5 5.
#'
#' @param vector_length The total length of the partition vector
#' @param equal_parts The size of each partition
#'
#' @return A vector of `vector_length` divided into `equal_parts` with possibly
#'   one additional vector of size less than `equal_parts`
#'
#' @examples
#' create_partitions(14, 3)
#' create_partitions(10, 4)
#'
#' @export
create_partitions = function(vector_length, equal_parts = 100){
  c(rep(seq(1,(vector_length/equal_parts)),
        each=equal_parts),
    rep(floor(vector_length/equal_parts)+1, vector_length%%equal_parts))
}

#' Calculate stable rank response for a dataframe
#'
#' Given a dataframe with binding signals and responsiveness, calculate the
#' response ratio for each group of records partitioned by their rank.
#'
#' @importFrom dplyr arrange mutate group_by summarise
#'
#' @param df A dataframe containing binding signals and responsiveness
#' @param binding_expr_source_string A string representing the source of
#'   binding expression
#' @param bin_size The number of records per group (partition)
#' @param separator A string separator used in the binding expression
#'   source string
#'
#' @return A dataframe with response ratios calculated for each group
#'
#' @examples
#' # Given a small example dataframe similar to your input
#' df <- data.frame(binding_signal = c(1, 2, 3, 4, 5),
#'                  responsive = c(TRUE, FALSE, TRUE, TRUE, FALSE),
#'                  experiment = "test_experiment",
#'                  source_expr = "test_source_expr")
#' stable_rank_response(df, "test_experiment;test_source_expr", 2, ";")
#'
#' @export
stable_rank_response = function(df, binding_expr_source_string,
                                bin_size, separator){

  # split the binding_expr_source_string
  binding_expr_source_split = strsplit(binding_expr_source_string,
                                       separator, perl = TRUE)
  experiment = binding_expr_source_split[[1]][[1]]
  expr_src = binding_expr_source_split[[1]][[2]]

  # if the number of records is less than the bin size, reset the bin size to
  # the number of records
  bin_size = min(nrow(df), bin_size)

  df %>%
    dplyr::arrange(binding_signal) %>%
    dplyr::mutate(rank = create_partitions(nrow(.),bin_size)*bin_size) %>%
    dplyr::group_by(rank) %>%
    dplyr::summarise(group_ratio = sum(responsive)) %>%
    dplyr::mutate(response_ratio = (cumsum(group_ratio)/rank)) %>%
    dplyr::mutate(binding_src = experiment,
                  source_expr = expr_src)
}


#' Calculate rank response ratio summary
#'
#' Given a dataframe with binding signals, effect expression, and p-values,
#' calculate the rank response ratio for each group and summarize the results.
#'
#' @import dplyr
#' @importFrom tidyr pivot_wider unite
#' @importFrom purrr map2
#'
#' @param df A dataframe with binding signals, effect expression, and p-values
#' @param effect_expr_thres A threshold for effect expression (default: 0)
#' @param p_expr_thres A threshold for p-values (default: 0.05)
#' @param normalize A boolean indicating whether to normalize the results (default: FALSE)
#' @param bin_size The number of records per group (partition)
#' @param separator A string separator used in the binding expression source string
#'
#' @return A list containing two dataframes: rr with response ratios calculated
#'   for each group, and random with the random expectation values

#' @examples
#' # Given a small example dataframe similar to your input
#' df <- data.frame(
#'   experiment = rep(c("exp1", "exp2"), each = 6),
#'   source_expr = rep(c("src1", "src2"), times = 6),
#'   binding_signal = runif(6,0,.07),
#'   effect_expr = c(0.1, 0.2, -0.3, -0.4, 0.5, 0.6, 0.1,
#'                   0.2, -0.3, -0.4, 0.5, 0.6),
#'   p_expr = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.01,
#'              0.02, 0.03, 0.04, 0.05, 0.06)
#' )
#' rank_response_ratio_summarize(df, bin_size = 3, separator = ";")
#'
#' @export
rank_response_ratio_summarize = function(df,
                                         effect_expr_thres = 0,
                                         p_expr_thres = 0.05,
                                         normalize=FALSE,
                                         bin_size=5,
                                         separator=';'){

  grouped_df = df %>%
    dplyr::group_by(experiment, source_expr)

  # if normalize is set to true, find the smallest number of responsive genes
  # given the thresholds
  min_responsive =
    if (normalize == TRUE){
      grouped_df %>%
        dplyr::filter(abs(effect_expr_thres) > effect_expr_thres,
                      p_expr < p_expr_thres) %>%
        dplyr::tally() %>%
        dplyr::pull(n) %>%
        min()
    } else{
      Inf
    }

  # add a field 'responsive' to the dataframe. A gene is called responsive if
  # it passes the threshold filters, and the group rank is less than the
  # min_responsive (Inf if normalized is false, meaning all genes passing
  # thresholds are marked responsive)
  grouped_df = grouped_df %>%
    dplyr::arrange(dplyr::desc(abs(effect_expr)), .by_group = TRUE)  %>%
    dplyr::mutate(responsive =
                    ifelse(
                      abs(effect_expr) > effect_expr_thres &
                        p_expr < p_expr_thres &
                        dplyr::n() <= min_responsive,
                      TRUE,
                      FALSE))

  df_split = grouped_df %>%
    droplevels() %>%
    dplyr::group_split()

  names(df_split) = dplyr::group_keys(grouped_df) %>%
    tidyr::unite('test',c(experiment,source_expr),
                 sep=separator) %>%
    dplyr::pull()

  rr_df = purrr::map2(df_split, names(df_split),
                      stable_rank_response,
                      bin_size=bin_size,
                      separator=separator) %>%
    do.call('rbind',.)

  random_expectation_df = grouped_df  %>%
    dplyr::group_by(experiment, source_expr, responsive) %>%
    dplyr::tally() %>%
    tidyr::pivot_wider(names_from = responsive,
                       values_from = n) %>%
    # if there are on TRUE/FALSE values for a given combination of
    # expr and binding, replace that with 0
    replace(is.na(.), 0) %>%
    dplyr::mutate(random = `TRUE` / (`FALSE`+`TRUE`)) %>%
    dplyr::select(source_expr, random)

  list(
    rr = rr_df,
    random = random_expectation_df)
}
