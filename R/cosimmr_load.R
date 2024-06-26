#' Function to load in cosimmr data and check for errors
#'
#' This function takes in the mixture data, food source means and standard
#' deviations, and (optionally) correction factor means and standard
#' deviations, and concentration proportions. It performs some (non-exhaustive)
#' checking of the data to make sure it will run through simmr. It outputs an
#' object of class \code{cosimmr_input}.
#'
#' For standard stable isotope mixture modelling, the mixture matrix will
#' contain a row for each individual and a column for each isotopic value.
#' \code{cosimmr} will allow for any number of isotopes and any number of
#' observations, within computational limits. The source means/sds should be
#' provided for each food source on each isotope. The correction means (usually
#' trophic enrichment factors) can be set as zero if required, and should be of
#' the same shape as the source values. The concentration dependence means
#' should be estimated values of the proportion of each element in the food
#' source in question and should be given in proportion format between 0 and 1.
#' At present there is no means to include concentration standard deviations.
#'
#' @param formula Formula giving in form y ~ x where y is a vector or matrix
#' of mixture values and x is a vector or matrix of covariates
#' @param source_names The names of the sources given as a character string
#' @param source_means The means of the source values, given as a matrix where
#' the number of rows is the number of sources and the number of columns is the
#' number of tracers
#' @param source_sds The standard deviations of the source values, given as a
#' matrix where the number of rows is the number of sources and the number of
#' columns is the number of tracers
#' @param correction_means The means of the correction values, given as a
#' matrix where the number of rows is the number of sources and the number of
#' columns is the number of tracers. If not provided these are set to 0.
#' @param correction_sds The standard deviations of the correction values,
#' given as a matrix where the number of rows is the number of sources and the
#' number of columns is the number of tracers. If not provided these are set to
#' 0.
#' @param concentration_means The means of the concentration values, given as a
#' matrix where the number of rows is the number of sources and the number of
#' columns is the number of tracers. These should be between 0 and 1. If not
#' provided these are all set to 1.
#' @param group A grouping variable. These can be a character or factor variable
#'
#' @import checkmate
#'
#'
#' @return An object of class \code{cosimmr_input} with the following elements:
#' \item{mixtures }{The mixture data} \item{source_names }{Source means}
#' \item{sources_sds }{Source standard deviations} \item{correction_means
#' }{Correction means} \item{correction_sds }{Correction standard deviations}
#' \item{concentration_means }{Concentration dependence means} \item{n_obs
#' }{The number of observations} \item{n_tracers }{The number of
#' tracers/isotopes} \item{n_sources }{The number of sources} \item{n_groups
#' }{The number of groups}
#' @author Emma Govan <emma.govan.2021@@mumail.ie>
#' @seealso See \code{\link{cosimmr_ffvb}} for complete examples.
#' @examples
#'
#' # A simple example with 10 observations, 2 tracers and 4 sources
#' data(geese_data_day1)
#' simmr_1 <- with(
#'   geese_data_day1,
#'   cosimmr_load(
#'     formula = mixtures ~ c(1,2,3,2,3,1,2,1,1),
#'     source_names = source_names,
#'     source_means = source_means,
#'     source_sds = source_sds,
#'     correction_means = correction_means,
#'     correction_sds = correction_sds,
#'     concentration_means = concentration_means,
#'     scale_x = TRUE
#'   )
#' )
#' 
#' Options: either input as matrix ~ whatever covariates
#' OR y1 + y2 ~ x1 + x2 + x3
#' And then use the `Formula` package
#'So then it would be 
#'mixtures = as.matrix(model.frame(Formula(formula))[,1:ncol(source_means)])
#'
#' print(simmr_1)
#' @export cosimmr_load
cosimmr_load <- function(formula,
                       source_names,
                       source_means,
                       source_sds,
                       correction_means = NULL,
                       correction_sds = NULL,
                       concentration_means = NULL,
                       scale_x = TRUE) {
  # Function to load in data for simmr and check whether it's appropriate for running through simmr_mcmc
  
  # Go through each object and check that it matches the requirements
  mixtures = as.matrix(model.frame(formula)[,1])
  if(nrow(mixtures) == 1){
    x_scaled == model.matrix(formula)
  } else{
  if(scale_x == TRUE){
    if(sd(model.matrix(formula)[,1]) == 0){
      intercept = TRUE
      scaled_mat = scale(model.matrix(formula)[,(2:ncol(model.matrix(formula)))])
  x_scaled = cbind(model.matrix(formula)[,1], 
                   scaled_mat)
  
  scaled_center = attr(scaled_mat, "scaled:center")
  
  scaled_scale = attr(scaled_mat, "scaled:scale")
    } else if(sd(model.matrix(formula)[,1]) != 0){
      intercept = FALSE
      x_scaled = scale(model.matrix(formula))
      
      scaled_center = attr(x_scaled, "scaled:center")
      
      scaled_scale = attr(x_scaled, "scaled:scale")
  }
    
  } else if(scale_x == FALSE){
    x_scaled = model.matrix(formula)
  }
  }
  # Write a function that generically tests for any 2D numeric data shape such as matrix, data frame or tibble
  assert_2D_numeric <- function(x,
                                nrows = NULL,
                                ncols = NULL,
                                null.ok = FALSE) {
    assert(
      test_data_frame(x,
                      types = c("double", "numeric"),
                      nrows = nrows,
                      ncols = ncols,
                      null.ok = null.ok
      ),
      test_matrix(x,
                  mode = "numeric",
                  nrows = nrows,
                  ncols = ncols,
                  null.ok = null.ok
      ),
      test_tibble(x,
                  types = c("double", "numeric"),
                  nrows = nrows,
                  ncols = ncols,
                  null.ok = null.ok
      )
    )
  }
  
  # Mixtures must be a matrix - the number of rows is the number of observations and the number of columns is the number of tracers
  # assert_matrix(mixtures)
  assert_2D_numeric(mixtures)
  n_obs <- nrow(mixtures)
  n_tracers <- ncol(mixtures)
  
  # Add column names if they're not there
  if (is.null(colnames(mixtures))) {
    colnames(mixtures) <- paste0("tracer", 1:n_tracers)
  }
  
  # source_names must be a character vector - the length of it is the number of sources
  assert_character(source_names)
  n_sources <- length(source_names)
  
  # source_means and source_sds must both be matrices where the number of rows is n_sources (in the same order as source_names) and the number of columns is n_tracers
  assert_2D_numeric(source_means,
                    nrows = n_sources,
                    ncols = n_tracers
  )
  # assert_matrix(source_means, nrows = n_sources, ncols = n_tracers)
  assert_2D_numeric(source_sds,
                    nrows = n_sources,
                    ncols = n_tracers
  )
  # assert_matrix(source_sds, nrows = n_sources, ncols = n_tracers)
  assert_2D_numeric(correction_means,
                    nrows = n_sources,
                    ncols = n_tracers,
                    null.ok = ifelse(is.null(correction_sds),
                                     TRUE, FALSE
                    )
  )
  # assert_matrix(correction_means,
  #   nrows = n_sources,
  #   ncols = n_tracers,
  #   null.ok = ifelse(is.null(correction_sds),
  #     TRUE, FALSE
  #   )
  # )
  assert_2D_numeric(correction_sds,
                    nrows = n_sources,
                    ncols = n_tracers,
                    null.ok = ifelse(is.null(correction_sds),
                                     TRUE, FALSE
                    )
  )
  # assert_matrix(correction_sds,
  #   nrows = n_sources,
  #   ncols = n_tracers,
  #   null.ok = ifelse(is.null(correction_means),
  #     TRUE, FALSE
  #   )
  # )
  assert_2D_numeric(concentration_means,
                    nrows = n_sources,
                    ncols = n_tracers,
                    null.ok = TRUE
  )
  # assert_matrix(concentration_means,
  #   nrows = n_sources,
  #   ncols = n_tracers, null.ok = TRUE
  # )
  
  # Fill in correction means
  if (is.null(correction_means)) {
    correction_means <- matrix(0, ncol = n_tracers, nrow = n_sources)
    correction_sds <- matrix(0, ncol = n_tracers, nrow = n_sources)
  }
  
  # concentration_means must be a matrix where all elements are less than 1
  if (is.null(concentration_means)) {
    concentration_means <- matrix(1, ncol = n_tracers, nrow = n_sources)
  } else {
    assert_true(all(concentration_means < 1) & all(concentration_means > 0))
  }
  
  # Check the groups are the right length and structure if given

  
  
  # Prepare output and give class
  out <- list(
    mixtures = mixtures,
    x_scaled = x_scaled,
    source_names = source_names,
    source_means = source_means,
    source_sds = source_sds,
    correction_means = correction_means,
    correction_sds = correction_sds,
    concentration_means = concentration_means,
    n_obs = n_obs,
    n_tracers = n_tracers,
    n_sources = n_sources,
    scale_x = scale_x,
    scaled_center = scaled_center,
    scaled_scale = scaled_scale,
    intercept = intercept
  )
  
  # Look through to see whether there are any missing values in anything
  if (any(unlist(lapply(out, "is.na")))) {
    warning("Missing values provided for some values. Check your inputs")
  }
  
  class(out) <- "cosimmr_input"
  
  return(out)
}
