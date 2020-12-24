#' Generate long format hazards data for pooled hazards estimation
#'
#' @details Generates a long-form dataset that represents each observation in
#'  terms of repeated measures across discretized bins derived from selecting
#'  break points over the support of A. This repeated measures dataset is
#'  suitable for estimating the hazard of failing in a particular bin over A
#'  using a highly adaptive lasso classification model.
#'
#' @param A The \code{numeric} vector or similar of the observed values of an
#'  intervention for a group of observational units of interest.
#' @param W A \code{data.frame}, \code{matrix}, or similar giving the values of
#'  baseline covariates (potential confounders) for the observed units whose
#'  observed intervention values are provided in the previous argument.
#' @param wts A \code{numeric} vector of observation-level weights. The default
#'  is to weight all observations equally.
#' @param grid_type A \code{character} indicating the strategy (or strategies)
#'  to be used in creating bins along the observed support of the intervention
#'  \code{A}. For bins of equal range, use "equal_range"; consult documentation
#'  of \code{\link[ggplot2]{cut_interval}} for more information. To ensure each
#'  bin has the same number of points, use "equal_mass"; consult documentation
#'  of \code{\link[ggplot2]{cut_number}} for details.
#' @param n_bins Only used if \code{grid_type} is set to \code{"equal_range"}
#'  or \code{"equal_mass"}. This \code{numeric} value indicates the number(s)
#'  of bins into which the support of \code{A} is to be divided.
#' @param breaks A \code{numeric} vector of break points to be used in dividing
#'  up the support of \code{A}. This is passed through the \code{...} argument
#'  to \code{\link[base]{cut.default}} by \code{\link[ggplot2]{cut_interval}}
#'  or \code{\link[ggplot2]{cut_number}}.
#'
#' @importFrom data.table as.data.table setnames
#' @importFrom ggplot2 cut_interval cut_number
#' @importFrom future.apply future_lapply
#' @importFrom assertthat assert_that
#'
#' @return A \code{list} containing the break points used in dividing the
#'  support of \code{A} into discrete bins, the length of each bin, and the
#'  reformatted data. The reformatted data is a \code{\link{data.table}} of
#'  repeated measures data, with an indicator for which bin an observation
#'  fails in, the bin ID, observation ID, values of \code{W} for each given
#'  observation, and observation-level weights.
format_long_hazards <- function(A, W, wts = rep(1, length(A)),
                                grid_type = c(
                                  "equal_range", "equal_mass"
                                ),
                                n_bins = NULL, breaks = NULL) {
  # clean up arguments
  grid_type <- match.arg(grid_type)

  # set grid along A and find interval membership of observations along grid
  if (is.null(breaks) & !is.null(n_bins)) {
    if (grid_type == "equal_range") {
      bins <- ggplot2::cut_interval(A, n_bins,
                                    right = FALSE,
                                    ordered_result = TRUE, dig.lab = 12
      )
    } else if (grid_type == "equal_mass") {
      bins <- ggplot2::cut_number(A, n_bins,
                                  right = FALSE,
                                  ordered_result = TRUE, dig.lab = 12
      )
    }
    # https://stackoverflow.com/questions/36581075/extract-the-breakpoints-from-cut
    breaks_left <- as.numeric(sub(".(.+),.+", "\\1", levels(bins)))
    breaks_right <- as.numeric(sub(".+,(.+).", "\\1", levels(bins)))
    bin_length <- round(breaks_right - breaks_left, 3)
    bin_id <- as.numeric(bins)
    all_bins <- matrix(seq_len(max(bin_id)), ncol = 1)
    # for predict method, only need to assign observations to existing intervals
  } else if (!is.null(breaks)) {
    # NOTE: findInterval() and cut() might return slightly different results...
    bin_id <- findInterval(A, breaks, all.inside = TRUE, rightmost.closed = T)
    all_bins <- matrix(seq_along(breaks), ncol = 1)
  } else {
    stop("Combination of arguments `breaks`, `n_bins` incorrectly specified.")
  }

  # loop over observations to create expanded set of records for each
  reformat_each_obs <- future.apply::future_lapply(seq_along(A), function(i) {
    # create indicator and "turn on" indicator for interval membership
    bin_indicator <- rep(0, nrow(all_bins))
    bin_indicator[bin_id[i]] <- 1
    id <- rep(i, nrow(all_bins))

    # get correct value of baseline variables and repeat along intervals
    if (is.null(dim(W))) {
      # assume vector
      obs_w <- rep(W[i], nrow(all_bins))
      names_w <- "W"
    } else {
      # assume two-dimensional array
      obs_w <- rep(as.numeric(W[i, ]), nrow(all_bins))
      obs_w <- matrix(obs_w, ncol = ncol(W), byrow = TRUE)

      # use names from array if present
      if (is.null(names(W))) {
        names_w <- paste("W", seq_len(ncol(W)), sep = "_")
      } else {
        names_w <- names(W)
      }
    }

    # get correct value of weights and repeat along intervals
    # NOTE: the weights are always a vector
    obs_wts <- rep(wts[i], nrow(all_bins))

    # create data table with membership indicator and interval limits
    suppressWarnings(
      hazards_df <- data.table::as.data.table(cbind(
        id, bin_indicator,
        all_bins, obs_w,
        obs_wts
      ))
    )

    # trim records to simply end at the failure time for a given observation
    hazards_df_reduced <- hazards_df[seq_len(bin_id[i]), ]

    # give explicit names and add to appropriate position in list
    hazards_df <-
      data.table::setnames(
        hazards_df_reduced,
        c("obs_id", "in_bin", "bin_id", names_w, "wts")
      )
    return(hazards_df)
  })

  # combine observation-level hazards data into larger structure
  reformatted_data <- do.call(rbind, reformat_each_obs)
  out <- list(
    data = reformatted_data,
    breaks =
      if (exists("breaks_left")) {
        c(breaks_left, last(breaks_right))
      } else {
        NULL
      },
    bin_length =
      if (exists("bin_length")) {
        bin_length
      } else {
        NULL
      }
  )
  return(out)
}

###############################################################################

#' Map a predicted hazard to a predicted density for a single observation
#'
#' @details For a single observation, map a predicted hazard of failure (as an
#'  occurrence in a particular bin, under a given partitioning of the support)
#'  to a density.
#'
#' @param hazard_pred_single_obs A \code{numeric} vector of predicted hazard of
#'  failure in a given bin (under a given partitioning of the support) for a
#'  single observational unit based on a long format data structure (from
#'  \code{\link{format_long_hazards}}). This is the probability that a given
#'  value falls in a corresponding bin, given that it has not yet failed
#'  (fallen in a preceding bin), as per \insertRef{diaz2011super}{haldensify}.
#'
#' @importFrom assertthat assert_that
#'
#' @return A \code{matrix} composed of a single row and a number of columns
#'  specified by the grid of penalization parameters used in fitting of the
#'  highly adaptive lasso. This is the predicted conditional density for a
#'  given observation, re-mapped from the hazard scale.
map_hazard_to_density <- function(hazard_pred_single_obs) {
  # number of records for the given observation
  n_records <- nrow(hazard_pred_single_obs)

  # NOTE: pred_hazard = (1 - pred) if 0 in this bin * pred if 1 in this bin
  if (n_records > 1) {
    hazard_prefailure <- matrix(1 - hazard_pred_single_obs[-n_records, ],
                                nrow = (n_records - 1)
    )
    hazard_at_failure <- hazard_pred_single_obs[n_records, ]
    hazard_predicted <- rbind(hazard_prefailure, hazard_at_failure)
    rownames(hazard_predicted) <- NULL
  } else {
    hazard_predicted <- hazard_pred_single_obs
  }

  # sanity check of dimensions
  assertthat::assert_that(all(dim(hazard_pred_single_obs) ==
                                dim(hazard_predicted)))

  # multiply hazards across rows to construct the individual-level density
  density_pred_from_hazards <- matrix(apply(hazard_predicted, 2, prod),
                                      nrow = 1
  )
  return(density_pred_from_hazards)
}









utils::globalVariables(c("wts"))

#' Prediction method for HAL-based conditional density estimation
#'
#' @details Method for computing and extracting predictions of the conditional
#'  density estimates based on the highly adaptive lasso estimator, returned as
#'  an S3 object of class \code{haldensify} from \code{\link{haldensify}}.
#'
#' @param object An object of class \code{\link{haldensify}}, containing the
#'  results of fitting the highly adaptive lasso for conditional density
#'  estimation, as produced by a call to \code{\link{haldensify}}.
#' @param ... Additional arguments passed to \code{predict} as necessary.
#' @param new_A The \code{numeric} vector or similar of the observed values for
#'  which a conditional density estimate is to be generated.
#' @param new_W A \code{data.frame}, \code{matrix}, or similar giving the
#'  values of baseline covariates (potential confounders) for the conditioning
#'  set of the observed values \code{A}.
#' @param cv_select A \code{logical} indicating whether to return the predicted
#'  density for the value of the regularization parameter selected by global
#'  cross-validation. The default is \code{TRUE}. When set to \code{FALSE}, a
#'  matrix of predicted densities is returned, with each column corresponding
#'  to a value of the regularization parameter less than or equal to the choice
#'  made by the global cross-validation selector.
#'
#' @importFrom stats predict
#' @importFrom data.table ":="
#'
#' @return A \code{numeric} vector of predicted conditional density values from
#'  a fitted \code{haldensify} object.
#'
#' @export
#'
#' @examples
#' # simulate data: W ~ U[-4, 4] and A|W ~ N(mu = W, sd = 0.5)
#' n_train <- 50
#' w <- runif(n_train, -4, 4)
#' a <- rnorm(n_train, w, 0.5)
#' # HAL-based density estimator of A|W
#' mod_haldensify <- haldensify(
#'   A = a, W = w, n_bins = 3,
#'   lambda_seq = exp(seq(-1, -10, length = 50))
#' )
#' # predictions to recover conditional density of A|W
#' new_a <- seq(-4, 4, by = 0.1)
#' new_w <- rep(0, length(new_a))
#' pred_dens <- predict(mod_haldensify, new_A = new_a, new_W = new_w)
predict.haldensify <- function(object, ..., new_A, new_W,
                               cv_select = TRUE) {
  # make long format data structure with new input data
  long_format_args <- list(
    A = new_A,
    W = new_W,
    grid_type = object$call$grid_type,
    breaks = object$breaks
  )
  reformatted_output <- do.call(format_long_hazards, long_format_args)
  long_data_pred <- reformatted_output$data
  long_data_pred[, wts := NULL]

  # predict conditional density estimate from HAL fit on new long format data
  # over the sequence of lambda less than or equal to CV-selected lambda
  hazard_pred <- stats::predict(
    object = object$hal_fit,
    new_data =
      long_data_pred[, 3:ncol(long_data_pred)]
  )

  # NOTE: we return hazard predictions for the loss minimizer and all lambda
  #       smaller than it, BUT if there are no such lambda, hazard_pred is only
  #       vector rather than the usually expected matrix
  if (!is.matrix(hazard_pred)) hazard_pred <- as.matrix(hazard_pred, ncol = 1)

  # estimate unscaled density for each observation and each lambda
  density_pred_rescaled <- apply(hazard_pred, 2, function(this_hazard_pred) {
    # coerce single column of predictions back to matrix
    this_hazard_pred <- matrix(this_hazard_pred, ncol = 1)

    # compute hazard for a given observation by looping over individuals
    dens_given_lambda <- lapply(unique(long_data_pred$obs_id), function(id) {
      # get predictions for the current observation only
      hazard_pred_this_obs <-
        matrix(this_hazard_pred[long_data_pred$obs_id == id, ])

      # map hazard to density for a single observation and return
      density_pred_this_obs <-
        map_hazard_to_density(hazard_pred_single_obs = hazard_pred_this_obs)

      # return density for a single observation
      return(as.numeric(density_pred_this_obs))
    })
    # aggregate predicted unscaled density at the level of observations
    density_pred_unscaled <- do.call(c, dens_given_lambda)

    # re-scale to densities by dividing by bin widths
    density_pred_scaled <- density_pred_unscaled /
      object$bin_sizes[long_data_pred[in_bin == 1, bin_id]]

    # return re-scaled densities
    return(density_pred_scaled)
  })

  # truncate predictions outside range of observed A
  outside_range <- new_A < object$range_a[1] | new_A > object$range_a[2]
  density_pred_rescaled[outside_range, ] <- 0

  # return predicted densities only for CV-selected lambda
  # NOTE: turn off for access to density estimates for all lambda >= CV-lambda
  if (cv_select) {
    density_pred_rescaled <- density_pred_rescaled[, 1]
  }

  # output
  return(density_pred_rescaled)
}








utils::globalVariables(c("in_bin", "bin_id"))

#' Conditional density estimation with HAL in a single cross-validation fold
#'
#' @details Estimates the conditional density of A|W for a subset of the full
#'  set of observations based on the inputted structure of the cross-validation
#'  folds. This is a helper function intended to be used to select the optimal
#'  value of the penalization parameter for the highly adaptive lasso estimates
#'  of the conditional hazard (via \code{\link[origami]{cross_validate}}). The
#'
#' @param fold Object specifying cross-validation folds as generated by a call
#'  to \code{\link[origami]{make_folds}}.
#' @param long_data A \code{data.table} or \code{data.frame} object containing
#'  the data in long format, as given in \insertRef{diaz2011super}{haldensify},
#'  as produced by \code{\link{format_long_hazards}}.
#' @param wts A \code{numeric} vector of observation-level weights, matching in
#'  its length the number of records present in the long format data. Default
#'  is to weight all observations equally.
#' @param lambda_seq A \code{numeric} sequence of values of the tuning
#'  parameter of the Lasso L1 regression passed to
#'  \code{\link[hal9001]{fit_hal}}.
#'
#' @importFrom stats aggregate plogis
#' @importFrom origami training validation fold_index
#' @importFrom assertthat assert_that
#' @importFrom hal9001 fit_hal
#' @importFrom Rdpack reprompt
#'
#' @return A \code{list}, containing density predictions, observations IDs,
#'  observation-level weights, and cross-validation indices for conditional
#'  density estimation on a single fold of the overall data.
cv_haldensify <- function(fold, long_data, wts = rep(1, nrow(long_data)),
                          lambda_seq = exp(seq(-1, -13, length = 100))) {
  # make training and validation folds
  train_set <- origami::training(long_data)
  valid_set <- origami::validation(long_data)

  # subset observation-level weights to the correct size
  wts_train <- wts[fold$training_set]
  wts_valid <- wts[fold$validation_set]

  # fit a HAL regression on the training set
  # NOTE: not selecting lambda by CV so no need to pass IDs for fold splitting
  hal_fit_train <- hal9001::fit_hal(
    X = as.matrix(train_set[, -c(1, 2)]),
    Y = as.numeric(train_set$in_bin),
    max_degree = NULL,
    fit_type = "glmnet",
    family = "binomial",
    lambda = lambda_seq,
    cv_select = FALSE,
    standardize = FALSE, # pass to glmnet
    weights = wts_train, # pass to glmnet
    yolo = FALSE
  )

  # get intercept and coefficient fits for this value of lambda from glmnet
  alpha_hat <- hal_fit_train$glmnet_lasso$a0
  betas_hat <- hal_fit_train$glmnet_lasso$beta
  coefs_hat <- rbind(alpha_hat, betas_hat)

  # make design matrix for validation set manually
  pred_x_basis <- hal9001::make_design_matrix(
    as.matrix(valid_set[, -c(1, 2)]),
    hal_fit_train$basis_list
  )
  pred_x_basis <- hal9001::apply_copy_map(
    pred_x_basis,
    hal_fit_train$copy_map
  )
  pred_x_basis <- cbind(rep(1, nrow(valid_set)), pred_x_basis)

  # manually predict along sequence of lambdas
  preds_logit <- pred_x_basis %*% coefs_hat
  preds <- stats::plogis(as.matrix(preds_logit))

  # compute hazard for a given observation by looping over individuals
  density_pred_each_obs <- lapply(unique(valid_set$obs_id), function(id) {
    # get predictions for the current observation only
    hazard_pred_this_obs <- matrix(preds[valid_set$obs_id == id, ],
                                   ncol = length(lambda_seq)
    )

    # map hazard to density for a single observation and return
    density_pred_this_obs <-
      map_hazard_to_density(hazard_pred_single_obs = hazard_pred_this_obs)

    return(density_pred_this_obs)
  })

  # aggregate predictions across observations
  density_pred <- do.call(rbind, as.list(density_pred_each_obs))

  # collapse weights to the observation level
  wts_valid_reduced <- stats::aggregate(
    wts_valid, list(valid_set$obs_id),
    unique
  )
  colnames(wts_valid_reduced) <- c("id", "weight")

  # construct output
  out <- list(
    preds = density_pred,
    ids = wts_valid_reduced$id,
    wts = wts_valid_reduced$weight,
    fold = origami::fold_index()
  )
  return(out)
}

###############################################################################

#' Cross-validated conditional density estimation with HAL
#'
#' @details Estimation of the conditional density A|W through using the highly
#'  adaptive lasso to estimate the conditional hazard of failure in a given
#'  bin over the support of A. Cross-validation is used to select the optimal
#'  value of the penalization parameters, based on minimization of the weighted
#'  log-likelihood loss for a density.
#'
#' @param A The \code{numeric} vector observed values.
#' @param W A \code{data.frame}, \code{matrix}, or similar giving the values of
#'  baseline covariates (potential confounders) for the observed units. These
#'  make up the conditioning set for the conditional density estimate.
#' @param wts A \code{numeric} vector of observation-level weights. The default
#'  is to weight all observations equally.
#' @param grid_type A \code{character} indicating the strategy to be used in
#'  creating bins along the observed support of \code{A}. For bins of equal
#'  range, use \code{"equal_range"}; consult the documentation of
#'  \code{\link[ggplot2]{cut_interval}} for more information. To ensure each
#'  bin has the same number of observations, use \code{"equal_mass"}; consult
#'  the documentation of \code{\link[ggplot2]{cut_number}} for details. The
#'  default is \code{"equal_range"} since this has been found to provide better
#'  performance in simulation experiments; however, both types may be specified
#'  (i.e., \code{c("equal_range", "equal_mass")}) together, in which case
#'  cross-validation will be used to select the optimal binning strategy.
#' @param n_bins This \code{numeric} value indicates the number(s) of bins into
#'  which the support of \code{A} is to be divided. As with \code{grid_type},
#'  multiple values may be specified, in which case a cross-validation selector
#'  will be used to choose the optimal number of bins. In fact, the default
#'  uses a cross-validation selector to choose between 10 and 25 bins.
#' @param cv_folds A \code{numeric} indicating the number of cross-validation
#'  folds to be used in fitting the sequence of HAL conditional density models.
#' @param lambda_seq A \code{numeric} sequence of values of the regularization
#'  parameter of Lasso regression; passed to \code{\link[hal9001]{fit_hal}}.
#' @param use_future A \code{logical} indicating whether to attempt to use
#'  parallelization based on \pkg{future} and \pkg{future.apply}. If set to
#'  \code{TRUE}, \code{\link[future.apply]{future_mapply}} will be used in
#'  place of \code{mapply}. When set to \code{TRUE}, a parallelization scheme
#'  must be specified externally by a call to \code{\link[future]{plan}}.
#'
#' @importFrom data.table ":="
#' @importFrom future.apply future_mapply
#' @importFrom hal9001 fit_hal
#'
#' @return Object of class \code{haldensify}, containing a fitted
#'  \code{hal9001} object; a vector of break points used in binning \code{A}
#'  over its support \code{W}; sizes of the bins used in each fit; the tuning
#'  parameters selected by cross-validation; the full sequence (in lambda) of
#'  HAL models for the CV-selected number of bins and binning strategy; and
#'  the range of \code{A}.
#'
#' @export
#'
#' @examples
#' # simulate data: W ~ U[-4, 4] and A|W ~ N(mu = W, sd = 0.5)
#' n_train <- 50
#' w <- runif(n_train, -4, 4)
#' a <- rnorm(n_train, w, 0.5)
#' # learn relationship A|W using HAL-based density estimation procedure
#' mod_haldensify <- haldensify(
#'   A = a, W = w, n_bins = 3,
#'   lambda_seq = exp(seq(-1, -10, length = 50))
#' )
haldensify <- function(A,
                       W,
                       wts = rep(1, length(A)),
                       grid_type = "equal_range",
                       n_bins = c(10, 25),
                       cv_folds = 5,
                       lambda_seq = exp(seq(-1, -13, length = 1000)),
                       use_future = FALSE) {
  # run CV-HAL for all combinations of n_bins and grid_type
  tune_grid <- expand.grid(
    grid_type = grid_type, n_bins = n_bins,
    stringsAsFactors = FALSE
  )

  # run procedure to select tuning parameters via cross-validation
  if (use_future) {
    select_out <- future.apply::future_mapply(
      FUN = fit_haldensify,
      grid_type = tune_grid$grid_type,
      n_bins = tune_grid$n_bins,
      MoreArgs = list(
        A = A, W = W, wts = wts,
        cv_folds = cv_folds,
        lambda_seq = lambda_seq
      ),
      SIMPLIFY = FALSE,
      future.seed = TRUE
    )
  } else {
    select_out <- mapply(
      FUN = fit_haldensify,
      grid_type = tune_grid$grid_type,
      n_bins = tune_grid$n_bins,
      MoreArgs = list(
        A = A, W = W, wts = wts,
        cv_folds = cv_folds,
        lambda_seq = lambda_seq
      ),
      SIMPLIFY = FALSE
    )
  }

  # extract n_bins/grid_type index that is empirical loss minimizer
  emp_risk_per_lambda <- lapply(select_out, `[[`, "emp_risks")
  min_loss_idx <- lapply(emp_risk_per_lambda, which.min)
  min_risk <- lapply(emp_risk_per_lambda, min)
  tune_select_params <- tune_grid[which.min(min_risk), , drop = FALSE]
  tune_select_fits <- select_out[[which.min(min_risk)]]


  # get index of CV-selected lambda; subset sequence to that + smaller lambdas
  lambda_selected_idx <- tune_select_fits$lambda_loss_min_idx
  # lambda_selected_idx <- 1
  temp_cut <- lambda_seq[lambda_selected_idx]
  # temp_cut <- max(lambda_seq[lambda_selected_idx], 0.01)
  # if(temp_cut < 0.01) {
  #   lambda_selected_idx <- max(lambda_selected_idx - 20, 1)
  #   temp_cut <- lambda_seq[lambda_selected_idx]
  # }
  lambda_seq_usm <- lambda_seq[lambda_seq <= temp_cut]
  lambda_seq_usm <- NULL


  # re-format input data into long hazards structure
  reformatted_output <- format_long_hazards(
    A = A, W = W, wts = wts,
    grid_type = tune_select_params$grid_type,
    n_bins = tune_select_params$n_bins
  )
  long_data <- reformatted_output$data
  breakpoints <- reformatted_output$breaks
  bin_sizes <- reformatted_output$bin_length

  # extract weights from long format data structure
  wts_long <- long_data$wts
  long_data[, wts := NULL]

  # fit a HAL regression on the full data set with the CV-selected lambda
  hal_fit <- hal9001::fit_hal(
    X = as.matrix(long_data[, -c(1, 2)]),
    Y = as.numeric(long_data$in_bin),
    max_degree = NULL,
    fit_type = "glmnet",
    family = "binomial",
    lambda = lambda_seq_usm,
    cv_select = T,
    standardize = FALSE, # passed to glmnet
    weights = wts_long, # passed to glmnet
    yolo = FALSE
  )

  # construct output
  out <- list(
    hal_fit = hal_fit,
    breaks = breakpoints,
    bin_sizes = bin_sizes,
    range_a = range(A),
    grid_type_tune_opt = tune_select_params$grid_type,
    n_bins_tune_opt = tune_select_params$n_bins,
    cv_hal_fits_tune_opt = tune_select_fits
  )
  class(out) <- "haldensify"
  return(out)
}



new_fit_hal <- function (X, Y, X_unpenalized = NULL, max_degree = 3, fit_type = c("glmnet",
                                                                   "lassi"), n_folds = 10, foldid = NULL, use_min = TRUE, reduce_basis = NULL,
          family = c("gaussian", "binomial", "cox"), return_lasso = TRUE,
          return_x_basis = FALSE, basis_list = NULL, lambda = NULL,
          id = NULL, offset = NULL, cv_select = TRUE, ..., yolo = TRUE)
{
  call <- match.call(expand.dots = TRUE)
  fit_type <- match.arg(fit_type)
  family <- match.arg(family)
  dot_args <- list(...)
  assertthat::assert_that(!("lambda.min.ratio" %in% names(dot_args) &
                              family == "binomial"), msg = paste("`glmnet` silently ignores",
                                                                 "`lambda.min.ratio` when", "`family = 'binomial'`."))
  assertthat::assert_that(!(fit_type == "lassi" && family ==
                              "binomial"), msg = paste("For binary outcomes, please set",
                                                       "argument 'fit_type' to 'glmnet'."))
  assertthat::assert_that(!(fit_type == "lassi" && family ==
                              "cox"), msg = paste("For Cox models, please set argument",
                                                  "'fit_type' to 'glmnet'."))
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (yolo)
    hal9000()
  if (is.null(foldid)) {
    if (is.null(id)) {
      foldid <- sample(seq_len(n_folds), length(Y), replace = TRUE)
    }
    else {
      unique_ids <- unique(id)
      id_foldid <- sample(seq_len(n_folds), length(unique_ids),
                          replace = TRUE)
      foldid <- id_foldid[match(id, unique_ids)]
    }
  }
  time_start <- proc.time()
  if (is.null(basis_list)) {
    basis_list <- enumerate_basis(X, max_degree)
  }
  col_lists <- unique(lapply(basis_list, `[[`, "cols"))
  col_names <- colnames(X)
  if (!is.null(colnames(X))) {
    col_lists <- lapply(col_lists, function(col_list) col_names[col_list])
  }
  col_lists <- sapply(col_lists, paste, collapse = ",")
  time_enumerate_basis <- proc.time()
  x_basis <- make_design_matrix(X, basis_list)
  time_design_matrix <- proc.time()
  copy_map <- make_copy_map(x_basis)
  unique_columns <- as.numeric(names(copy_map))
  x_basis <- x_basis[, unique_columns]
  penalty_factor <- rep(1, ncol(x_basis))
  unpenalized_covariates <- ifelse(test = is.null(X_unpenalized),
                                   yes = 0, no = {
                                     assertthat::assert_that(is.matrix(X_unpenalized))
                                     assertthat::assert_that(nrow(X_unpenalized) == nrow(x_basis))
                                     ncol(X_unpenalized)
                                   })
  if (unpenalized_covariates > 0) {
    x_basis <- cbind(x_basis, X_unpenalized)
    penalty_factor <- c(penalty_factor, rep(0, ncol(X_unpenalized)))
  }
  time_rm_duplicates <- proc.time()
  if (!is.null(reduce_basis) && is.numeric(reduce_basis)) {
    reduced_basis_map <- make_reduced_basis_map(x_basis,
                                                reduce_basis)
    x_basis <- x_basis[, reduced_basis_map]
  }
  time_reduce_basis <- proc.time()
  if (family == "cox") {
    x_basis <- as.matrix(x_basis)
  }
  if (fit_type == "lassi") {
    message(paste("'lassi' is experimental:", "fit_type='glmnet' is recommended in nearly all cases."))
    hal_lasso <- cv_lasso(x_basis = x_basis, y = Y, n_folds = n_folds)
    if (use_min) {
      lambda_star <- hal_lasso$lambda_min
      coefs <- hal_lasso$betas_mat[, "lambda_min"]
    }
    else {
      lambda_star <- hal_lasso$lambda_1se
      coefs <- hal_lasso$betas_mat[, "lambda_1se"]
    }
  }
  else if (fit_type == "glmnet") {
    if (!cv_select) {
      hal_lasso <- glmnet::glmnet(x = x_basis, y = Y, family = family,
                                  lambda = lambda, penalty.factor = penalty_factor,
                                  ...)
      lambda_star <- hal_lasso$lambda
      coefs <- stats::coef(hal_lasso)
    }
    else {
      hal_lasso <- glmnet::cv.glmnet(x = x_basis, y = Y,
                                     nfolds = n_folds, family = family, lambda = lambda,
                                     foldid = foldid, penalty.factor = penalty_factor,
                                     ...)
      if (use_min) {
        lambda_type <- "lambda.min"
        lambda_star <- hal_lasso$lambda.min
      }
      else {
        lambda_type <- "lambda.1se"
        lambda_star <- hal_lasso$lambda.1se
      }
      coefs <- stats::coef(hal_lasso, lambda_type)
    }
  }
  time_lasso <- proc.time()
  time_final <- proc.time()
  times <- rbind(enumerate_basis = time_enumerate_basis - time_start,
                 design_matrix = time_design_matrix - time_enumerate_basis,
                 remove_duplicates = time_rm_duplicates - time_design_matrix,
                 reduce_basis = time_reduce_basis - time_rm_duplicates,
                 lasso = time_lasso - time_rm_duplicates, total = time_final -
                   time_start)
  fit <- list(call = call, x_basis = if (return_x_basis) {
    x_basis
  } else {
    NULL
  }, basis_list = basis_list, col_lists = col_lists, copy_map = copy_map,
  coefs = as.matrix(coefs), times = times, lambda_star = lambda_star,
  reduce_basis = reduce_basis, family = family, hal_lasso = if (cv_select &
                                                                return_lasso) {
    hal_lasso
  } else {
    NULL
  }, glmnet_lasso = if (!cv_select & return_lasso) {
    hal_lasso
  } else if (cv_select & return_lasso) {
    hal_lasso$glmnet.fit
  } else {
    NULL
  }, unpenalized_covariates = unpenalized_covariates)
  class(fit) <- "hal9001"
  return(fit)
}







###############################################################################

#' Fit conditional density estimation for a sequence of HAL models
#'
#' @details Estimation of the conditional density A|W via a cross-validated
#'  highly adaptive lasso, used to estimate the conditional hazard of failure
#'  in a given bin over the support of A.
#'
#' @param A The \code{numeric} vector of observed values.
#' @param W A \code{data.frame}, \code{matrix}, or similar giving the values of
#'  baseline covariates (potential confounders) for the observed units. These
#'  make up the conditioning set for the conditional density estimate.
#' @param wts A \code{numeric} vector of observation-level weights. The default
#'  is to weight all observations equally.
#' @param grid_type A \code{character} indicating the strategy to be used in
#'  creating bins along the observed support of \code{A}. For bins of equal
#'  range, use \code{"equal_range"}; consult the documentation of
#'  \code{\link[ggplot2]{cut_interval}} for more information. To ensure each
#'  bin has the same number of observations, use \code{"equal_mass"}; consult
#'  the documentation of \code{\link[ggplot2]{cut_number}} for details.
#' @param n_bins This \code{numeric} value indicates the number(s) of bins into
#'  which the support of \code{A} is to be divided.
#' @param cv_folds A \code{numeric} indicating the number of cross-validation
#'  folds to be used in fitting the sequence of HAL conditional density models.
#' @param lambda_seq A \code{numeric} sequence of values of the regularization
#'  parameter of Lasso regression; passed to \code{\link[hal9001]{fit_hal}}.
#'
#' @importFrom data.table ":="
#' @importFrom matrixStats colMeans2
#' @importFrom origami make_folds cross_validate
#'
#' @return A \code{list}, containing density predictions for the sequence of
#'  fitted HAL models; the index and value of the L1 regularization parameter
#'  minimizing the density loss; and the sequence of empirical risks for the
#'  sequence of fitted HAL models.
#'
#' @export
#'
#' @examples
#' # simulate data: W ~ U[-4, 4] and A|W ~ N(mu = W, sd = 0.5)
#' n_train <- 50
#' w <- runif(n_train, -4, 4)
#' a <- rnorm(n_train, w, 0.5)
#' # fit cross-validated HAL-based density estimator of A|W
#' fit_cv_haldensify <- fit_haldensify(
#'   A = a, W = w, n_bins = 3,
#'   lambda_seq = exp(seq(-1, -10, length = 50))
#' )
fit_haldensify <- function(A, W,
                           wts = rep(1, length(A)),
                           grid_type = "equal_range",
                           n_bins = 20,
                           cv_folds = 5,
                           lambda_seq = exp(seq(-1, -13, length = 1000))) {
  # re-format input data into long hazards structure
  reformatted_output <- format_long_hazards(
    A = A, W = W, wts = wts,
    grid_type = grid_type, n_bins = n_bins
  )
  long_data <- reformatted_output$data
  bin_sizes <- reformatted_output$bin_length

  # extract weights from long format data structure
  wts_long <- long_data$wts
  long_data[, wts := NULL]

  # make folds with origami
  folds <- origami::make_folds(long_data,
                               V = cv_folds,
                               cluster_ids = long_data$obs_id
  )

  # call cross_validate on cv_density function
  haldensity <- origami::cross_validate(
    cv_fun = cv_haldensify,
    folds = folds,
    long_data = long_data,
    wts = wts_long,
    lambda_seq = lambda_seq,
    use_future = FALSE,
    .combine = FALSE
  )

  # re-organize output cross-validation procedure
  density_pred_unscaled <- do.call(rbind, as.list(haldensity$preds))

  # re-scale predictions by multiplying by bin width for each failure bin
  density_pred_scaled <- apply(density_pred_unscaled, 2, function(x) {
    pred <- x / bin_sizes[long_data[in_bin == 1, bin_id]]
    return(pred)
  })
  obs_wts <- do.call(c, as.list(haldensity$wts))

  # compute loss for the given individual
  density_loss <- apply(density_pred_scaled, 2, function(x) {
    pred_weighted <- x * obs_wts
    loss_weighted <- -log(pred_weighted)
    return(loss_weighted)
  })

  # take column means to have average loss across sequence of lambdas
  emp_risks_density_loss <- matrixStats::colMeans2(density_loss)

  # find minimizer of loss in lambda sequence
  lambda_loss_min_idx <- which.min(emp_risks_density_loss)
  lambda_loss_min <- lambda_seq[lambda_loss_min_idx]

  # return loss minimizer in lambda, Pn losses, and all density estimates
  out <- list(
    lambda_loss_min_idx = lambda_loss_min_idx,
    lambda_loss_min = lambda_loss_min,
    emp_risks = emp_risks_density_loss,
    density_pred = density_pred_scaled,
    lambda_seq = lambda_seq
  )
  return(out)
}
