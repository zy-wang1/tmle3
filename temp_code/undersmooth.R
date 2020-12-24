n_mode <- 1
# generate true discrete density
n_sim <- 1e6
# set.seed(1234)
data_out <- simulate_data(n_sim = n_sim, n_mode = n_mode)
# truth <- mean(data_out$f_d(data_out$x))
# truth
# data_out$x %>% hist(freq = F, ylim = c(0, 0.5))
data_out$x %>% density %>% plot(ylim = c(0, 0.8))
cdf_obj <- data_out$x %>% ecdf
bin_width <- 0.5
support_points <- seq(-5, 5, by = bin_width)
target_density <- sapply(1:length(support_points), function(i) {
  if (i == 1) cdf_obj(support_points[i]) else cdf_obj(support_points[i]) - cdf_obj(support_points[i - 1])
})
target_density[length(target_density)] <- 1 - sum(target_density[1:(length(target_density) - 1)])
true_pmf <- data.frame(x = support_points,
                       p = target_density)
true_pdf <- data.frame(x = support_points,
                       p = target_density/bin_width)
true_pmf_1_mode <- true_pmf


B <- 300
library(parallel)
set.seed(123)

report_list <- mclapply(1:8, mc.cores = 8, function(loc_bias) {
  delta_bias <- seq(0.01, 0.08, by = 0.01)[loc_bias]
  n_mode <- 4
  # generate true discrete density
  n_sim <- 1e6
  # set.seed(1234)
  data_out <- simulate_data(n_sim = n_sim, n_mode = n_mode)
  # truth <- mean(data_out$f_d(data_out$x))
  # truth
  # data_out$x %>% hist(freq = F, ylim = c(0, 0.5))
  data_out$x %>% density %>% plot(ylim = c(0, 0.8))
  cdf_obj <- data_out$x %>% ecdf
  bin_width <- 0.5
  support_points <- seq(-5, 5, by = bin_width)
  target_density <- sapply(1:length(support_points), function(i) {
    if (i == 1) cdf_obj(support_points[i]) else cdf_obj(support_points[i]) - cdf_obj(support_points[i - 1])
  })
  target_density[length(target_density)] <- 1 - sum(target_density[1:(length(target_density) - 1)])
  true_pmf <- data.frame(x = support_points,
                         p = target_density)
  true_pdf <- data.frame(x = support_points,
                         p = target_density/bin_width)

n_iteration <- 100
record <- list()
for (iteration in 1:n_iteration) {

obs_x <- sample(support_points, B, replace = T, prob = target_density)  # simulated data

# empirical
np_pmf <- data.frame(
  x = support_points,
  p = sapply(support_points, function(x) sum(obs_x == x))/B
  # (cut_interval(obs_x, n = n_bins, right = F) %>% table)/B
)

# hal
n_bins <- diff(range(obs_x))/bin_width+1  # make sure each HAL bin contains one support point
temp_cut <- ggplot2::cut_interval(obs_x, n = n_bins, right = F, ordered_result = TRUE, dig.lab = 12)
# temp_cut %>% levels %>% length
# diff(range(obs_x))/bin_width+1
example_data <- data.frame(W = 0, A = obs_x)
density_hal <- haldensify(
  A = example_data$A,
  W = example_data$W,
  n_bins = n_bins,
  grid_type = "equal_range",
  lambda_seq = exp(seq(-1, -10, length = 100)), use_future = F  # lmabda_seq is not used; let cv.glmnet search lambda
)

# density_hal$hal_fit$hal_lasso$lambda.min
# density_hal$hal_fit$lambda_star

undersmooth_seq <- seq(density_hal$hal_fit$lambda_star, 0, length.out = 11)
undersmooth_seq <- undersmooth_seq[-length(undersmooth_seq)]
temp_D1_vec <- sapply(undersmooth_seq, function(temp_lambda) {
  left_breaks <- density_hal$breaks[-length(density_hal$breaks)]
  support_points_used <- support_points[support_points <= range(obs_x)[2] & support_points >= range(obs_x)[1]]
  # get preds for all support_points_used at temp_lambda
  {
    new_A <- support_points_used
    new_W <- rep(0, length(new_A))
    object <- density_hal
    hal_fit_object <- density_hal$hal_fit
    offset <- NULL

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

    new_data =
      long_data_pred[, 3:ncol(long_data_pred)]
    new_data <- as.matrix(new_data)

    # generate design matrix
    pred_x_basis <- hal9001::make_design_matrix(new_data, hal_fit_object$basis_list)
    group <- hal_fit_object$copy_map[[1]]

    # reduce matrix of basis functions
    pred_x_basis <- hal9001::apply_copy_map(pred_x_basis, hal_fit_object$copy_map)

    new_X_unpenalized = NULL
    # add unpenalized covariates
    new_unpenalized_covariates <- ifelse(
      test = is.null(new_X_unpenalized),
      yes = 0,
      no = {
        assert_that(is.matrix(new_X_unpenalized))
        assert_that(nrow(new_X_unpenalized) == nrow(new_data))
        ncol(new_X_unpenalized)
      }
    )

    # column rank of X_unpenalized should be consistent between the prediction
    # and training phases
    assertthat::assert_that(hal_fit_object$unpenalized_covariates ==
                              new_unpenalized_covariates)
    if (new_unpenalized_covariates > 0) {
      pred_x_basis <- cbind(pred_x_basis, new_X_unpenalized)
    }

    # hal_fit_object$coefs
    # coef(hal_fit_object$glmnet_lasso, s = object$hal_fit$lambda_star)
    temp_coef <- coef(hal_fit_object$glmnet_lasso, s = temp_lambda)


    # generate predictions
    if (hal_fit_object$family != "cox") {
      if (ncol(temp_coef) > 1) {
        preds <- apply(temp_coef, 2, function(hal_coefs) {
          as.vector(Matrix::tcrossprod(
            x = pred_x_basis,
            y = hal_coefs[-1]
          ) + temp_coef[1])
        })
      } else {
        preds <- as.vector(Matrix::tcrossprod(
          x = pred_x_basis,
          y = temp_coef[-1]
        ) + temp_coef[1])
      }
    } else {
      # Note: there is no intercept in the Cox model (built into the baseline
      #       hazard and would cancel in the partial likelihood).
      # message(paste("The Cox Model is not commonly used for prediction,",
      # "proceed with caution."))
      if (ncol(hal_fit_object$coefs) > 1) {
        preds <- apply(hal_fit_object$coefs, 2, function(hal_coefs) {
          as.vector(Matrix::tcrossprod(
            x = pred_x_basis,
            y = hal_coefs
          ))
        })
      } else {
        preds <- as.vector(Matrix::tcrossprod(
          x = pred_x_basis,
          y = as.vector(hal_fit_object$coefs)
        ))
      }
    }

    if (!is.null(offset)) {
      preds <- preds + offset
    }

    # apply logit transformation for logistic regression predictions
    if (hal_fit_object$family == "binomial") {
      preds <- stats::plogis(preds)
    } else if (hal_fit_object$family == "cox") {
      preds <- exp(preds)
    }

    hazard_pred <- preds

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
    # if (cv_select) {
    density_pred_rescaled <- density_pred_rescaled[, 1]
    # }
  }

  temp_hal_pmf <- data.frame(x = support_points, p = 0)
  temp_hal_pmf[support_points %in% support_points_used, "p"] <- density_pred_rescaled * density_hal$bin_sizes

  temp_int_mat <- temp_hal_pmf %>% mutate(D1 = (2 * (p - sum(p^2))),
                          p_obs = np_pmf$p,
                          to_sum = D1*p_obs)
  sum(temp_int_mat$to_sum) %>% return
})
undersmoothed_lambda <- undersmooth_seq[first(which(abs(temp_D1_vec) < 1/B))]

{
  left_breaks <- density_hal$breaks[-length(density_hal$breaks)]
  support_points_used <- support_points[support_points <= range(obs_x)[2] & support_points >= range(obs_x)[1]]
  # get preds for all support_points_used at temp_lambda
  temp_lambda <- undersmoothed_lambda
  {
    new_A <- support_points_used
    new_W <- rep(0, length(new_A))
    object <- density_hal
    hal_fit_object <- density_hal$hal_fit
    offset <- NULL

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

    new_data =
      long_data_pred[, 3:ncol(long_data_pred)]
    new_data <- as.matrix(new_data)

    # generate design matrix
    pred_x_basis <- hal9001::make_design_matrix(new_data, hal_fit_object$basis_list)
    group <- hal_fit_object$copy_map[[1]]

    # reduce matrix of basis functions
    pred_x_basis <- hal9001::apply_copy_map(pred_x_basis, hal_fit_object$copy_map)

    new_X_unpenalized = NULL
    # add unpenalized covariates
    new_unpenalized_covariates <- ifelse(
      test = is.null(new_X_unpenalized),
      yes = 0,
      no = {
        assert_that(is.matrix(new_X_unpenalized))
        assert_that(nrow(new_X_unpenalized) == nrow(new_data))
        ncol(new_X_unpenalized)
      }
    )

    # column rank of X_unpenalized should be consistent between the prediction
    # and training phases
    assertthat::assert_that(hal_fit_object$unpenalized_covariates ==
                              new_unpenalized_covariates)
    if (new_unpenalized_covariates > 0) {
      pred_x_basis <- cbind(pred_x_basis, new_X_unpenalized)
    }

    # hal_fit_object$coefs
    # coef(hal_fit_object$glmnet_lasso, s = object$hal_fit$lambda_star)
    temp_coef <- coef(hal_fit_object$glmnet_lasso, s = temp_lambda)


    # generate predictions
    if (hal_fit_object$family != "cox") {
      if (ncol(temp_coef) > 1) {
        preds <- apply(temp_coef, 2, function(hal_coefs) {
          as.vector(Matrix::tcrossprod(
            x = pred_x_basis,
            y = hal_coefs[-1]
          ) + temp_coef[1])
        })
      } else {
        preds <- as.vector(Matrix::tcrossprod(
          x = pred_x_basis,
          y = temp_coef[-1]
        ) + temp_coef[1])
      }
    } else {
      # Note: there is no intercept in the Cox model (built into the baseline
      #       hazard and would cancel in the partial likelihood).
      # message(paste("The Cox Model is not commonly used for prediction,",
      # "proceed with caution."))
      if (ncol(hal_fit_object$coefs) > 1) {
        preds <- apply(hal_fit_object$coefs, 2, function(hal_coefs) {
          as.vector(Matrix::tcrossprod(
            x = pred_x_basis,
            y = hal_coefs
          ))
        })
      } else {
        preds <- as.vector(Matrix::tcrossprod(
          x = pred_x_basis,
          y = as.vector(hal_fit_object$coefs)
        ))
      }
    }

    if (!is.null(offset)) {
      preds <- preds + offset
    }

    # apply logit transformation for logistic regression predictions
    if (hal_fit_object$family == "binomial") {
      preds <- stats::plogis(preds)
    } else if (hal_fit_object$family == "cox") {
      preds <- exp(preds)
    }

    hazard_pred <- preds

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
    # if (cv_select) {
    density_pred_rescaled <- density_pred_rescaled[, 1]
    # }
  }

  temp_hal_pmf <- data.frame(x = support_points, p = 0)
  temp_hal_pmf[support_points %in% support_points_used, "p"] <- density_pred_rescaled * density_hal$bin_sizes

}
hal_pmf <- temp_hal_pmf  # now it is undersmoothed
hal_pdf <- data.frame(x = hal_pmf$x,
                      p = hal_pmf$p / bin_width)

# biased initial
temp_p <- np_pmf$p * (8-loc_bias) / 8 + true_pmf_1_mode$p * loc_bias / 8
temp_p <- temp_p/sum(temp_p)
biased_pmf <- data.frame(x = support_points, p = temp_p)

# Pn update
updated_pmf <- biased_pmf
record_epsilon <- c()
record_D <- c()
for (step in 1:300) {
  loss_empirical <- function(epsilon) {
    temp_df_update <- data.frame(np_pmf, initial_p = updated_pmf$p) %>%
      mutate(D1 = 2 * (initial_p - sum(initial_p^2))) %>%
      mutate(updated_p = (1 + epsilon * D1) * initial_p) %>%
      mutate(to_sum = ifelse(updated_p != 0, p * (-log(updated_p)), 0) )
    sum(temp_df_update$to_sum) %>% return()
  }
  optim_fit <- optim(
    par = list(epsilon = 0.1), fn = loss_empirical,
    lower = -0.1, upper = 0.1,
    method = "Brent"
  )
  epsilon <- optim_fit$par
  temp_df_update <- data.frame(np_pmf, initial_p = updated_pmf$p) %>%
    mutate(D1 = 2 * (initial_p - sum(initial_p^2))) %>%
    mutate(updated_p = (1 + epsilon * D1) * initial_p) %>%
    mutate(to_sum = ifelse(updated_p != 0, p * (-log(updated_p)), 0) )
  updated_pmf <- data.frame(x = temp_df_update$x, p = temp_df_update$updated_p)
  # decide if to stop
  record_epsilon[step] <- epsilon
  record_D[step] <- sum(temp_df_update$p * temp_df_update$D1)
  if (!(all(record_epsilon > 0) | all(record_epsilon < 0)) | abs(record_D[step]) < 1/length(obs_x)) break
}
step_Pn <- step
updated_pmf_Pn <- updated_pmf

# tilde Pn update
updated_pmf <- biased_pmf
record_epsilon <- c()
record_D <- c()
for (step in 1:300) {
  loss <- function(epsilon) {
    temp_df_update <- data.frame(hal_pmf, initial_p = updated_pmf$p) %>%
      mutate(D1 = 2 * (initial_p - sum(initial_p^2))) %>%
      mutate(updated_p = (1 + epsilon * D1) * initial_p) %>%
      mutate(to_sum = ifelse(updated_p != 0, p * (-log(updated_p)), 0) )
    sum(temp_df_update$to_sum) %>% return()  # expectation w.r.t. tilde Pn
  }
  optim_fit <- optim(
    par = list(epsilon = 0.1), fn = loss,
    lower = -0.1, upper = 0.1,
    method = "Brent"
  )
  epsilon <- optim_fit$par
  temp_df_update <- data.frame(hal_pmf, initial_p = updated_pmf$p) %>%
    mutate(D1 = 2 * (initial_p - sum(initial_p^2))) %>%
    mutate(updated_p = (1 + epsilon * D1) * initial_p) %>%
    mutate(to_sum = ifelse(updated_p != 0, p * (-log(updated_p)), 0) )
  updated_pmf <- data.frame(x = temp_df_update$x, p = temp_df_update$updated_p)
  # decide if to stop
  record_epsilon[step] <- epsilon
  record_D[step] <- sum(temp_df_update$p * temp_df_update$D1)
  if (!(all(record_epsilon > 0) | all(record_epsilon < 0)) | abs(record_D[step]) < 1/length(obs_x)) break
}
step_tildePn <- step
updated_pmf_tildePn <- updated_pmf

# 2TMLE
{
  # tilde Pn update
  updated_pmf_1 <- updated_pmf_2 <- updated_pmf <- biased_pmf
  record_epsilon <- record_epsilon_2 <- c()
  record_D_1 <- record_D_2 <- c()
  for (step in 1:300) {
    # first order update
    loss <- function(epsilon) {
      temp_df_update <- data.frame(hal_pmf, initial_p = updated_pmf_1$p) %>%
        mutate(D1 = 2 * (initial_p - sum(initial_p^2))) %>%
        mutate(updated_p = (1 + epsilon * D1) * initial_p) %>%
        mutate(to_sum = ifelse(updated_p != 0, p * (-log(updated_p)), 0) )
      sum(temp_df_update$to_sum) %>% return()  # expectation w.r.t. tilde Pn
    }
    optim_fit <- optim(
      par = list(epsilon = 0.1), fn = loss,
      lower = -0.1, upper = 0.1,
      method = "Brent"
    )
    epsilon <- optim_fit$par
    epsilon_1 <- epsilon
    temp_df_update <- data.frame(hal_pmf, initial_p = updated_pmf_1$p) %>%
      mutate(D1 = 2 * (initial_p - sum(initial_p^2))) %>%
      mutate(updated_p = (1 + epsilon_1 * D1) * initial_p) %>%
      mutate(to_sum = ifelse(updated_p != 0, p * (-log(updated_p)), 0) )
    updated_pmf_1 <- data.frame(x = temp_df_update$x, p = temp_df_update$updated_p)  # might be used below to get D2
    # second order update

    # everything not depending on epsilon here
    # get cnP
    dat_int_tilde <- temp_df_update %>%
      mutate(value = D1^2 / (1 + epsilon_1 * D1)^2) %>%
      mutate(cnP_to_sum = p * value)
    cnP <- sum(dat_int_tilde$cnP_to_sum)
    # get dP
    dat_int_P <- data.frame(updated_pmf_2, initial_p = updated_pmf_2$p) %>%
      mutate(D1 = 2 * (initial_p - sum(initial_p^2))) %>%
      mutate(updated_p = (1 + epsilon_1 * D1) * initial_p) %>%
      mutate(D1P1 = 2 * (updated_p - sum(updated_p^2))) %>%
      mutate(to_sum = p*D1*D1P1)
    dP <- sum(dat_int_P$to_sum) / cnP
    # deside if there is need to do a 2nd order update
    temp_df_update_2 <- temp_df_update %>%
      mutate(D1P1 = 2 * (updated_p - sum(updated_p^2))) %>%
      mutate(line1 = D1P1 * (1 + epsilon_1 * D1)) %>%
      mutate(to_sum_denominator = p / (1 + epsilon_1 * D1)^2) %>%
      mutate(line2 = 2 * dP * (p / (1 + epsilon_1 * D1)^2 - 2*initial_p*sum(to_sum_denominator)) ) %>%
      mutate(to_sum_line3 = initial_p * p / (1 + epsilon_1 * D1)^2 ) %>%
      mutate(line3 = -2 * dP * sum(to_sum_line3)) %>%
      mutate(line4 = 4 * dP * sum(initial_p^2) * sum(to_sum_denominator)) %>%
      mutate(to_sum_line5 = initial_p * D1P1) %>%
      mutate(line5 = 2 * epsilon_1 * initial_p * (D1P1 - 2 * sum(to_sum_line5))) %>%
      mutate(to_sum_line6 = initial_p * initial_p * (D1P1 - 2 * sum(to_sum_line5))) %>%
      mutate(line6 = -2 * epsilon_1 * sum(to_sum_line6)) %>%
      mutate(D2 = line1+line2+line3+line4+line5+line6)
    if (abs(sum(temp_df_update_2$p * temp_df_update_2$D2)) < 1/length(obs_x)) {
      epsilon_2 <- 0
      temp_df_update_2 <- temp_df_update_2 %>%
        mutate(updated_p_2 =
                 # (1 + epsilon_2 * D2) *
                 initial_p) %>%
        mutate(to_sum_final = ifelse(updated_p_2 != 0, p * (-log(updated_p_2)), 0) )
      updated_pmf_2 <- data.frame(x = temp_df_update_2$x,
                                  p = temp_df_update_2$updated_p_2)
    } else {
      # new loss
      loss_2 <- function(epsilon) {
        temp_df_update_2 <- temp_df_update_2 %>%
          mutate(updated_p_2 = (1 + epsilon * D2) * initial_p) %>%
          mutate(to_sum_final = ifelse(updated_p_2 != 0, p * (-log(updated_p_2)), 0) )
        sum(temp_df_update_2$to_sum_final) %>% return()  # expectation w.r.t. tilde Pn
      }
      optim_fit <- optim(
        par = list(epsilon = 0.1), fn = loss_2,
        lower = -0.1, upper = 0.1,
        method = "Brent"
      )
      epsilon_2 <- optim_fit$par
      temp_df_update_2 <- temp_df_update_2 %>%
        mutate(updated_p_2 = (1 + epsilon_2 * D2) * initial_p) %>%
        mutate(to_sum_final = ifelse(updated_p_2 != 0, p * (-log(updated_p_2)), 0) )
      updated_pmf_2 <- data.frame(x = temp_df_update_2$x,
                                  p = temp_df_update_2$updated_p_2)
    }
    record_epsilon_2[step] <- epsilon_2

    # final update
    updated_pmf <- updated_pmf_2
    loss <- function(epsilon) {
      temp_df_update <- data.frame(hal_pmf, initial_p = updated_pmf$p) %>%
        mutate(D1 = 2 * (initial_p - sum(initial_p^2))) %>%
        mutate(updated_p = (1 + epsilon * D1) * initial_p) %>%
        mutate(to_sum = ifelse(updated_p != 0, p * (-log(updated_p)), 0) )
      sum(temp_df_update$to_sum) %>% return()  # expectation w.r.t. tilde Pn
    }
    optim_fit <- optim(
      par = list(epsilon = 0.1), fn = loss,
      lower = -0.1, upper = 0.1,
      method = "Brent"
    )
    epsilon <- optim_fit$par
    temp_df_update <- data.frame(hal_pmf, initial_p = updated_pmf$p) %>%
      mutate(D1 = 2 * (initial_p - sum(initial_p^2))) %>%
      mutate(updated_p = (1 + epsilon * D1) * initial_p) %>%
      mutate(to_sum = ifelse(updated_p != 0, p * (-log(updated_p)), 0) )
    updated_pmf <- data.frame(x = temp_df_update$x, p = temp_df_update$updated_p)  # might be used below to get D2

    # stop if both D1 and D2 equations are solved
    record_epsilon[step] <- epsilon
    record_D_2[step] <- sum(temp_df_update_2$p * temp_df_update_2$D2)
    record_D_1[step] <- sum(temp_df_update$p * temp_df_update$D1)
    if (
      # !(all(record_epsilon > 0) | all(record_epsilon < 0)) |
      (abs(record_D_1[step]) < 1/length(obs_x) & abs(record_D_2[step]) < 1/length(obs_x)) ) break else
        updated_pmf_1 <- updated_pmf_2 <- updated_pmf
  }
  step_tildePn_2 <- step
  updated_pmf_tildePn_2 <- updated_pmf
  ratio_2 <- sum(record_epsilon_2 != 0)/length(record_epsilon_2)
}

record[[iteration]] <- c(sum(np_pmf$p^2)/bin_width,
                         sum(hal_pmf$p^2)/bin_width,
                         sum(biased_pmf$p^2)/bin_width,
                         sum(updated_pmf_Pn$p^2)/bin_width,
                         sum(updated_pmf_tildePn$p^2)/bin_width,
                         sum(updated_pmf_tildePn_2$p^2)/bin_width,
                         step_Pn,
                         step_tildePn,
                         step_tildePn_2,
                         ratio_2
)

}

truth <- sum(true_pmf$p^2)/bin_width

record_mat <- do.call(rbind, record)
report <- data.frame(Bias = colMeans(record_mat[,1:6] - truth),
                     SD = record_mat[,1:6] %>% apply(2, sd),
                     MSE = colMeans((record_mat[,1:6] - truth)^2),
                     Step = c(rep(0, 3), colMeans(record_mat[,7:9])),
                     Ratio = c(rep(0, 5), mean(record_mat[,10]))
)
rownames(report) <- c("NP-MLE", "HAL-MLE", "Biased initial",  "Pn 1-TMLE", "tilde Pn 1-TMLE", "tilde Pn 2-TMLE")

return(report)
})


# saveRDS(report_list, "./temp_output/higher_order_linear_combination_undersmooth_ZW-20201223.rds")






plot_df <- report_list[-8] %>% lapply(function(x) x[, 1]) %>% abind::abind(along = 0)
colnames(plot_df) <- c("NP-MLE", "HAL-MLE", "Biased Initial", "Pn 1TMLE", "Tilde 1TMLE", "Tilde 2TMLE")
plot_df <- data.frame(plot_df)
plot_df[["TV"]] <- record_TV[3:9, 4]
plot_df[["Bias.Mass"]] <- seq(0.01, 0.07, 0.01) %>% as.factor()
plot_df <- plot_df[
  # (1:4)*2-1
  , -c(1:2, 7)]

library(reshape2)
library(ggrepel)
test_data_long <- reshape2::melt(plot_df, id="Bias.Mass")  # convert to long format

plot1 <- ggplot(data=test_data_long,
                aes(x=Bias.Mass, y=value, colour=variable, label = variable, group = variable)) +
  geom_line() + geom_text_repel(data = test_data_long %>% filter(Bias.Mass == 0.01), show.legend = FALSE) + geom_point() + theme_classic() + geom_hline(yintercept = 0, linetype="dashed", color = "black") +
  labs(y = "Bias", x = "Bias Mass", color = "Estimator") +
  scale_x_discrete(name ="Bias Mass",
                   breaks = plot_df$Bias.Mass)
plot1
ggsave("./temp_output/bias_linear_combination_undersmooth.jpg", plot1, width = 8, height = 8, units = "in", dpi = 300)




plot_df <- report_list[-8] %>% lapply(function(x) abs(x[, 1])/x[, 2]) %>% abind::abind(along = 0)
colnames(plot_df) <- c("NP-MLE", "HAL-MLE", "Biased Initial", "Pn 1TMLE", "Tilde 1TMLE", "Tilde 2TMLE")
plot_df <- data.frame(plot_df)
plot_df[["TV"]] <- record_TV[3:9, 4]
plot_df[["Bias.Mass"]] <- seq(0.01, 0.07, 0.01) %>% as.factor()
plot_df <- plot_df[
  # (1:4)*2-1
  , -c(1, 2, 3, 7)]

library(reshape2)
library(ggrepel)
test_data_long <- reshape2::melt(plot_df, id="Bias.Mass")  # convert to long format

plot1 <- ggplot(data=test_data_long,
                aes(x=Bias.Mass, y=value, colour=variable, label = variable, group = variable)) +
  geom_line() + geom_text_repel(data = test_data_long %>% filter(Bias.Mass == 0.01), show.legend = FALSE) + geom_point() + theme_classic() + geom_hline(yintercept = 0, linetype="dashed", color = "black") +
  labs(y = "Bias/SD Ratio", x = "Bias Mass", color = "Estimator") +
  scale_x_discrete(name ="Bias Mass",
                   breaks = plot_df$Bias.Mass)
plot1
ggsave("./temp_output/ratio_linear_combination_undersmooth.jpg", plot1, width = 12, height = 8, units = "in", dpi = 300)








plot_df <- report_list[-8] %>% lapply(function(x) x[, 3]) %>% abind::abind(along = 0)
colnames(plot_df) <- c("NP-MLE", "HAL-MLE", "Biased Initial", "Pn 1TMLE", "Tilde 1TMLE", "Tilde 2TMLE")
plot_df <- data.frame(plot_df)
plot_df[["TV"]] <- record_TV[3:9, 4]
plot_df[["Bias.Mass"]] <- seq(0.01, 0.07, 0.01) %>% as.factor()
plot_df <- plot_df[
  # (1:4)*2-1
  , -c(1:3, 7)]

library(reshape2)
library(ggrepel)
test_data_long <- reshape2::melt(plot_df, id="Bias.Mass")  # convert to long format

plot1 <- ggplot(data=test_data_long,
                aes(x=Bias.Mass, y=value, colour=variable, label = variable, group = variable)) +
  geom_line() + geom_text_repel(data = test_data_long %>% filter(Bias.Mass == 0.01), show.legend = FALSE) + geom_point() + theme_classic() + geom_hline(yintercept = 0, linetype="dashed", color = "black") +
  labs(y = "MSE", x = "Bias Mass", color = "Estimator") +
  scale_x_discrete(name ="Bias Mass",
                   breaks = plot_df$Bias.Mass)
plot1
ggsave("./temp_output/mse_linear_combination_undersmooth.jpg", plot1, width = 8, height = 8, units = "in", dpi = 300)

