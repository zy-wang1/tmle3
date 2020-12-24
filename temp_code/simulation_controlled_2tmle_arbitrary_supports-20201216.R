# library(R6)  # R6class
library(data.table)  # setDT
# library(sl3)
# library(digest)
# library(uuid)  # UUIDgenerate
# library(delayed)  # bundle_delayed
# library(assertthat)  # assert_that
# library(speedglm)  # speedglm
# library(methods)  # is

library(dplyr)  # dplyr::select, to mask other packages
# library(haldensify)  # haldensify
library(ggplot2)  # ggplot
library(xtable)

# code_list <- list.files("./R", full.names = T)
# for (code in code_list) source(code)
source("./temp_code/haldensify.R")

###
# generate target density by discretizing Cai's example
###
n_mode <- 10
{
  simulate_data <- function(n_sim, n_mode) {
    modes <- seq(from = -4, to = 4, length.out = n_mode)
    sigma <- 10 / n_mode / 6
    # sample object
    x <- rnorm(n = n_sim, mean = sample(x = modes, size = n_sim, replace = TRUE), sd = sigma)
    # function object
    dmixture <- function(x, modes, sigma) {
      all_d <- c()
      for (mode_here in modes) all_d <- c(all_d, dnorm(x, mean = mode_here, sd = sigma))
      return(mean(all_d))
    }
    foo2 <- Vectorize(function(x) dmixture(x = x, modes = modes, sigma = sigma))
    foo <- Vectorize(function(x) dmixture(x = x, modes = modes, sigma = sigma)^2)
    return(list(x = x, f_d = foo2, f_dsquare = foo))
  }

  # generate true discrete density
  n_sim <- 1e6
  set.seed(1234)
  data_out <- simulate_data(n_sim = n_sim, n_mode = n_mode)
  # truth <- mean(data_out$f_d(data_out$x))
  # truth
  # data_out$x %>% hist(freq = F, ylim = c(0, 0.5))
  data_out$x %>% density %>% plot(ylim = c(0, 0.8))
  cdf_obj <- data_out$x %>% ecdf
  bin_width <- 1
  support_points <- seq(-5, 5, by = bin_width)
  target_density <- sapply(1:length(support_points), function(i) {
    if (i == 1) cdf_obj(support_points[i]) else cdf_obj(support_points[i]) - cdf_obj(support_points[i - 1])
  })
  target_density[length(target_density)] <- 1 - sum(target_density[1:(length(target_density) - 1)])
  true_pmf <- data.frame(x = support_points,
                         p = target_density)
}
# truth is the expectation of p(x)
# true_pmf
truth <- sum(true_pmf$p^2)
sum(true_pmf$p^2)
truth
true_pmf
support_points

# manually created bias
{
  set.seed(123)
  B <- 300
  n_iteration <- 50
  record <- matrix(0, n_iteration, 3)
  for (i in 1:n_iteration) {
    obs_x <- sample(support_points, B, replace = T, prob = target_density)  # simulated data
    # n_bins <- length(unique(obs_x))
    n_bins <- diff(obs_x %>% range) / bin_width + 1

    ggplot2::cut_interval(obs_x, n_bins,
                          right = FALSE,
                          ordered_result = TRUE, dig.lab = 12
    ) %>% levels
    example_data <- data.frame(W = 0, A = obs_x)
    density_hal <- haldensify(
      A = example_data$A,
      W = example_data$W,
      n_bins = n_bins,
      grid_type = "equal_range",
      lambda_seq = exp(seq(1, -13, length = 100)), use_future = T
    )
    # density_hal$breaks
    # density_hal$hal_fit
    # predict(density_hal, new_A = support_points, new_W = 0)
    fit_width <- diff(obs_x %>% range) / n_bins
    temp_p <- predict(density_hal, new_A = support_points, new_W = 0) * fit_width
    sum(temp_p)
    # temp_p <- temp_p/sum(temp_p)
    sum(predict(density_hal, new_A = support_points, new_W = 0) * fit_width)
    record[i, 1] <- sum((temp_p)^2)  # each predicted value is pdf; need to be transformed to pmf
    record[i, 2] <- sum(((cut_interval(obs_x, n = n_bins, right = F, ordered_result = TRUE, dig.lab = 12) %>% table)/B)^2)
    temp_p <- predict(density_hal, new_A = support_points, new_W = 0) * fit_width
    temp_p <- temp_p + 0.1
    temp_p <- temp_p/sum(temp_p)
    record[i, 3] <- sum(temp_p^2)
  }
  report <- data.frame(bias = (record - truth) %>% colMeans(),
                       sd = record %>% apply(2, sd),
                       MSE = ((record - truth)^2) %>% colMeans()
  )
  rownames(report) <- c("HAL", "Empirical", "Biased")
  report
  report %>%
    xtable(type = "latex",
           caption = paste0("Sample size ", 300, "; iteration: ", 100),
           digits = 5)
}

# one example of hal fit and biased version
{
  support_points <- seq(-5, 6, 1)
  B <- 300
  obs_x <- sample(support_points, B, replace = T, prob = target_density)  # simulated data
  n_bins <- length(unique(obs_x))-1
  # # empirical
  # n_bins
  # cut_interval(obs_x, n = n_bins, right = F) %>% levels
  # data.frame(
  #   est = (cut_interval(obs_x, n = n_bins, right = F) %>% table)/B,
  #   true = true_pmf$p[-length(true_pmf$p)]
  # )

  example_data <- data.frame(W = 0, A = obs_x)
  density_hal <- haldensify(
    A = example_data$A,
    W = example_data$W,
    n_bins = n_bins,
    grid_type = "equal_range",
    lambda_seq = exp(seq(-1, -10, length = 100)), use_future = T
  )
  data_out$x %>% density %>% plot(ylim = c(0, 0.8), main = "")
  true_pmf_transformed %>% points(pch = 20)
  points(seq(-5, 5, 1),
         predict(density_hal, new_A = seq(-5, 5, 1), new_W = 0),
         col = "blue"
  )
  temp_p <- predict(density_hal, new_A = seq(-5, 5, 1), new_W = 0)
  temp_p[length(temp_p)] <- min(1, 3*temp_p[length(temp_p)])
  ratio <- (1 - temp_p[length(temp_p)]) / sum(temp_p[-length(temp_p)])
  temp_p[-length(temp_p)] <- temp_p[-length(temp_p)] * ratio
  points(seq(-5, 5, 1),
         temp_p,
         col = "red")
  legend("topleft", legend=c("Discretized true pmf", "HAL", "Biased estimation"),
         col=c("black", "blue", "red"), pch=c(20, 1, 1))

}



# 1-TMLE w.r.t. Pn and tilde Pn
{
  set.seed(123)
  record <- list()
  for (iteration in 1:100) {
    support_points <- seq(-5, 6, 1)  # the last support point is pseudo
    B <- 300
    obs_x <- sample(support_points, B, replace = T, prob = target_density)  # simulated data
    actual_data <- obs_x
    actual_data[actual_data == 6] <- 5
    # empirical
    np_pmf <- data.frame(
      x = seq(-5, 5, 1),
      p = sapply(seq(-5, 5, 1), function(x) sum(actual_data == x))/B
      # (cut_interval(obs_x, n = n_bins, right = F) %>% table)/B
    )
    # hal
    # n_bins <- length(support_points) - 1
    n_bins <- diff(range(actual_data)) + 1  # make sure they are integer intervals, and the last one is [5, 6]
    cut_interval(obs_x, n = n_bins, right = F) %>% table
    example_data <- data.frame(W = 0, A = obs_x)
    density_hal <- haldensify(
      A = example_data$A,
      W = example_data$W,
      n_bins = n_bins,
      grid_type = "equal_range",
      lambda_seq = exp(seq(-1, -10, length = 100)), use_future = T
    )
    hal_pmf <- data.frame(x = seq(-5, 5, 1),
                          p = predict(density_hal, new_A = seq(-5, 5, 1), new_W = 0) * (diff(range(obs_x)) / n_bins)
    )
    # biased initial
    temp_p <- predict(density_hal, new_A = seq(-5, 5, 1), new_W = 0)
    temp_p[length(temp_p)] <- min(1, 3*temp_p[length(temp_p)])
    ratio <- (1 - temp_p[length(temp_p)]) / sum(temp_p[-length(temp_p)])
    temp_p[-length(temp_p)] <- temp_p[-length(temp_p)] * ratio
    biased_pmf <- data.frame(x = seq(-5, 5, 1), p = temp_p)
    # Pn update
    updated_pmf <- biased_pmf
    record_epsilon <- c()
    record_D <- c()
    for (step in 1:100) {
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
    for (step in 1:100) {
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

    record[[iteration]] <- c(sum(biased_pmf$p^2),
                             sum(updated_pmf_Pn$p^2),
                             sum(updated_pmf_tildePn$p^2),
                             step_Pn,
                             step_tildePn)
  }
  record_mat <- do.call(rbind, record)
  report <- data.frame(Bias = colMeans(record_mat[,1:3] - truth),
                       SD = record_mat[,1:3] %>% apply(2, sd),
                       MSE = colMeans((record_mat[,1:3] - truth)^2),
                       Step = c(0, colMeans(record_mat[,4:5]))
  )
  rownames(report) <- c("Biased initial", "Pn update", "tilde Pn update")
  report %>%
    xtable(type = "latex",
           caption = paste0("Sample size ", 300, "; iteration: ", 100),
           digits = 5)

}

# 600 sample

{
  set.seed(1234)
  record <- list()
  for (iteration in 1:100) {
    support_points <- seq(-5, 6, 1)  # the last support point is pseudo
    B <- 600
    obs_x <- sample(support_points, B, replace = T, prob = target_density)  # simulated data
    actual_data <- obs_x
    actual_data[actual_data == 6] <- 5
    # empirical
    np_pmf <- data.frame(
      x = seq(-5, 5, 1),
      p = sapply(seq(-5, 5, 1), function(x) sum(actual_data == x))/B
      # (cut_interval(obs_x, n = n_bins, right = F) %>% table)/B
    )
    # hal
    # n_bins <- length(support_points) - 1
    n_bins <- diff(range(actual_data)) + 1  # make sure they are integer intervals, and the last one is [5, 6]
    cut_interval(obs_x, n = n_bins, right = F) %>% table
    example_data <- data.frame(W = 0, A = obs_x)
    density_hal <- haldensify(
      A = example_data$A,
      W = example_data$W,
      n_bins = n_bins,
      grid_type = "equal_range",
      lambda_seq = exp(seq(-1, -10, length = 100)), use_future = T
    )
    hal_pmf <- data.frame(x = seq(-5, 5, 1),
                          p = predict(density_hal, new_A = seq(-5, 5, 1), new_W = 0) * (diff(range(obs_x)) / n_bins)
    )
    # biased initial
    temp_p <- predict(density_hal, new_A = seq(-5, 5, 1), new_W = 0)
    temp_p[length(temp_p)] <- min(1, 3*temp_p[length(temp_p)])
    ratio <- (1 - temp_p[length(temp_p)]) / sum(temp_p[-length(temp_p)])
    temp_p[-length(temp_p)] <- temp_p[-length(temp_p)] * ratio
    biased_pmf <- data.frame(x = seq(-5, 5, 1), p = temp_p)
    # Pn update
    updated_pmf <- biased_pmf
    record_epsilon <- c()
    record_D <- c()
    for (step in 1:100) {
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
    for (step in 1:100) {
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

    record[[iteration]] <- c(sum(biased_pmf$p^2),
                             sum(updated_pmf_Pn$p^2),
                             sum(updated_pmf_tildePn$p^2),
                             step_Pn,
                             step_tildePn)
  }
  record_mat <- do.call(rbind, record)
  report <- data.frame(Bias = colMeans(record_mat[,1:3] - truth),
                       SD = record_mat[,1:3] %>% apply(2, sd),
                       MSE = colMeans((record_mat[,1:3] - truth)^2),
                       Step = c(0, colMeans(record_mat[,4:5]))
  )
  rownames(report) <- c("Biased initial", "Pn update", "tilde Pn update")
  report %>%
    xtable(type = "latex",
           caption = paste0("Sample size ", 300, "; iteration: ", 100),
           digits = 5)

}



# 2-TMLE w.r.t. tilde Pn
{
  set.seed(123)
  record <- list()
  for (iteration in 1:100) {
    support_points <- seq(-5, 6, 1)  # the last support point is pseudo
    B <- 300
    obs_x <- sample(support_points, B, replace = T, prob = target_density)  # simulated data
    actual_data <- obs_x
    actual_data[actual_data == 6] <- 5
    # empirical
    np_pmf <- data.frame(
      x = seq(-5, 5, 1),
      p = sapply(seq(-5, 5, 1), function(x) sum(actual_data == x))/B
      # (cut_interval(obs_x, n = n_bins, right = F) %>% table)/B
    )
    # hal
    # n_bins <- length(support_points) - 1
    n_bins <- diff(range(actual_data)) + 1  # make sure they are integer intervals, and the last one is [5, 6]
    cut_interval(obs_x, n = n_bins, right = T) %>% table
    example_data <- data.frame(W = 0, A = obs_x)
    density_hal <- haldensify(
      A = example_data$A,
      W = example_data$W,
      n_bins = n_bins,
      grid_type = "equal_range",
      lambda_seq = exp(seq(-1, -10, length = 100)), use_future = T
    )
    hal_pmf <- data.frame(x = seq(-5, 5, 1),
                          p = predict(density_hal, new_A = seq(-5, 5, 1), new_W = 0) * (diff(range(obs_x)) / n_bins)
    )
    # biased initial
    temp_p <- predict(density_hal, new_A = seq(-5, 5, 1), new_W = 0)
    temp_p[length(temp_p)] <- min(1, 3*temp_p[length(temp_p)])
    ratio <- (1 - temp_p[length(temp_p)]) / sum(temp_p[-length(temp_p)])
    temp_p[-length(temp_p)] <- temp_p[-length(temp_p)] * ratio
    biased_pmf <- data.frame(x = seq(-5, 5, 1), p = temp_p)
    # tilde Pn update
    updated_pmf_1 <- updated_pmf_2 <- updated_pmf <- biased_pmf
    record_epsilon <- record_epsilon_2 <- c()
    record_D_1 <- record_D_2 <- c()
    for (step in 1:100) {
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

    record[[iteration]] <- c(sum(biased_pmf$p^2),
                             # sum(updated_pmf_Pn$p^2),
                             sum(updated_pmf_tildePn_2$p^2),
                             # step_Pn,
                             step_tildePn_2,
                             ratio_2)
  }
  record_mat <- do.call(rbind, record)
  report <- data.frame(Bias = colMeans(record_mat[,1:2] - truth),
                       SD = record_mat[,1:2] %>% apply(2, sd),
                       MSE = colMeans((record_mat[,1:2] - truth)^2),
                       Step = c(0, mean(record_mat[,3])),
                       Ratio = c(0, mean(record_mat[,4]))
  )
  rownames(report) <- c("Biased initial", "tilde Pn 2-TMLE")
  report
  report %>%
    xtable(type = "latex",
           caption = paste0("Sample size ", 300, "; iteration: ", 100),
           digits = 5)

}

# 600 sample
{
  set.seed(1234)
  record <- list()
  for (iteration in 1:100) {
    support_points <- seq(-5, 6, 1)  # the last support point is pseudo
    B <- 600
    obs_x <- sample(support_points, B, replace = T, prob = target_density)  # simulated data
    actual_data <- obs_x
    actual_data[actual_data == 6] <- 5
    # empirical
    np_pmf <- data.frame(
      x = seq(-5, 5, 1),
      p = sapply(seq(-5, 5, 1), function(x) sum(actual_data == x))/B
      # (cut_interval(obs_x, n = n_bins, right = F) %>% table)/B
    )
    # hal
    # n_bins <- length(support_points) - 1
    n_bins <- diff(range(actual_data)) + 1  # make sure they are integer intervals, and the last one is [5, 6]
    cut_interval(obs_x, n = n_bins, right = T) %>% table
    example_data <- data.frame(W = 0, A = obs_x)
    density_hal <- haldensify(
      A = example_data$A,
      W = example_data$W,
      n_bins = n_bins,
      grid_type = "equal_range",
      lambda_seq = exp(seq(-1, -10, length = 100)), use_future = T
    )
    hal_pmf <- data.frame(x = seq(-5, 5, 1),
                          p = predict(density_hal, new_A = seq(-5, 5, 1), new_W = 0) * (diff(range(obs_x)) / n_bins)
    )
    # biased initial
    temp_p <- predict(density_hal, new_A = seq(-5, 5, 1), new_W = 0)
    temp_p[length(temp_p)] <- min(1, 3*temp_p[length(temp_p)])
    ratio <- (1 - temp_p[length(temp_p)]) / sum(temp_p[-length(temp_p)])
    temp_p[-length(temp_p)] <- temp_p[-length(temp_p)] * ratio
    biased_pmf <- data.frame(x = seq(-5, 5, 1), p = temp_p)
    # tilde Pn update
    updated_pmf_1 <- updated_pmf_2 <- updated_pmf <- biased_pmf
    record_epsilon <- record_epsilon_2 <- c()
    record_D_1 <- record_D_2 <- c()
    for (step in 1:100) {
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

    record[[iteration]] <- c(sum(biased_pmf$p^2),
                             # sum(updated_pmf_Pn$p^2),
                             sum(updated_pmf_tildePn_2$p^2),
                             # step_Pn,
                             step_tildePn_2,
                             ratio_2)
  }
  record_mat <- do.call(rbind, record)
  report <- data.frame(Bias = colMeans(record_mat[,1:2] - truth),
                       SD = record_mat[,1:2] %>% apply(2, sd),
                       MSE = colMeans((record_mat[,1:2] - truth)^2),
                       Step = c(0, mean(record_mat[,3])),
                       Ratio = c(0, mean(record_mat[,4]))
  )
  rownames(report) <- c("Biased initial", "tilde Pn 2-TMLE")
  report
  report %>%
    xtable(type = "latex",
           caption = paste0("Sample size ", 300, "; iteration: ", 100),
           digits = 5)

}

# 2-TMLE w.r.t. tilde Pn; initial biased universally
{
  set.seed(123)
  record <- list()
  for (iteration in 1:10) {
    support_points <- seq(-5, 6, 1)  # the last support point is pseudo
    B <- 300
    obs_x <- sample(support_points, B, replace = T, prob = target_density)  # simulated data
    actual_data <- obs_x
    actual_data[actual_data == 6] <- 5
    # empirical
    np_pmf <- data.frame(
      x = seq(-5, 5, 1),
      p = sapply(seq(-5, 5, 1), function(x) sum(actual_data == x))/B
      # (cut_interval(obs_x, n = n_bins, right = F) %>% table)/B
    )
    # hal
    # n_bins <- length(support_points) - 1
    n_bins <- diff(range(actual_data)) + 1  # make sure they are integer intervals, and the last one is [5, 6]
    cut_interval(obs_x, n = n_bins, right = T) %>% table
    example_data <- data.frame(W = 0, A = obs_x)
    density_hal <- haldensify(
      A = example_data$A,
      W = example_data$W,
      n_bins = n_bins,
      grid_type = "equal_range",
      lambda_seq = exp(seq(-1, -10, length = 100)), use_future = T
    )
    hal_pmf <- data.frame(x = seq(-5, 5, 1),
                          p = predict(density_hal, new_A = seq(-5, 5, 1), new_W = 0) * (diff(range(obs_x)) / n_bins)
    )
    # biased initial
    temp_p <- hal_pmf$p
    temp_p <- temp_p + 0.2
    temp_p <- temp_p/sum(temp_p)
    biased_pmf <- data.frame(x = seq(-5, 5, 1), p = temp_p)
    # tilde Pn update
    updated_pmf_1 <- updated_pmf_2 <- updated_pmf <- biased_pmf
    record_epsilon <- record_epsilon_2 <- c()
    record_D_1 <- record_D_2 <- c()
    for (step in 1:100) {
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

    record[[iteration]] <- c(sum(biased_pmf$p^2),
                             # sum(updated_pmf_Pn$p^2),
                             sum(updated_pmf_tildePn_2$p^2),
                             # step_Pn,
                             step_tildePn_2,
                             ratio_2)
  }
  record_mat <- do.call(rbind, record)
  report <- data.frame(Bias = colMeans(record_mat[,1:2] - truth),
                       SD = record_mat[,1:2] %>% apply(2, sd),
                       MSE = colMeans((record_mat[,1:2] - truth)^2),
                       Step = c(0, mean(record_mat[,3])),
                       Ratio = c(0, mean(record_mat[,4]))
  )
  rownames(report) <- c("Biased initial", "tilde Pn 2-TMLE")
  report
  report %>%
    xtable(type = "latex",
           caption = paste0("Sample size ", 300, "; iteration: ", 100),
           digits = 5)

}
# 1-TMLE w.r.t. Pn and tilde Pn; initial biased universally
{
  set.seed(123)
  record <- list()
  for (iteration in 1:100) {
    support_points <- seq(-5, 6, 1)  # the last support point is pseudo
    B <- 300
    obs_x <- sample(support_points, B, replace = T, prob = target_density)  # simulated data
    actual_data <- obs_x
    actual_data[actual_data == 6] <- 5
    # empirical
    np_pmf <- data.frame(
      x = seq(-5, 5, 1),
      p = sapply(seq(-5, 5, 1), function(x) sum(actual_data == x))/B
      # (cut_interval(obs_x, n = n_bins, right = F) %>% table)/B
    )
    # hal
    # n_bins <- length(support_points) - 1
    n_bins <- diff(range(actual_data)) + 1  # make sure they are integer intervals, and the last one is [5, 6]
    cut_interval(obs_x, n = n_bins, right = F) %>% table
    example_data <- data.frame(W = 0, A = obs_x)
    density_hal <- haldensify(
      A = example_data$A,
      W = example_data$W,
      n_bins = n_bins,
      grid_type = "equal_range",
      lambda_seq = exp(seq(-1, -10, length = 100)), use_future = T
    )
    hal_pmf <- data.frame(x = seq(-5, 5, 1),
                          p = predict(density_hal, new_A = seq(-5, 5, 1), new_W = 0) * (diff(range(obs_x)) / n_bins)
    )
    # biased initial
    temp_p <- hal_pmf$p
    temp_p <- temp_p + 0.2
    temp_p <- temp_p/sum(temp_p)
    biased_pmf <- data.frame(x = seq(-5, 5, 1), p = temp_p)
    # Pn update
    updated_pmf <- biased_pmf
    record_epsilon <- c()
    record_D <- c()
    for (step in 1:100) {
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
    for (step in 1:100) {
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

    record[[iteration]] <- c(sum(biased_pmf$p^2),
                             sum(updated_pmf_Pn$p^2),
                             sum(updated_pmf_tildePn$p^2),
                             step_Pn,
                             step_tildePn)
  }
  record_mat <- do.call(rbind, record)
  report <- data.frame(Bias = colMeans(record_mat[,1:3] - truth),
                       SD = record_mat[,1:3] %>% apply(2, sd),
                       MSE = colMeans((record_mat[,1:3] - truth)^2),
                       Step = c(0, colMeans(record_mat[,4:5]))
  )
  rownames(report) <- c("Biased initial", "Pn update", "tilde Pn update")
  report %>%
    xtable(type = "latex",
           caption = paste0("Sample size ", 300, "; iteration: ", 100),
           digits = 5)

}

# n=600 with universal bias
{
  # 2-TMLE w.r.t. tilde Pn; initial biased universally
  {
    set.seed(123)
    record <- list()
    for (iteration in 1:10) {
      support_points <- seq(-5, 6, 1)  # the last support point is pseudo
      B <- 300
      obs_x <- sample(support_points, B, replace = T, prob = target_density)  # simulated data
      actual_data <- obs_x
      actual_data[actual_data == 6] <- 5
      # empirical
      np_pmf <- data.frame(
        x = seq(-5, 5, 1),
        p = sapply(seq(-5, 5, 1), function(x) sum(actual_data == x))/B
        # (cut_interval(obs_x, n = n_bins, right = F) %>% table)/B
      )
      # hal
      # n_bins <- length(support_points) - 1
      n_bins <- diff(range(actual_data)) + 1  # make sure they are integer intervals, and the last one is [5, 6]
      cut_interval(obs_x, n = n_bins, right = T) %>% table
      example_data <- data.frame(W = 0, A = obs_x)
      density_hal <- haldensify(
        A = example_data$A,
        W = example_data$W,
        n_bins = n_bins,
        grid_type = "equal_range",
        lambda_seq = exp(seq(-1, -10, length = 100)), use_future = T
      )
      hal_pmf <- data.frame(x = seq(-5, 5, 1),
                            p = predict(density_hal, new_A = seq(-5, 5, 1), new_W = 0) * (diff(range(obs_x)) / n_bins)
      )
      # biased initial
      temp_p <- predict(density_hal, new_A = seq(-5, 5, 1), new_W = 0)
      temp_p <- temp_p + 0.2
      temp_p <- temp_p/sum(temp_p)
      biased_pmf <- data.frame(x = seq(-5, 5, 1), p = temp_p)
      # tilde Pn update
      updated_pmf_1 <- updated_pmf_2 <- updated_pmf <- biased_pmf
      record_epsilon <- record_epsilon_2 <- c()
      record_D_1 <- record_D_2 <- c()
      for (step in 1:100) {
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

      record[[iteration]] <- c(sum(biased_pmf$p^2),
                               # sum(updated_pmf_Pn$p^2),
                               sum(updated_pmf_tildePn_2$p^2),
                               # step_Pn,
                               step_tildePn_2,
                               ratio_2)
    }
    record_mat <- do.call(rbind, record)
    report <- data.frame(Bias = colMeans(record_mat[,1:2] - truth),
                         SD = record_mat[,1:2] %>% apply(2, sd),
                         MSE = colMeans((record_mat[,1:2] - truth)^2),
                         Step = c(0, mean(record_mat[,3])),
                         Ratio = c(0, mean(record_mat[,4]))
    )
    rownames(report) <- c("Biased initial", "tilde Pn 2-TMLE")
    report
    report %>%
      xtable(type = "latex",
             caption = paste0("Sample size ", 300, "; iteration: ", 100),
             digits = 5)

  }
  # 1-TMLE w.r.t. Pn and tilde Pn; initial biased universally
  {
    set.seed(123)
    record <- list()
    for (iteration in 1:100) {
      support_points <- seq(-5, 6, 1)  # the last support point is pseudo
      B <- 300
      obs_x <- sample(support_points, B, replace = T, prob = target_density)  # simulated data
      actual_data <- obs_x
      actual_data[actual_data == 6] <- 5
      # empirical
      np_pmf <- data.frame(
        x = seq(-5, 5, 1),
        p = sapply(seq(-5, 5, 1), function(x) sum(actual_data == x))/B
        # (cut_interval(obs_x, n = n_bins, right = F) %>% table)/B
      )
      # hal
      # n_bins <- length(support_points) - 1
      n_bins <- diff(range(actual_data)) + 1  # make sure they are integer intervals, and the last one is [5, 6]
      cut_interval(obs_x, n = n_bins, right = F) %>% table
      example_data <- data.frame(W = 0, A = obs_x)
      density_hal <- haldensify(
        A = example_data$A,
        W = example_data$W,
        n_bins = n_bins,
        grid_type = "equal_range",
        lambda_seq = exp(seq(-1, -10, length = 100)), use_future = T
      )
      hal_pmf <- data.frame(x = seq(-5, 5, 1),
                            p = predict(density_hal, new_A = seq(-5, 5, 1), new_W = 0) * (diff(range(obs_x)) / n_bins)
      )
      # biased initial
      temp_p <- predict(density_hal, new_A = seq(-5, 5, 1), new_W = 0)
      temp_p <- temp_p + 0.2
      temp_p <- temp_p/sum(temp_p)
      biased_pmf <- data.frame(x = seq(-5, 5, 1), p = temp_p)
      # Pn update
      updated_pmf <- biased_pmf
      record_epsilon <- c()
      record_D <- c()
      for (step in 1:100) {
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
      for (step in 1:100) {
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

      record[[iteration]] <- c(sum(biased_pmf$p^2),
                               sum(updated_pmf_Pn$p^2),
                               sum(updated_pmf_tildePn$p^2),
                               step_Pn,
                               step_tildePn)
    }
    record_mat <- do.call(rbind, record)
    report <- data.frame(Bias = colMeans(record_mat[,1:3] - truth),
                         SD = record_mat[,1:3] %>% apply(2, sd),
                         MSE = colMeans((record_mat[,1:3] - truth)^2),
                         Step = c(0, colMeans(record_mat[,4:5]))
    )
    rownames(report) <- c("Biased initial", "Pn update", "tilde Pn update")
    report %>%
      xtable(type = "latex",
             caption = paste0("Sample size ", 300, "; iteration: ", 100),
             digits = 5)

  }

}

# predictions to recover conditional density of A|W
new_a <- seq(-7, 7, by = 0.01)
new_dat <- as.data.table(list(a = new_a,
                              w_zero = rep(0, length(new_a))
))
new_dat[, pred_w_zero := predict(density_hal, new_A = new_dat$a,
                                 new_W = new_dat$w_zero)]
# visualize results
dens_dat <-  melt(new_dat, id = c("a"),
                  measure.vars = c("pred_w_zero"))
p_dens <- ggplot(dens_dat, aes(x = a, y = value, colour = variable)) +
  geom_point() +
  geom_line() +
  stat_function(fun = data_out$f_d, args = list(),
                colour = "darkgreen", linetype = "dashed") +
  xlab("Observed value") +
  ylab("Predicted probability") +
  ggtitle("Conditional density p(A|W)") +
  theme_bw() +
  theme(legend.position = "none")
p_dens



support_points <- seq(-7, 7, 0.01)
new_data <- data.table(W = 0,
                       A = support_points)
new_dens <- data.frame(outcome = new_data$A,
                       cont_dens = predict(density_hal, new_A = new_data$A, new_W = new_data$W))

temp <- c(1)
for (i in 2:nrow(new_data)) {
  if (new_dens$cont_dens[i] != new_dens$cont_dens[i-1]) temp <- c(temp, last(temp)+1) else temp <- c(temp, last(temp))
}
new_dens[["group"]] <- temp
new_dens$cont_dens %>% sum * 0.01

new_pmf <- sapply(new_dens$group %>% unique, function(which_group) {
  temp <- which(new_dens$group == which_group)
  a1 <- first(temp)
  a2 <- last(temp)
  c(new_dens$outcome[a1],
    (new_dens$outcome[a2] - new_dens$outcome[a1]) *
      new_dens$cont_dens[a1],
    new_dens$cont_dens[a1]
  )
}) %>% t
new_pmf
new_pmf[, 2] %>% sum
predict(density_hal, new_A = seq(-5, 5, 1), new_W = 0) %>% sum
(cut_interval(obs_x, n = n_bins, right = F) %>% table)/B



sum(predict(density_hal, new_A = seq(-5, 5, 1), new_W = 0)^2)




B <- 800
n_mode <- 2
bin_width <- 1e-2

record <- c()
for (i in 1:100) {
  data_obs <- simulate_data(n_sim = B, n_mode = n_mode)
  example_data <- as.data.table(list(W = 0, A = data_obs$x))
  density_hal <- haldensify(
    A = example_data$A,
    W = example_data$W,
    n_bins = c(10),
    grid_type = "equal_mass",
    lambda_seq = exp(seq(1, -10, length = 100)), use_future = T
  )

  new_a <- seq(-10, 10, by = bin_width)
  new_dat <- as.data.table(list(a = new_a,
                                w_zero = rep(0, length(new_a))
  ))
  record[i] <- sum(predict(density_hal, new_A = new_dat$a, new_W = new_dat$w_zero)^2 * bin_width)
}
mean(record) - truth_sim
record %>% sd
mean((record - truth_sim)^2)



B <- 800
n_mode <- 2
bin_width <- 1e-2

data_obs <- simulate_data(n_sim = B, n_mode = n_mode)
example_data <- as.data.table(list(W = 0, A = data_obs$x))
density_hal <- haldensify(
  A = example_data$A,
  W = example_data$W,
  n_bins = c(10),
  grid_type = "equal_range",
  cv_folds = 10,
  lambda_seq = exp(seq(1, -10, length = 100)), use_future = T
)



# predictions to recover conditional density of A|W
new_a <- seq(-10, 10, by = bin_width)
new_dat <- as.data.table(list(a = new_a,
                              w_zero = rep(0, length(new_a))
))
new_dat[, pred_w_zero := predict(density_hal, new_A = new_dat$a,
                                 new_W = new_dat$w_zero)]

# visualize results
dens_dat <-  melt(new_dat, id = c("a"),
                  measure.vars = c("pred_w_zero"))
p_dens <- ggplot(dens_dat, aes(x = a, y = value, colour = variable)) +
  geom_point() +
  geom_line() +
  stat_function(fun = data_obs$f_d, args = list(),
                colour = "darkgreen", linetype = "dashed") +
  xlab("Observed value") +
  ylab("Predicted probability") +
  ggtitle("Conditional density p(A|W)") +
  theme_bw() +
  theme(legend.position = "none")
p_dens

