library(data.table)  # setDT
library(dplyr)  # dplyr::select, to mask other packages
library(ggplot2)  # ggplot
library(xtable)

# library(haldensify)  # haldensify
source("./temp_code/haldensify.R")

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



# discretized density: see how truth changes as numbers of peaks and supports increase
record <- matrix(0, 10, 10)
for (i in 1:10) for (j in 1:10) {
  n_mode <- i
  n_actual_bins <-  j*10

  {

    # generate true discrete density
    n_sim <- 1e6
    set.seed(1234)
    data_out <- simulate_data(n_sim = n_sim, n_mode = n_mode)
    # truth <- mean(data_out$f_d(data_out$x))
    # truth
    # data_out$x %>% hist(freq = F, ylim = c(0, 0.5))
    data_out$x %>% density %>% plot(ylim = c(0, 0.8))
    cdf_obj <- data_out$x %>% ecdf
    bin_width <- 20/n_actual_bins
    support_points <- seq(-10, 10, by = bin_width)
    target_density <- sapply(1:length(support_points), function(i) {
      if (i == 1) cdf_obj(support_points[i]) else cdf_obj(support_points[i]) - cdf_obj(support_points[i - 1])
    })
    target_density[length(target_density)] <- 1 - sum(target_density[1:(length(target_density) - 1)])
    true_pmf <- data.frame(x = support_points,
                           p = target_density)
    true_pdf <- data.frame(x = support_points,
                           p = target_density/bin_width)
  }
  record[i, j] <- sum(true_pdf$p^2 * bin_width)
}
rownames(record) <- paste(1:10, "peaks")
colnames(record) <- paste((1:10)*10, "supports")
record %>% round(3)

truth_vec <- rep(0, 10)
for (n_mode in 1:10) {
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
  truth <- mean(data_out$f_d(data_out$x))
  truth_vec[n_mode] <- truth
}
record <- cbind(record, truth_vec)
colnames(record)[ncol(record)] <- "continuous truth"

record[,] %>% round(4)
record[, c((1:5)*2, 11)] %>% round(4) %>%
  xtable(type = "latex",
         caption = paste0("Sample size ", 300, "; iteration: ", 100),
         digits = 5)




# discretization; see how TV changes
record_TV <- matrix(0, 10, 10)
true_TV <- rep(0, 10)
for (i in 1:10) for (j in 1:10) {
  n_mode <- i
  n_actual_bins <-  j*10

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
    temp_grid <- seq(-10, 10, 0.001)
    temp_f <- c(0, data_out$f_d(temp_grid), 0)
    temp_diff <- sapply(1:(length(temp_f)-1), function(temp_loc) abs(temp_f[temp_loc] - temp_f[temp_loc+1]))
    true_TV[i] <- sum(temp_diff)

    data_out$x %>% density %>% plot(ylim = c(0, 0.8))

    cdf_obj <- data_out$x %>% ecdf
    bin_width <- 20/n_actual_bins
    support_points <- seq(-10, 10, by = bin_width)
    target_density <- sapply(1:length(support_points), function(i) {
      if (i == 1) cdf_obj(support_points[i]) else cdf_obj(support_points[i]) - cdf_obj(support_points[i - 1])
    })
    target_density[length(target_density)] <- 1 - sum(target_density[1:(length(target_density) - 1)])
    true_pmf <- data.frame(x = support_points,
                           p = target_density)
    true_pdf <- data.frame(x = support_points,
                           p = target_density/bin_width)
    temp_p <- c(0, true_pdf$p, 0)
    temp_diff <- sapply(1:(length(temp_p)-1), function(temp_loc) abs(temp_p[temp_loc] - temp_p[temp_loc+1]))
    # temp_p <- c(0, true_pmf$p, 0)
    # temp_diff <- sapply(1:(length(temp_p)-1), function(temp_loc) abs(temp_p[temp_loc] - temp_p[temp_loc+1]))

    record_TV[i, j] <- sum(temp_diff)
  }
}
rownames(record_TV) <- paste(1:10, "peaks")
colnames(record_TV) <- paste((1:10)*10, "supports")
record_TV <- cbind(record_TV, true_TV)

record_TV[, c((1:5)*2, 11)] %>% round(3)%>%
  xtable(type = "latex",
         caption = paste0("Sample size ", 300, "; iteration: ", 100),
         digits = 5)




# truth [-5, 5]
# discretized density: see how truth changes as numbers of peaks and supports increase
record <- matrix(0, 10, 10)
for (i in 1:10) for (j in 1:10) {
  n_mode <- i
  n_actual_bins <-  j*10

  {

    # generate true discrete density
    n_sim <- 1e6
    set.seed(1234)
    data_out <- simulate_data(n_sim = n_sim, n_mode = n_mode)
    # truth <- mean(data_out$f_d(data_out$x))
    # truth
    # data_out$x %>% hist(freq = F, ylim = c(0, 0.5))
    data_out$x %>% density %>% plot(ylim = c(0, 0.8))
    cdf_obj <- data_out$x %>% ecdf
    bin_width <- 10/n_actual_bins
    support_points <- seq(-5, 5, by = bin_width)
    target_density <- sapply(1:length(support_points), function(i) {
      if (i == 1) cdf_obj(support_points[i]) else cdf_obj(support_points[i]) - cdf_obj(support_points[i - 1])
    })
    target_density[length(target_density)] <- 1 - sum(target_density[1:(length(target_density) - 1)])
    true_pmf <- data.frame(x = support_points,
                           p = target_density)
    true_pdf <- data.frame(x = support_points,
                           p = target_density/bin_width)
  }
  record[i, j] <- sum(true_pdf$p^2 * bin_width)
}
rownames(record) <- paste(1:10, "peaks")
colnames(record) <- paste((1:10)*10, "supports")
record %>% round(3)

truth_vec <- rep(0, 10)
for (n_mode in 1:10) {

  # generate true discrete density
  n_sim <- 1e6
  set.seed(1234)
  data_out <- simulate_data(n_sim = n_sim, n_mode = n_mode)
  truth <- mean(data_out$f_d(data_out$x))
  truth_vec[n_mode] <- truth
}
record <- cbind(record, truth_vec)
colnames(record)[ncol(record)] <- "continuous truth"

record[4:10,] %>% round(4)
record[4:10, c((1:5)*2, 11)] %>% round(4) %>%
  xtable(type = "latex",
         caption = paste0("Sample size ", 300, "; iteration: ", 100),
         digits = 5)


# TV [-5, 5]
# discretization; see how TV changes
record_TV <- matrix(0, 10, 10)
true_TV <- rep(0, 10)
for (i in 1:10) for (j in 1:10) {
  n_mode <- i
  n_actual_bins <-  j*10

  {
    # generate true discrete density
    n_sim <- 1e6
    set.seed(1234)
    data_out <- simulate_data(n_sim = n_sim, n_mode = n_mode)
    temp_grid <- seq(-10, 10, 0.001)
    temp_f <- c(0, data_out$f_d(temp_grid), 0)
    temp_diff <- sapply(1:(length(temp_f)-1), function(temp_loc) abs(temp_f[temp_loc] - temp_f[temp_loc+1]))
    true_TV[i] <- sum(temp_diff)

    data_out$x %>% density %>% plot(ylim = c(0, 0.8))

    cdf_obj <- data_out$x %>% ecdf
    bin_width <- 10/n_actual_bins
    support_points <- seq(-5, 5, by = bin_width)
    target_density <- sapply(1:length(support_points), function(i) {
      if (i == 1) cdf_obj(support_points[i]) else cdf_obj(support_points[i]) - cdf_obj(support_points[i - 1])
    })
    target_density[length(target_density)] <- 1 - sum(target_density[1:(length(target_density) - 1)])
    true_pmf <- data.frame(x = support_points,
                           p = target_density)
    true_pdf <- data.frame(x = support_points,
                           p = target_density/bin_width)
    temp_p <- c(0, true_pdf$p, 0)
    temp_diff <- sapply(1:(length(temp_p)-1), function(temp_loc) abs(temp_p[temp_loc] - temp_p[temp_loc+1]))
    # temp_p <- c(0, true_pmf$p, 0)
    # temp_diff <- sapply(1:(length(temp_p)-1), function(temp_loc) abs(temp_p[temp_loc] - temp_p[temp_loc+1]))

    record_TV[i, j] <- sum(temp_diff)
  }
}
rownames(record_TV) <- paste(1:10, "peaks")
colnames(record_TV) <- paste((1:10)*10, "supports")
record_TV <- cbind(record_TV, true_TV)

record_TV[4:10, c((1:5)*2, 11)] %>% round(3)%>%
  xtable(type = "latex",
         caption = paste0("Sample size ", 300, "; iteration: ", 100),
         digits = 5)




# discretization; see second order remainder

record_R <- matrix(0, 10, 10)
for (i in 1:10) for (j in 1:10) {
  n_mode <- i
  n_actual_bins <-  j*10

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
    temp_grid <- seq(-10, 10, 0.001)
    temp_f <- data_out$f_d(temp_grid)

    data_out$x %>% density %>% plot(ylim = c(0, 0.8))

    cdf_obj <- data_out$x %>% ecdf
    bin_width <- 20/n_actual_bins
    support_points <- seq(-10, 10, by = bin_width)
    target_density <- sapply(1:length(support_points), function(i) {
      if (i == 1) cdf_obj(support_points[i]) else cdf_obj(support_points[i]) - cdf_obj(support_points[i - 1])
    })
    target_density[length(target_density)] <- 1 - sum(target_density[1:(length(target_density) - 1)])
    true_pmf <- data.frame(x = support_points,
                           p = target_density)
    true_pdf <- data.frame(x = support_points,
                           p = target_density/bin_width)
    temp_p <- sapply(1:length(temp_grid), function(temp_loc) {
      true_pdf$p[first(which(true_pdf$x >= temp_grid[temp_loc]))]
    })
    temp_diff <- -(temp_f - temp_p)^2 * 0.001
    record_R[i, j] <- sum(temp_diff)
  }
}
rownames(record_R) <- paste(1:10, "peaks")
colnames(record_R) <- paste((1:10)*10, "supports")
(-record_R) %>% round(3)





# 4-10 peaks, [-5, 5] 30 grids; biased, HAL, empirical, 1-TMLE, 2-TMLE comparison
{
  B <- 300
library(parallel)
  set.seed(123)

  report_list <- mclapply(3:10, mc.cores = 8, function(n_mode) {
    # generate true discrete density
    n_sim <- 1e6
    # set.seed(1234)
    data_out <- simulate_data(n_sim = n_sim, n_mode = n_mode)
    # truth <- mean(data_out$f_d(data_out$x))
    # truth
    # data_out$x %>% hist(freq = F, ylim = c(0, 0.5))
    data_out$x %>% density %>% plot(ylim = c(0, 0.8))
    cdf_obj <- data_out$x %>% ecdf
    bin_width <- 0.25
    support_points <- seq(-5, 5, by = bin_width)
    target_density <- sapply(1:length(support_points), function(i) {
      if (i == 1) cdf_obj(support_points[i]) else cdf_obj(support_points[i]) - cdf_obj(support_points[i - 1])
    })
    target_density[length(target_density)] <- 1 - sum(target_density[1:(length(target_density) - 1)])
    true_pmf <- data.frame(x = support_points,
                           p = target_density)
    true_pdf <- data.frame(x = support_points,
                           p = target_density/bin_width)


    n_iteration <- 72
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
      left_breaks <- density_hal$breaks[-length(density_hal$breaks)]
      support_points_used <- support_points[support_points <= range(obs_x)[2] & support_points >= range(obs_x)[1]]
      # predict(density_hal, new_A = support_points_used, new_W = 0)
      # predict(density_hal, new_A = left_breaks, new_W = 0)
      hal_pmf <- data.frame(x = support_points, p = 0)
      hal_pmf[support_points %in% support_points_used, "p"] <- predict(density_hal, new_A = support_points_used, new_W = 0) * density_hal$bin_sizes
      # sum(hal_pmf$p)
      hal_pdf <- data.frame(x = hal_pmf$x,
                            p = hal_pmf$p / bin_width)

      # biased initial
      temp_p <- hal_pmf$p
      temp_p <- temp_p + 0.1
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

}

# this is saved for different modes
# report_list %>% saveRDS("./temp_output/higher_order_ZW-20201223.rds")


report %>%
  xtable(type = "latex",
         caption = paste0("Sample size ", 300, "; iteration: ", 100),
         digits = 5)






plot_df <- report_list %>% lapply(function(x) x[, 1]) %>% abind::abind(along = 0)
colnames(plot_df) <- c("NP-MLE", "HAL-MLE", "Biased", "Pn 1TMLE", "tilde 1TMLE", "tilde 2TMLE")
plot_df <- data.frame(plot_df)
plot_df[["TV"]] <- record_TV[3:10, 4]
plot_df[["peaks"]] <- 3:10

library(reshape2)
library(ggrepel)
test_data_long <- reshape2::melt(plot_df[, -7], id="peaks")  # convert to long format

ggplot(data=test_data_long,
       aes(x=peaks, y=value, colour=variable, label = variable)) +
  geom_line() + geom_text_repel(data = test_data_long %>% filter(peaks == 3)) + geom_point() + theme_classic() + geom_hline(yintercept = 0, linetype="dashed", color = "black") +
  labs(y = "Bias", x = "Number of Peaks")




plot_df <- report_list %>% lapply(function(x) x[, 3]) %>% abind::abind(along = 0)
colnames(plot_df) <- c("NP-MLE", "HAL-MLE", "Biased", "Pn 1TMLE", "tilde 1TMLE", "tilde 2TMLE")
plot_df <- data.frame(plot_df)
plot_df[["TV"]] <- record_TV[3:10, 4]
plot_df[["peaks"]] <- 3:10
plot_df


test_data_long <- reshape2::melt(plot_df[, -7], id="peaks")  # convert to long format

ggplot(data=test_data_long,
       aes(x=peaks, y=value, colour=variable, label = variable)) +
  geom_line() + geom_text_repel(data = test_data_long %>% filter(peaks == 3)) + geom_point() + theme_classic() + geom_hline(yintercept = 0, linetype="dashed", color = "black") +
  labs(y = "MSE", x = "Number of Peaks")





