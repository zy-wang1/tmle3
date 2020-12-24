library(dplyr)
library(haldensify)  # haldensify
library(ggplot2)  # ggplot

# example from Cai's paper
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

n_sim <- 1e6
n_mode <- 3
set.seed(1234)
data_out <- simulate_data(n_sim = n_sim, n_mode = n_mode)
data_out$x %>% hist  # mixture of n_mode=3 Gaussian densities





###
# issue 1, irregular high tail
###

{
  B <- 800  # sample size
  n_mode <- 2  # number of peaks of the Gaussian mixture
  bin_width <- 1e-2

  set.seed(123)
  data_obs <- simulate_data(n_sim = B, n_mode = n_mode)
  example_data <- as.data.table(list(W = 0, A = data_obs$x))
  density_hal <- haldensify(
    A = example_data$A,
    W = example_data$W,
    n_bins = c(10),
    grid_type = "equal_range",
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
}


###
# issue 2, fail to converge with large b_bins
###


{
  # 2.1, equal_range, n_bins=15, fail to converge
  B <- 800
  n_mode <- 3
  bin_width <- 1e-2

  set.seed(123)
  data_obs <- simulate_data(n_sim = B, n_mode = n_mode)
  data_obs$x %>% density %>% plot
  example_data <- as.data.table(list(W = 0, A = data_obs$x))
  density_hal <- haldensify(
    A = example_data$A,
    W = example_data$W,
    n_bins = c(15),
    grid_type = "equal_range",
    lambda_seq = exp(seq(1, -10, length = 100)),
    use_future = T
  )
}

{
  # 2.2, it converges with equal_mass (but still with the high tail)
  B <- 800
  n_mode <- 3
  bin_width <- 1e-2

  set.seed(123)
  data_obs <- simulate_data(n_sim = B, n_mode = n_mode)
  data_obs$x %>% density %>% plot
  example_data <- as.data.table(list(W = 0, A = data_obs$x))
  density_hal <- haldensify(
    A = example_data$A,
    W = example_data$W,
    n_bins = c(15),
    grid_type = "equal_mass",
    lambda_seq = exp(seq(1, -10, length = 100)),
    use_future = T
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
}


{
  # 2.3, equal_mass also fails with larger n_bins=20
  B <- 800
  n_mode <- 3
  bin_width <- 1e-2

  set.seed(123)
  data_obs <- simulate_data(n_sim = B, n_mode = n_mode)
  data_obs$x %>% density %>% plot
  example_data <- as.data.table(list(W = 0, A = data_obs$x))
  density_hal <- haldensify(
    A = example_data$A,
    W = example_data$W,
    n_bins = c(20),
    grid_type = "equal_mass",
    lambda_seq = exp(seq(1, -10, length = 100)),
    use_future = T
  )
}
