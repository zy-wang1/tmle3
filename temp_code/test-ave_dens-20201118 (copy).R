library(R6)  # R6class
library(data.table)  # setDT
library(sl3)
library(digest)
library(uuid)  # UUIDgenerate
library(delayed)  # bundle_delayed
library(assertthat)  # assert_that
library(speedglm)  # speedglm
# library(methods)  # is

library(dplyr)  # dplyr::select, to mask other packages

code_list <- list.files("./R", full.names = T)
for (code in code_list) source(code)
source("./temp_code/generate_data.R")

###
# generate target density by discretizing Cai's example
###

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
n_mode <- 10
# bin_width <- 5e-1
set.seed(1234)
data_out <- simulate_data(n_sim = n_sim, n_mode = n_mode)
truth <- mean(data_out$f_d(data_out$x))
truth

data_out$x %>% hist
cdf_obj <- data_out$x %>% ecdf

support_points <- seq(-20, 20, length.out = 70)
target_density <- sapply(1:length(support_points), function(i) {
  if (i == 1) cdf_obj(support_points[i]) else cdf_obj(support_points[i]) - cdf_obj(support_points[i - 1])
})
target_density[length(target_density)] <- 1 - sum(target_density[1:(length(target_density) - 1)])

# support_points
# target_density %>% sum
# truth <- sum(target_density^2)
# truth
true_dens <- data.frame(x = support_points, p = target_density)
true_dens %>% lines

###
# simulate B observations
###

B <- 1000

# set.seed(1234)
data_obs <- simulate_data(n_sim = B, n_mode = n_mode)
obs_x <- data_obs$x
# obs_x <- sample(support_points, B, replace = T, prob = target_density)  # simulated data
one <- rep(0, B)
# obs_x <- as.factor(obs_x)
# one <- as.factor(one)
obs_data <- data.frame(
  one = one,
  outcome = obs_x)
# sum((table(obs_x) / B)^2)  # empirical Pn plug-in; for categorical
truth

###
# generate tmle tasks and initial likelihoods
###

node_list <- lapply(names(obs_data), function(x) x)
names(node_list) <- names(obs_data)
variable_types <- NULL
# variable_types <- list("one" = "categorical",
#                        "x" = "categorical")
npsem <- lapply(1:length(node_list), function(k) {
  if (k > 1) {
    define_node(names(node_list)[k],
                node_list[[k]],
                names(node_list)[1:(k-1)],
                variable_type = variable_types[[ names(node_list)[k] ]])
  } else {
    define_node(name = names(node_list)[k],
                variables = node_list[[ names(node_list)[k] ]],
                variable_type = variable_types[[ names(node_list)[k] ]])
  }
})
names(npsem) <- names(node_list)

# setDT(obs_data)
tmle_task <- tmle3_Task$new(obs_data, npsem = npsem)

# lrnr_mean <- make_learner(Lrnr_mean)
# glm_learner <- Lrnr_glm$new()
# doesn't work
{
  hal <- Lrnr_hal9001$new(lambda_seq = exp(seq(-1, -13, length = 100)))
  lrnr_pooled <- Lrnr_pooled_hazards$new(hal)
  density_learner <- Lrnr_density_discretize$new(lrnr_pooled,
                                                 type = "equal_mass",
                                                 n_bins = 10
  )
  learner_list <- lapply(1:length(npsem), function(i) {
    density_learner
  })
}

# simple version
{
  hal_dens <- Lrnr_haldensify$new(
    grid_type = "equal_range",
    n_bins = 11,
    lambda_seq = exp(seq(-1, -10, length = 100))
  )
  sl_learner_density <- Lrnr_sl$new(learners = hal_dens
                                    # , metalearner = Lrnr_solnp_density$new()
                                    )
}

{
  lrnr_ranger50 <- make_learner(Lrnr_ranger, num.trees = 50)
  learner_list <- lapply(1:length(tmle_task$npsem), function(s) Lrnr_sl$new(
    learners = list(
      # lrnr_mean,
      # lrnr_glm,
      # lrnr_glm_fast
      # ,
      lrnr_ranger50
      # ,
      # lrnr_hal_simple
      # ,
      # lrnr_lasso, lrnr_ridge, lrnr_elasticnet
    )
  ))
}

# haldensify not working for now
{
  # SL learners for conditional densities to be used for the propensity score fit
  lrnr_dens_lib <- make_learner_stack(
    # list("Lrnr_haldensify", n_bins = 15,
    #                                        grid_type = "equal_range",
    #                                        lambda_seq = exp(seq(-1, -10,
    # length = 100)))
    # ,
    list("Lrnr_haldensify", n_bins = 10,
         bin_method = "equal_mass",
         lambda_seq = exp(seq(-1, -10,
                              length = 100)))
  )
  sl_learner_density <- Lrnr_sl$new(learners = lrnr_dens_lib,
                                    metalearner = Lrnr_solnp_density$new())
}

lrn_1 <- Lrnr_condensier$new(
  nbins = 25,
  bin_method = "equal.len",
  pool = FALSE,
  bin_estimator = condensier::speedglmR6$new()
)
sl_learner_density <- Lrnr_sl$new(learners = lrn_1,
                                  metalearner = Lrnr_solnp_density$new())

learner_list <- c(0, sl_learner_density)
names(learner_list) <- names(npsem)

# make likelihood

factor_list <- list()
temp_names <- names(npsem)
factor_list <- lapply(1:length(npsem), function(k) {
  if (k == 1)
    define_lf(LF_emp, names(npsem)[k])
  else
    # define_lf(LF_emp, names(npsem)[k])
    LF_fit$new(temp_names[k], learner = learner_list[[ temp_names[k] ]], type = "density")
})
names(factor_list) <- names(npsem)

likelihood_def <- Likelihood$new(factor_list)
initial_likelihood <- likelihood_def$train(tmle_task)

# data.frame(obs_x, initial_likelihood$get_likelihood(tmle_task, "outcome")) %>% plot
# tmle_task$data$outcome %>% density %>% lines
#
# data.frame(support_points, target_density) %>% points(col = "red")
#
# unique(initial_likelihood$get_likelihood(new_task, "outcome") )^2 %>% sum
# truth
# unique(initial_likelihood$get_likelihood(new_task, "outcome") ) %>% sum



###
# design param and non-targeting est
###

tmle_params <- list(
  Param_ave_dens$new(observed_likelihood = initial_likelihood)
)

suppressMessages(
  nontargeting <- tmle_params[[1]]$estimates(tmle_task)
)
temp_IC <- nontargeting$IC
var_D <- var(temp_IC)
n <- length(temp_IC)
se <- sqrt(var_D / n)
CI2 <- nontargeting$psi + 1.96 * se
CI1 <- nontargeting$psi - 1.96 * se
nontargeting$psi


obs_x %>% density %>% plot(ylim = c(0, 1))
{
  self <- tmle_params[[1]]
  support_points <- seq(-10, 10, by = 0.1)
  new_data <- data.table(one = 0,
                         outcome = support_points)
  new_task <- tmle3_Task$new(new_data, npsem)
  new_dens <- data.frame(outcome = new_data$outcome,
                         cont_dens = self$observed_likelihood$get_likelihood(new_task, "outcome", fold_number))

  temp <- c(1)
  for (i in 2:nrow(new_data)) {
    if (new_dens$cont_dens[i] != new_dens$cont_dens[i-1]) temp <- c(temp, last(temp)+1) else temp <- c(temp, last(temp))
  }
  new_dens[["group"]] <- temp

  new_pmf <- sapply(new_dens$group %>% unique, function(which_group) {
    temp <- which(new_dens$group == which_group)
    a1 <- first(temp)
    a2 <- last(temp)
    c(new_dens$outcome[a1],
      # (new_dens$outcome[a2] - new_dens$outcome[a1]) *
        new_dens$cont_dens[a1]
    )
  }) %>% t
  # temp_loc <- last(which(new_pmf[, 2] !=0))
  # new_pmf[temp_loc, 2] <- 1-sum(new_pmf[-temp_loc, 2])
  new_pmf[, 2] %>% sum
  # psi <- sum(new_pmf[, 2]^2)
  # new_pmf %>% points()
  new_pmf %>% lines(type = "s", col = "red")
}


###
# 1-TMLE
###

updater <- tmle3_Update$new(convergence_type = "scaled_var",
                            constrain_step = T,
                            optim_delta_epsilon = F,
                            one_dimensional = F,
                            delta_epsilon = function(x) {
                              ifelse(abs(mean(x %>% as.vector)) < sqrt(var(x %>% as.vector)/tmle_task$nrow)/log(tmle_task$nrow),
                                     0,
                                     ifelse(mean(x %>% as.vector) > 0, 0.01, -0.01)
                              )
                            },
                            maxit=1000
                            ,
                            cvtmle=F)
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood,
                                               submodel_type_by_node = "EIC" ,
                                               updater = updater)
tmle_params <- list(
  Param_ave_dens$new(observed_likelihood = targeted_likelihood)
)
updater$tmle_params <- tmle_params

tmle_params[[1]]$estimates()$IC %>% mean
# abs(mean(tmle_params[[1]]$estimates()$IC %>% as.vector))
sqrt(var(tmle_params[[1]]$estimates()$IC %>% as.vector)/tmle_task$nrow)/log(tmle_task$nrow)
targeted_likelihood$get_likelihood(tmle_task, "outcome", "full") %>% log %>% sum

suppressMessages(
  fit_onestep_100 <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
)
fit_onestep_100

fold_number <- "full"
{
  self <- tmle_params[[1]]
  support_points <- seq(-20, 20, by = 0.01)
  new_data <- data.table(one = 0,
                         outcome = support_points)
  new_task <- tmle3_Task$new(new_data, npsem)
  new_dens <- data.frame(outcome = new_data$outcome,
                         cont_dens = self$observed_likelihood$get_likelihood(new_task, "outcome", fold_number))

  temp <- c(1)
  for (i in 2:nrow(new_data)) {
    if (new_dens$cont_dens[i] != new_dens$cont_dens[i-1]) temp <- c(temp, last(temp)+1) else temp <- c(temp, last(temp))
  }
  new_dens[["group"]] <- temp

  new_pmf <- sapply(new_dens$group %>% unique, function(which_group) {
    temp <- which(new_dens$group == which_group)
    a1 <- first(temp)
    a2 <- last(temp)
    c(new_dens$outcome[a1],
      # (new_dens$outcome[a2] - new_dens$outcome[a1]) *
      new_dens$cont_dens[a1]
    )
  }) %>% t
  # temp_loc <- last(which(new_pmf[, 2] !=0))
  # new_pmf[temp_loc, 2] <- 1-sum(new_pmf[-temp_loc, 2])
  new_pmf[, 2] %>% sum
  # psi <- sum(new_pmf[, 2]^2)
  # new_pmf %>% points()
  new_pmf %>% lines(type = "s", col = "blue")
}

updater$epsilons %>% unlist %>% sum

tmle_params[[1]]$estimates()$IC %>% mean
# abs(mean(tmle_params[[1]]$estimates()$IC %>% as.vector))
sqrt(var(tmle_params[[1]]$estimates()$IC %>% as.vector)/tmle_task$nrow)/log(tmle_task$nrow)
targeted_likelihood$get_likelihood(tmle_task, "outcome", "full") %>% log %>% sum




updater_2 <- tmle3_Update$new(convergence_type = "scaled_var",
                            constrain_step = T,
                            optim_delta_epsilon = F,
                            one_dimensional = F,
                            delta_epsilon = function(x) {
                              ifelse(abs(mean(x %>% as.vector)) < sqrt(var(x %>% as.vector)/tmle_task$nrow)/log(tmle_task$nrow),
                                     0,
                                     ifelse(mean(x %>% as.vector) > 0, 0.01, -0.01)
                              )
                            },
                            maxit=1000
                            ,
                            cvtmle=F)
targeted_likelihood_2 <- Targeted_Likelihood$new(initial_likelihood,
                                               submodel_type_by_node = "EIC" ,
                                               updater = updater_2)
tmle_param_2 <- Param_ave_dens_2$new(observed_likelihood = targeted_likelihood_2, targeted_likelihood_1 = targeted_likelihood, updater_1 = updater, tilde_likelihood = targeted_likelihood)
updater_2$tmle_params <- list(tmle_param_2)

suppressMessages(
  fit_2 <- fit_tmle3(tmle_task, targeted_likelihood_2, updater_2$tmle_params, updater_2)
)

# tmle_param_2$estimates()$IC %>% mean
# # abs(mean(tmle_params[[1]]$estimates()$IC %>% as.vector))
# sqrt(var(tmle_param_2$estimates()$IC %>% as.vector)/tmle_task$nrow)/log(tmle_task$nrow)
# targeted_likelihood_2$get_likelihood(tmle_task, "outcome", "full") %>% log %>% sum

# updater_2$epsilons

fit_2







updater_final <- tmle3_Update$new(convergence_type = "scaled_var",
                            constrain_step = T,
                            optim_delta_epsilon = F,
                            one_dimensional = F,
                            delta_epsilon = function(x) {
                              ifelse(abs(mean(x %>% as.vector)) < sqrt(var(x %>% as.vector)/tmle_task$nrow)/log(tmle_task$nrow),
                                     0,
                                     ifelse(mean(x %>% as.vector) > 0, 0.01, -0.01)
                              )
                            },
                            maxit=1000
                            ,
                            cvtmle=F)
targeted_likelihood_final <- Targeted_Likelihood$new(targeted_likelihood_2,
                                               submodel_type_by_node = "EIC" ,
                                               updater = updater_final)
tmle_params <- list(
  Param_ave_dens$new(observed_likelihood = targeted_likelihood_final)
)
updater$tmle_params <- tmle_params

tmle_params[[1]]$estimates()$IC %>% mean
# abs(mean(tmle_params[[1]]$estimates()$IC %>% as.vector))
sqrt(var(tmle_params[[1]]$estimates()$IC %>% as.vector)/tmle_task$nrow)/log(tmle_task$nrow)
targeted_likelihood$get_likelihood(tmle_task, "outcome", "full") %>% log %>% sum

suppressMessages(
  fit_final <- fit_tmle3(tmle_task, targeted_likelihood_final, tmle_params, updater_final)
)
fit_final


nontargeting$psi
fit_onestep_100
fit_final
{
  self <- tmle_params[[1]]
  support_points <- seq(-20, 20, by = 0.1)
  new_data <- data.table(one = 0,
                         outcome = support_points)
  new_task <- tmle3_Task$new(new_data, npsem)
  new_dens <- data.frame(outcome = new_data$outcome,
                         cont_dens = self$observed_likelihood$get_likelihood(new_task, "outcome", fold_number))

  temp <- c(1)
  for (i in 2:nrow(new_data)) {
    if (new_dens$cont_dens[i] != new_dens$cont_dens[i-1]) temp <- c(temp, last(temp)+1) else temp <- c(temp, last(temp))
  }
  new_dens[["group"]] <- temp

  new_pmf <- sapply(new_dens$group %>% unique, function(which_group) {
    temp <- which(new_dens$group == which_group)
    a1 <- first(temp)
    a2 <- last(temp)
    c(new_dens$outcome[a1],
      # (new_dens$outcome[a2] - new_dens$outcome[a1]) *
      new_dens$cont_dens[a1]
    )
  }) %>% t
  # temp_loc <- last(which(new_pmf[, 2] !=0))
  # new_pmf[temp_loc, 2] <- 1-sum(new_pmf[-temp_loc, 2])
  new_pmf[, 2] %>% sum
  # psi <- sum(new_pmf[, 2]^2)
  # new_pmf %>% points()
  new_pmf %>% lines(type = "s", col = "green")
}
