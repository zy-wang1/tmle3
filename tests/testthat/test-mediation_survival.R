context("Mediation target parameters, v2; survival")

library(R6)  # R6class
library(data.table)  # setDT
library(sl3)
library(digest)
library(uuid)  # UUIDgenerate
library(delayed)  # bundle_delayed
library(assertthat)  # assert_that
# library(speedglm)  # avoid speedglm; sometimes mess up column ordering
# library(methods)  # is
library(dplyr)  # dplyr::select, to mask other packages
library(purrr)  # map_ functions

code_list <- list.files("./R", full.names = T)
for (code in code_list) source(code)

timepoint <- 2
if_LY_misspec <- F
if_A_misspec <- F

B <- 1E6
data_truth <- generate_Zheng_data_survival(sample_size = B, tau = timepoint, seed = 202008, setAM = c(1, 0),
                                                if_LY_misspec = if_LY_misspec, event_label = 0)
truth <- ((data_truth[[timepoint + 1]]$Y) %>% sum(na.rm = T)) / B  # prob of survival
1-truth  # we use decreasing event process; this is the endpoint death probability


sample_size <- 400
set.seed(123)
data_sim <- generate_Zheng_data_survival(sample_size = sample_size, tau = timepoint, if_LY_misspec = if_LY_misspec, if_A_misspec = if_A_misspec)
data_wide <- data.frame(data_sim)

# match the name for projection
node_list <- lapply(1:timepoint, function(t) paste0(c("A_C", "A_E", "R", "Z", "L", "Y"), "_", t)) %>% unlist %>% as.list
# node_list <- lapply(1:timepoint, function(t) paste0(c("A_C", "A_E", "R", "Z", "L1", "Y"), "_", t)) %>% unlist %>% as.list  # raw names is find for non-projection

names(node_list) <- lapply(1:timepoint, function(t) paste0(c("A_C", "A_E", "R", "Z", "L", "Y"), "_", t)) %>% unlist %>% as.list  # node names
node_list <- c(list(L_0 = c("L1_0", "L2_0")), node_list)  # baseline node

# match the name for projection
names(data_wide) <- as.vector(unlist(node_list))

mediation_spec <- tmle_mediation(
  treatment_level = 1,
  control_level = 0
)
tmle_task <- mediation_spec$make_tmle_task(data_wide, node_list, if_drop_censored = T)

# choose base learners
lrnr_mean <- make_learner(Lrnr_mean)
lrnr_glm <- make_learner(Lrnr_glm, outcome_type = "binomial")
learner_list <- lapply(1:length(tmle_task$npsem), function(s) {
  if (length(grep("Y", names(tmle_task$npsem)[s])) > 0) {
    # learners = list(lrnr_mean)
    learners = list(lrnr_glm)
  } else {
    learners = list(lrnr_glm)
  }
  Lrnr_sl$new(learners = learners)
})
names(learner_list) <- names(tmle_task$npsem)  # the first will be ignored; empirical dist. will be used for covariates

# define obs data initial likelihood
initial_likelihood <- mediation_spec$make_initial_likelihood(
  tmle_task,
  learner_list
)


# no update, analytic EIC
{
  tlik <- initial_likelihood
  tmle_params <- mediation_spec$make_params_survival(tmle_task, tlik, options = list("tc"))
  suppressMessages(
    nontargeting_analytic <- tmle_params[[1]]$estimates(tmle_task)
  )
  temp_IC <- nontargeting_analytic$IC
  var_D <- var(temp_IC)
  n <- length(temp_IC)
  se <- sqrt(var_D / n)
  CI2_no_an <- nontargeting_analytic$psi + 1.96 * se
  CI1_no_an <- nontargeting_analytic$psi - 1.96 * se
}

# onestep update, analytic EIC
{
  # test update
  n_subject <- nrow(tmle_task$data)
  tlik <- Targeted_Likelihood$new(initial_likelihood,
                                  submodel_type_by_node = "EIC" ,
                                  updater = list(convergence_type = "scaled_var",
                                                 constrain_step = T,
                                                 optim_delta_epsilon = T,  # fixed small step_size
                                                 one_dimensional=F,
                                                 delta_epsilon=function(x) {
                                                   ifelse(abs(mean(x %>% as.vector)) < sqrt(var(x %>% as.vector)/n_subject)/log(n_subject),
                                                          0,
                                                          ifelse(mean(x %>% as.vector) > 0, 0.01, -0.01)
                                                   )
                                                 },
                                                 maxit=100
                                                 ,
                                                 cvtmle=F
                                  ))
  tmle_params <- mediation_spec$make_params_survival(tmle_task, tlik)
  ic_0 <- tmle_params[[1]]$clever_covariates()$"IC" %>% colMeans()
  ic_0
  suppressMessages(
    tlik$updater$update(tlik, tmle_task)
  )
  ic_1 <- tmle_params[[1]]$clever_covariates()$"IC" %>% colMeans()
  ic_1
  suppressWarnings(suppressMessages(
    new_est <- tmle_params[[1]]$estimates()$psi
  ))
  new_est

  onestep_an <- tmle_params[[1]]$estimates()
  onestep_an_est <- onestep_an$psi
  temp_IC <- onestep_an$IC
  var_D <- var(temp_IC)
  n <- length(temp_IC)
  se <- sqrt(var_D / n)
  CI2_onestep_an <- onestep_an_est + 1.96 * se
  CI1_onestep_an <- onestep_an_est - 1.96 * se

}

# no update, projected EIC
{
  tlik <- initial_likelihood
  tmle_params_no_re <- mediation_spec$make_params_survival(tmle_task, tlik, options = list("tc")
                                                     , static_likelihood = initial_likelihood,
                                                     if_projection = T
                                                     , n_resampling = 50000
  )  # new est
  suppressMessages(
    nontargeting_re <- tmle_params_no_re[[1]]$estimates(tmle_task)
  )
  temp_IC <- nontargeting_re$IC
  var_D <- var(temp_IC)
  n <- length(temp_IC)
  se <- sqrt(var_D / n)
  CI2_re <- nontargeting_re$psi + 1.96 * se
  CI1_re <- nontargeting_re$psi - 1.96 * se
}

# onestep update, projected EIC
{
  updater <- tmle3_Update$new(convergence_type = "scaled_var",
                              constrain_step = T,
                              optim_delta_epsilon = T,
                              one_dimensional=F,
                              delta_epsilon=function(x) {
                                ifelse(abs(sum(x %>% as.vector)/tmle_task$nrow) <
                                         (sqrt(
                                           (sum(x^2 %>% as.vector)/tmle_task$nrow - (sum(x %>% as.vector)/tmle_task$nrow)^2) * tmle_task$nrow / (tmle_task$nrow - 1)
                                         )/tmle_task$nrow)/log(tmle_task$nrow),
                                       0,
                                       ifelse(mean(x %>% as.vector) > 0, 0.01, -0.01)
                                )
                              },
                              maxit=100
                              ,
                              cvtmle=F)
  tlik <- Targeted_Likelihood$new(initial_likelihood,
                                  submodel_type_by_node = "EIC" ,
                                  updater = updater)
  tmle_params <- mediation_spec$make_params_survival(tmle_task, tlik, options = list("tc")
                                               , static_likelihood = initial_likelihood,
                                               if_projection = T
                                               , n_resampling = 50000
  )

  # projected-EIC, resampling
  new_fit <- fit_tmle3(tmle_task, tlik, tmle_params, updater)
  new_fit
  # tlik$cache$tasks %>% lapply(function(x) x$data %>% dim)
  new_psi <- new_fit$summary$tmle_est
  CI1_new <- new_fit$summary$lower
  CI2_new <- new_fit$summary$upper
  updater$step_number
}








list(c(nontargeting_analytic$psi, CI1_no_an, CI2_no_an),
     c(onestep_an_est, CI1_onestep_an, CI2_onestep_an),
     c(nontargeting_re$psi, CI1_re, CI2_re),
     c(new_psi, CI1_new, CI2_new)
)
