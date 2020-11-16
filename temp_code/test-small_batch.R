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

library(parallel)

n_sim <- 8 * 2
nCores <- 8

timepoint <- 1
if_misspec <- T

data_truth <- generate_Zheng_data(B = 100000, tau = timepoint, seed = 202008, setAM = c(1, 0), if_LY_misspec = if_misspec)
truth <- data_truth[[timepoint + 1]]$Y %>% mean
truth

results <- mclapply(X = 1:n_sim, mc.cores = nCores, FUN = function(i) {

  # set.seed(1234)

  data_sim <- generate_Zheng_data(B = 1000, tau = timepoint, if_LY_misspec = if_misspec)
  data_wide <- data.frame(data_sim)

  node_list <- list(L_0 = c("L1_0", "L2_0"),
                    A_1 = "A_1",
                    R_1 = "R_1",
                    Z_1 = "Z_1",
                    L_1 = "L1_1",
                    Y_1 = "Y_1"
                    # ,
                    # A_2 = "A_2",
                    # R_2 = "R_2",
                    # Z_2 = "Z_2",
                    # L_2 = "L1_2",
                    # Y_2 = "Y_2"
                    # ,
                    # A_3 = "A_3",
                    # R_3 = "R_3",
                    # Z_3 = "Z_3",
                    # L_3 = "L1_3",
                    # Y_3 = "Y_3"
  )
  node_list$L_1 <- "L_1"
  names(data_sim[[2]])[grep("L1_1", names(data_sim[[2]]))] <- "L_1"
  names(data_wide)[grep("L1_1", names(data_wide))] <- "L_1"
  # node_list$L_2 <- "L_2"
  # names(data_sim[[3]])[grep("L1_2", names(data_sim[[3]]))] <- "L_2"
  # names(data_wide)[grep("L1_2", names(data_wide))] <- "L_2"

  med_spec <- tmle_med(
    treatment_level = 1,
    control_level = 0
  )

  tmle_task1 <- med_spec$make_tmle_task(data_wide, node_list)

  # choose base learners
  lrnr_glm_fast <- Lrnr_glm_fast$new(outcome_type = "binomial")
  learner_list <- lapply(1:length(tmle_task1$npsem), function(s) Lrnr_sl$new(
    learners = list(
      lrnr_glm_fast
    )
  ))
  names(learner_list) <- names(tmle_task1$npsem)  # the first will be ignored; empirical dist. will be used for covariates

  initial_likelihood <- med_spec$make_initial_likelihood(
    tmle_task1,
    learner_list
  )

  # # analytic EIC
  # updater <- tmle3_Update$new(convergence_type = "scaled_var",
  #                             constrain_step = T,
  #                             optim_delta_epsilon = T,
  #                             one_dimensional=F,
  #                             delta_epsilon=function(x) {
  #                               ifelse(abs(mean(x %>% as.vector)) < sqrt(var(x %>% as.vector)/tmle_task$nrow)/log(tmle_task$nrow),
  #                                      0.00000001,
  #                                      ifelse(mean(x %>% as.vector) > 0, 0.01, -0.01)
  #                               )
  #                             },
  #                             maxit=20
  #                             ,
  #                             cvtmle=F)
  # tlik <- Targeted_Likelihood$new(initial_likelihood,
  #                                 submodel_type_by_node = "EIC" ,
  #                                 updater = updater)

  # # logistic
  # updater <- tmle3_Update$new(convergence_type = "scaled_var",
  #                             # constrain_step = T,
  #                             # optim_delta_epsilon = T,
  #                             # one_dimensional=F,
  #                             # delta_epsilon=function(x) {
  #                             #   ifelse(abs(mean(x %>% as.vector)) < sqrt(var(x %>% as.vector)/tmle_task$nrow)/log(tmle_task$nrow),
  #                             #          0.00000001,
  #                             #          ifelse(mean(x %>% as.vector) > 0, 1, -1)
  #                             #   )
  #                             # },
  #                             maxit=4
  #                             ,
  #                             cvtmle=F)
  # tlik <- Targeted_Likelihood$new(initial_likelihood,
  #                                 submodel_type_by_node = "logistic" ,
  #                                 updater = updater)

  # projected-EIC, resampling
  updater <- tmle3_Update$new(convergence_type = "scaled_var",
                              constrain_step = T,
                              optim_delta_epsilon = T,
                              one_dimensional=F,
                              delta_epsilon=function(x) {
                                ifelse(abs(mean(x %>% as.vector)) < sqrt(var(x %>% as.vector)/tmle_task1$nrow)/log(tmle_task1$nrow),
                                       0.00000001,
                                       ifelse(mean(x %>% as.vector) > 0, 0.01, -0.01)
                                )
                              },
                              maxit=100
                              ,
                              cvtmle=F)
  tlik <- Targeted_Likelihood$new(initial_likelihood,
                                  submodel_type_by_node = "EIC" ,
                                  updater = updater)


  tmle_params_no <- med_spec$make_params(tmle_task1, tlik, options = list("tc"), static_likelihood = initial_likelihood,
                                         if_projection = T
                                         # , n_resampling = 0
  )

  # projection param
  tmle_params <- med_spec$make_params(tmle_task1, tlik, options = list("tc"), static_likelihood = initial_likelihood,
                                      if_projection = T
                                      , n_resampling = 50000
                                      )  # new est

  # tmle_params_raw <- med_spec$make_params(tmle_task1, initial_likelihood, options = list("tc")
  #                                     # ,
  #                                     # if_projection = T,
  #                                     # if_resampling = T
  #                                     )

  # tmle_params <- med_spec$make_params(tmle_task, tlik, options = list("tc")
  #                                     # ,
  #                                     # if_projection = T,
  #                                     # if_resampling = T
  #                                     )

  updater$tmle_params <- tmle_params

  suppressMessages(
    nontargeting <- tmle_params_no[[1]]$estimates(tmle_task1)
  )
  temp_IC <- nontargeting$IC
  var_D <- var(temp_IC)
  n <- length(temp_IC)
  se <- sqrt(var_D / n)
  CI2 <- nontargeting$psi + 1.96 * se
  CI1 <- nontargeting$psi - 1.96 * se

  suppressMessages(
    new <- tmle_params[[1]]$estimates(tmle_task1)
  )
  temp_IC <- new$IC
  var_D <- var(temp_IC)
  n <- length(temp_IC)
  se <- sqrt(var_D / n)
  CI2_new <- new$psi + 1.96 * se
  CI1_new <- new$psi - 1.96 * se
  new_psi <- new$psi

  # new_fit <- fit_tmle3(tmle_task, tlik, tmle_params, updater)
  # new_fit
  # tlik$cache$tasks %>% lapply(function(x) x$data %>% dim)
  # new_psi <- new_fit$summary$tmle_est
  # CI1_new <- new_fit$summary$lower
  # CI2_new <- new_fit$summary$upper
  # updater$step_number

  # updater$update_step(likelihood = tlik, tmle_task, fold_number = "full")
  # tmle_params[[1]]$clever_covariates()$IC %>% rowSums %>% mean
  # sd(tmle_params[[1]]$clever_covariates()$IC %>% rowSums) / sqrt(tmle_task$nrow) / log(tmle_task$nrow)
  # var(tmle_params[[1]]$clever_covariates()$IC %>% rowSums)
  # var_D
  # tmle_params[[1]]$estimates()$psi
  # nontargeting$psi
  # tmle_params[[1]]$estimates()$IC %>% hist
  # tmle_params_no[[1]]$estimates()$IC %>% hist


  # capture.output(
  #   updater$update(tlik, tmle_task)
  # )
  # new_est <- updater$tmle_params[[1]]$estimates()
  # new_psi <- new_est$psi
  # nontargeting$psi
  # new_psi
  # temp_IC <- new_est$IC
  # var_D <- var(temp_IC)
  # n <- length(temp_IC)
  # se <- sqrt(var_D / n)
  # CI2_new <- new_est$psi + 1.96 * se
  # CI1_new <- new_est$psi - 1.96 * se

  return(
    list(non = nontargeting$psi,
         CI1 = CI1,
         CI2 = CI2,
         new = new_psi,
         CI1_new = CI1_new,
         CI2_new = CI2_new,
         step = updater$step_number
         # ,
         # ed = test$ED,
         # threshold = threshold1
         # ,
         # one_cvtmle = test2$estimates[[1]]$psi,
         # step_cvtmle = test2$steps,
         # ed_cvtmle = test2$ED,
         # threshold_cvtmle = sd(tmle_params[[1]]$estimates()$IC) / sqrt(1000) / log(1000)
         # ,
         # last_dir = updater$record_direction %>% last %>% unlist,
         # dir10 = ifelse_vec(length(updater$record_direction) >= 10, updater$record_direction[[10]] %>% unlist, 0)
    )
  )
})


(do.call(rbind, results) %>% as.matrix)
report <- sapply(1:2, function(i) {
  c(mean(( unlist((do.call(rbind, results) %>% as.matrix)[, i * 3 - 2]) - truth)^2),
  mean(( unlist((do.call(rbind, results) %>% as.matrix)[, i * 3 - 2]) - truth)),
  sd( unlist((do.call(rbind, results) %>% as.matrix)[, i * 3 - 2])),
  mean(( unlist((do.call(rbind, results) %>% as.matrix)[, i * 3 - 1]) <= truth & unlist((do.call(rbind, results) %>% as.matrix)[, i * 3]) >= truth)),
  mean(( unlist((do.call(rbind, results) %>% as.matrix)[, i * 3]) - unlist((do.call(rbind, results) %>% as.matrix)[, i * 3 - 1]) ))
  )
})
rownames(report) <- c("MSE", "bias", "sd", "coverage", "width")
report
(do.call(rbind, results) %>% as.matrix)[, 7] %>% unlist

