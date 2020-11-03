library(parallel)

n_sim <- 8 *2
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
  )

  middle_spec <- tmle_middle(
    treatment_level = 1,
    control_level = 0
  )

  tmle_task <- middle_spec$make_tmle_task(data_wide, node_list)

  # choose base learners
  lrnr_glm_fast <- Lrnr_glm_fast$new(outcome_type = "binomial")
  learner_list <- lapply(1:length(tmle_task$npsem), function(s) Lrnr_sl$new(
    learners = list(
      lrnr_glm_fast
    )
  ))
  names(learner_list) <- names(tmle_task$npsem)  # the first will be ignored; empirical dist. will be used for covariates

  initial_likelihood <- middle_spec$make_initial_likelihood(
    tmle_task,
    learner_list
  )




  # updater <- tmle3_Update_middle$new(maxit = 1, convergence_type = "sample_size",
  #                                    fluctuation_type = "standard", submodel_type = "logistic"
  #                                    # , cvtmle = T
  # )
  # targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
  # tmle_params <- middle_spec$make_params(tmle_task, targeted_likelihood)
  # updater$tmle_params <- tmle_params
  #
  # tmle_params[[1]]$estimates()
  # tmle_params[[1]]$list_D %>% lapply(mean)
  # suppressMessages(
  #   test <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
  # )
  # test
  # tmle_params[[1]]$list_D %>% lapply(mean)




  updater <- tmle3_Update_middle$new(maxit = 100, convergence_type = "scaled_var",
                                     fluctuation_type = "standard", submodel_type = "onestep", d_epsilon = 0.01,
                                     cvtmle = F
                                     , if_direction = T
  )
  targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
  tmle_params <- middle_spec$make_params(tmle_task, targeted_likelihood)
  updater$tmle_params <- tmle_params

  suppressMessages(
    tmle_params[[1]]$estimates()$psi
  )
  list1 <- tmle_params[[1]]$list_D %>% lapply(mean) %>% compact %>% unlist

  suppressMessages(
    test <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
  )
  test
  test$ED
  threshold1 <- sd(tmle_params[[1]]$estimates()$IC) / sqrt(1000) / log(1000)
  test$steps


  updater <- tmle3_Update_middle$new(maxit = 100, convergence_type = "scaled_var",
                                     fluctuation_type = "standard", submodel_type = "onestep", d_epsilon = 0.01,
                                     cvtmle = T
                                     , if_direction = T
  )
  targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
  tmle_params <- middle_spec$make_params(tmle_task, targeted_likelihood)
  updater$tmle_params <- tmle_params

  suppressMessages(
    tmle_params[[1]]$estimates()$psi
  )
  list1 <- tmle_params[[1]]$list_D %>% lapply(mean) %>% compact %>% unlist

  suppressMessages(
    test2 <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
  )

  return(
    list(non = test$initial_psi,
         one = test$estimates[[1]]$psi,
         step = test$steps,
         ed = test$ED,
         threshold = threshold1
         ,
         one_cvtmle = test2$estimates[[1]]$psi,
         step_cvtmle = test2$steps,
         ed_cvtmle = test2$ED,
         threshold_cvtmle = sd(tmle_params[[1]]$estimates()$IC) / sqrt(1000) / log(1000)
         # ,
         # last_dir = updater$record_direction %>% last %>% unlist,
         # dir10 = ifelse_vec(length(updater$record_direction) >= 10, updater$record_direction[[10]] %>% unlist, 0)
         )
  )
})

data.frame(
  step = do.call(rbind, results)[, 3] %>% unlist,
  ed = do.call(rbind, results)[, 4] %>% unlist,
  threshold = do.call(rbind, results)[, 5] %>% unlist,
  step_cvtmle = do.call(rbind, results)[, 7] %>% unlist,
  ed_cvtmle = do.call(rbind, results)[, 8] %>% unlist,
  threshold_cvtmle = do.call(rbind, results)[, 9] %>% unlist
)

(do.call(rbind, results) %>% as.matrix)[, -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth)^2 %>% mean)
(do.call(rbind, results) %>% as.matrix)[, -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth) %>% sd)
(do.call(rbind, results) %>% as.matrix)[, -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth) %>% abs %>% mean)

(do.call(rbind, results) %>% as.matrix)[, c(1, 2, 6)] %>% apply(2, function(x) (as.numeric(x) - truth)^2 %>% mean)


(do.call(rbind, results) %>% as.matrix)[(do.call(rbind, results)[, 3] %>% unlist) < 150, -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth)^2 %>% mean)
(do.call(rbind, results) %>% as.matrix)[(do.call(rbind, results)[, 3] %>% unlist) < 150, -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth) %>% sd)
(do.call(rbind, results) %>% as.matrix)[(do.call(rbind, results)[, 3] %>% unlist) < 150, -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth) %>% abs %>% mean)


results
do.call(rbind, results)
results[[1]]$last_dir %>% is.null()

if_Z <- do.call(rbind, map(results, ~ifelse_vec(is.null(.x$last_dir), rep(0, 6), .x$last_dir)))[, 4] > 0
if_Z_at_step_10 <- do.call(rbind, map(results, ~ifelse_vec(is.null(.x$dir10), rep(0, 6), .x$dir10)))[, 4] > 0
(do.call(rbind, results) %>% as.matrix)[if_Z, -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth)^2 %>% mean)
(do.call(rbind, results) %>% as.matrix)[if_Z, -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth) %>% sd)
(do.call(rbind, results) %>% as.matrix)[if_Z, -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth) %>% abs %>% mean)
(do.call(rbind, results) %>% as.matrix)[!if_Z, -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth)^2 %>% mean)
(do.call(rbind, results) %>% as.matrix)[!if_Z, -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth) %>% sd)
(do.call(rbind, results) %>% as.matrix)[!if_Z, -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth) %>% abs %>% mean)


if_converge <- (do.call(rbind, results)[, 3] %>% unlist) < 50
(do.call(rbind, results) %>% as.matrix)[if_converge, -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth)^2 %>% mean)
(do.call(rbind, results) %>% as.matrix)[if_converge, -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth) %>% sd)
(do.call(rbind, results) %>% as.matrix)[if_converge, -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth) %>% abs %>% mean)
(do.call(rbind, results) %>% as.matrix)[!if_converge, -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth)^2 %>% mean)
(do.call(rbind, results) %>% as.matrix)[!if_converge, -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth) %>% sd)
(do.call(rbind, results) %>% as.matrix)[!if_converge, -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth) %>% abs %>% mean)
(do.call(rbind, results) %>% as.matrix)[!if_converge, -c(3:length(results))]

sum(!if_converge)
sum((if_Z) & (!if_converge))
sum((!if_Z) & (!if_converge))
(do.call(rbind, results) %>% as.matrix)[(!if_Z) & (!if_converge), -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth)^2 %>% mean)
(do.call(rbind, results) %>% as.matrix)[(!if_Z) & (!if_converge), -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth) %>% sd)
(do.call(rbind, results) %>% as.matrix)[(!if_Z) & (!if_converge), -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth) %>% abs %>% mean)
(do.call(rbind, results) %>% as.matrix)[(if_Z) & (!if_converge), -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth)^2 %>% mean)
(do.call(rbind, results) %>% as.matrix)[(if_Z) & (!if_converge), -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth) %>% sd)
(do.call(rbind, results) %>% as.matrix)[(if_Z) & (!if_converge), -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth) %>% abs %>% mean)


(do.call(rbind, results) %>% as.matrix)[(if_Z) & (if_converge), -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth)^2 %>% mean)
(do.call(rbind, results) %>% as.matrix)[(if_Z) & (if_converge), -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth) %>% sd)
(do.call(rbind, results) %>% as.matrix)[(if_Z) & (if_converge), -c(3:length(results))] %>% apply(2, function(x) (as.numeric(x) - truth) %>% abs %>% mean)
