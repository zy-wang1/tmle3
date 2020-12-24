results <- mclapply(X = 1:n_sim, mc.cores = nCores, FUN = function(i) {

  set.seed(1234)
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
  node_list$L_1 <- "L_1"
  names(data_sim[[2]])[grep("L1_1", names(data_sim[[2]]))] <- "L_1"
  names(data_wide)[grep("L1_1", names(data_wide))] <- "L_1"
  # node_list$L_2 <- "L_2"
  # names(data_sim[[3]])[grep("L1_2", names(data_sim[[3]]))] <- "L_2"
  # names(data_wide)[grep("L1_2", names(data_wide))] <- "L_2"
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
  tmle_params_ini <- middle_spec$make_params(tmle_task, initial_likelihood, static_likelihood = NULL)
  suppressMessages(
    nontargeting <- tmle_params_ini[[1]]$estimates(tmle_task)
  )
  # tmle_params_ini[[1]]$estimates()$IC %>% mean
  # (tmle_params_ini[[1]]$estimates()$full_IC %>% lapply(mean) %>% compact() %>% unlist)[-c(1, 2)]


  tlik <- Targeted_Likelihood$new(initial_likelihood,
                                  submodel_type_by_node = "EIC" ,
                                  updater = list(convergence_type = "scaled_var",
                                                 constrain_step = T,
                                                 optim_delta_epsilon = T,
                                                 one_dimensional=T,
                                                 delta_epsilon=function(x) {
                                                   ifelse(abs(mean(x %>% as.vector)) < sqrt(var(x %>% as.vector)/1000)/log(1000),
                                                          0.00000001,
                                                          ifelse(mean(x %>% as.vector) > 0, 0.01, -0.01)
                                                   )
                                                   # ifelse(mean(x %>% as.vector) > 0, 0.0001, -0.0001)
                                                 },
                                                 maxit=100
                                                 ,
                                                 cvtmle=F
                                  ))
  tmle_params <- middle_spec$make_params(tmle_task, tlik, if_projection = T, initial_likelihood)
  est_ini <- tmle_params[[1]]$estimates()$psi
  list1 <- tmle_params[[1]]$clever_covariates()$IC %>% colMeans

  capture.output(
    tlik$updater$update(tlik, tmle_task)
  )
  list2 <- tmle_params[[1]]$clever_covariates()$IC %>% colMeans
  tlik$updater$step_number
  onestep_projected <- tmle_params[[1]]$estimates()
  onestep_projected_est <- onestep_projected$psi
  onestep_projected_est
  nontargeting$psi

  temp_IC <- onestep_projected$IC
  var_D <- var(temp_IC)
  n <- length(temp_IC)
  se <- sqrt(var_D / n)
  CI2_onestep_projected <- onestep_projected_est + 1.96 * se
  CI1_onestep_projected <- onestep_projected_est - 1.96 * se


  return(
    list(non = est_ini,
         one = onestep_projected_est,
         step = tlik$updater$step_number
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
