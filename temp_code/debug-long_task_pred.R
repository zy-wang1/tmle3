# to recreate the issue; setwd() to the tmle3 folder
{
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

  timepoint <- 1
  if_misspec <- T
  n_subject <- 70

  set.seed(12)
  data_sim <- generate_Zheng_data(B = n_subject, tau = timepoint, if_LY_misspec = if_misspec)
  data_wide <- data.frame(data_sim)
  node_list <- list(L_0 = c("L1_0", "L2_0"),
                    A_1 = "A_1",
                    R_1 = "R_1",
                    Z_1 = "Z_1",
                    L_1 = "L1_1",
                    Y_1 = "Y_1"
  )
  node_list$L_1 <- "L_1"
  names(data_sim[[2]])[grep("L1_1", names(data_sim[[2]]))] <- "L_1"
  names(data_wide)[grep("L1_1", names(data_wide))] <- "L_1"
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
  tlik <- Targeted_Likelihood$new(initial_likelihood,
                                  submodel_type_by_node = "EIC" ,
                                  updater = list(convergence_type = "scaled_var",
                                                 constrain_step = T,
                                                 optim_delta_epsilon = T,
                                                 one_dimensional=T,
                                                 delta_epsilon=function(x) {
                                                   ifelse(abs(mean(x %>% as.vector)) < sqrt(var(x %>% as.vector)/n_subject)/log(n_subject),
                                                          0.00000001,
                                                          ifelse(mean(x %>% as.vector) > 0, 0.01, -0.01)
                                                   )
                                                 },
                                                 maxit=100
                                                 ,
                                                 cvtmle=F
                                  ))
  tmle_params <- middle_spec$make_params(tmle_task, tlik, if_projection = T, initial_likelihood)
}





# create a clean library of inputs and outputs of pY
self <- tmle_params[[1]]
loc_node <- length(temp_node_names)
node <- temp_node_names[loc_node]
current_variable <- tmle_task$npsem[[loc_node]]$variables
temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)])  # all possible inputs
temp_task <- tmle3_Task$new(temp_input, tmle_task$npsem[1:loc_node])
likelihood_factor <- tlik$factor_list[[node]]
clean_library <- data.frame(
  temp_task$data %>% select(-c(id, t)),
  value = likelihood_factor$get_likelihood(temp_task, fold_number, expand = expand, node = node)
)

# create a library with long_task as inputs, and the outputs from get_likelihood
long_task <- gradient$expand_task(tmle_task, node = "Y_1")
to_verify_library <- data.frame(
  long_task$data %>% select(-c(id, t, trueid)),
  value = likelihood_factor$get_likelihood(long_task, fold_number, expand = expand, node = node)
)
# on cluster the two lines below agree; on R3.6.3 and R3.6.1 on my laptop they do not
left_join(long_task$data %>% select(-c(id, t, trueid)) %>% as.data.frame(), clean_library %>% as.data.frame())$value
to_verify_library$value

# there is no issue with non-expended tasks, such as the observed tmle_task
to_verify_library <- data.frame(
  tmle_task$data %>% select(-c(id, t)),
  value = likelihood_factor$get_likelihood(tmle_task, fold_number, expand = expand, node = node)
)
# these two always agree
left_join(tmle_task$data %>% select(-c(id, t)) %>% as.data.frame(), clean_library %>% as.data.frame())$value
to_verify_library$value
