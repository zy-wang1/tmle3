timepoint <- 1
if_misspec <- T
# data_truth <- generate_Zheng_data(B = 100000, tau = timepoint, seed = 202008, setAM = c(1, 0), if_LY_misspec = if_misspec)
# truth <- data_truth[[timepoint + 1]]$Y %>% mean
# truth

# set.seed(1234)

# saveRDS(data_sim, file = "./temp_code/example_data.rds")
data_sim <- readRDS("./temp_code/example_data.rds")


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
                                   fluctuation_type = "standard", submodel_type = "onestep", d_epsilon = 0.01, cvtmle = F
                                   , if_direction = T
                                   )
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
tmle_params <- middle_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params

suppressMessages(
  tmle_params[[1]]$estimates()$psi
)
list1 <- tmle_params[[1]]$list_D %>% lapply(mean) %>% compact %>% unlist

# tmle_params[[1]]$clever_covariates() %>% last %>% head
# suppressMessages(
# tmle_params[[1]]$clever_covariates(update = T) %>% last %>% head(10)
# )
# observed <- tmle_task$get_tmle_node("Z_1")
# values <- targeted_likelihood$get_likelihood(tmle_task = targeted_likelihood$cache$tasks[[5]], fold_number = "full", node = "Z_1")
# values <- ifelse(observed == 1, values, 1 - values)
# values %>% log %>% sum

# observed <- tmle_task$get_tmle_node("Y_1")
# values <- targeted_likelihood$get_likelihood(tmle_task = targeted_likelihood$cache$tasks[[6]], fold_number = "full", node = "Y_1")
# values <- ifelse(observed == 1, values, 1 - values)
# values %>% log %>% sum


suppressMessages(
  test <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
)
# observed <- tmle_task$get_tmle_node("Z_1")
# values <- targeted_likelihood$get_likelihood(tmle_task = targeted_likelihood$cache$tasks[[5]], fold_number = "full", node = "Z_1")
# values <- ifelse(observed == 1, values, 1 - values)
# values %>% log %>% sum

# observed <- tmle_task$get_tmle_node("Y_1")
# values <- targeted_likelihood$get_likelihood(tmle_task = targeted_likelihood$cache$tasks[[6]], fold_number = "full", node = "Y_1")
# values <- ifelse(observed == 1, values, 1 - values)
# values %>% log %>% sum



test
test$ED
sd(tmle_params[[1]]$estimates()$IC) / sqrt(1000) / log(1000)

# tmle_params[[1]]$clever_covariates() %>% last %>% head

do.call(rbind, updater$record_direction) %>% tail

list2 <- tmle_params[[1]]$list_D %>% lapply(mean) %>% compact %>% unlist
cbind(list1[!is.na(list1)], list2[!is.na(list2)], updater$list_directions_ini %>% unlist)

updater$list_directions_ini

#
#
#
#
# # self$updater$update(self$likelihood, self$tmle_task)
# self <- updater
# update_fold <- self$update_fold
# maxit <- 3
# for (steps in seq_len(maxit)) {
#   self$update_step(targeted_likelihood, tmle_task, update_fold)
#   if (self$check_convergence(tmle_task, update_fold)) {
#     break
#   }
# }

update_fold <- "full"

record <- list()
tprob_record <- list()

updater$update_step(targeted_likelihood, tmle_task, update_fold)
record[[updater$step_number + 1]] <- ED_from_estimates(tmle_params %>% lapply(function(x) x$estimates()))

# tprob_record[[updater$step_number + 1]] <- data.frame(list_tail = left_join(tmle_task$data[, 1:6], tmle_params[[1]]$list_all_predicted_lkd$L_1)$output %>% tail %>% as.vector,
#                                                       tlik_tail = targeted_likelihood$get_likelihood(tmle_params[[1]]$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]], "L_1") %>% tail)

record
# tprob_record
sd(tmle_params[[1]]$estimates()$IC) / sqrt(1000) / log(1000)
tmle_params[[1]]$estimates()$psi


tmle_params[[1]]$list_D %>% lapply(mean) %>% compact() %>% unlist %>% sum(na.rm = T)


initial_likelihood$get_likelihood(tmle_task, "R_1") %>% tail




obs_update_nodes <- do.call(cbind, update_nodes %>% lapply(tmle_task$get_tmle_node))

# tmle_params[[1]]$estimates()$psi
list1 <- tmle_params[[1]]$list_D %>% lapply(function(x) mean(log(1 + x * 0.01)))
loglkd1 <- targeted_likelihood$get_likelihoods(tmle_task, node = update_nodes) %>% log %>% colSums(na.rm = T)
  # sapply(1:ncol(obs_update_nodes), function(i) {
  #   (targeted_likelihood$get_likelihoods(tmle_task, node = update_nodes[i]) %>% as.vector)[obs_update_nodes[, i] == 1] %>% log %>% sum
  # })
updater$update_step(targeted_likelihood, tmle_task, update_fold)
loglkd2 <-
#   sapply(1:ncol(obs_update_nodes), function(i) {
#   (targeted_likelihood$get_likelihoods(tmle_task, node = update_nodes[i]) %>% as.vector)[obs_update_nodes[, i] == 1] %>% log %>% sum
# })
  targeted_likelihood$get_likelihoods(tmle_task, node = update_nodes) %>% log %>% colSums(na.rm = T)
list2 <- tmle_params[[1]]$list_D %>% lapply(function(x) mean(log(1 + x * 0.01)))
updater$update_step(targeted_likelihood, tmle_task, update_fold)
loglkd3 <-
#   sapply(1:ncol(obs_update_nodes), function(i) {
#   (targeted_likelihood$get_likelihoods(tmle_task, node = update_nodes[i]) %>% as.vector)[obs_update_nodes[, i] == 1] %>% log %>% sum
# })
  targeted_likelihood$get_likelihoods(tmle_task, node = update_nodes) %>% log %>% colSums(na.rm = T)
list3 <- tmle_params[[1]]$list_D %>% lapply(function(x) mean(log(1 + x * 0.01)))

updater$list_directions_ini %>% unlist
# list1 %>% unlist
# list2 %>% unlist
# list3 %>% unlist
data.frame(list1 %>% unlist, list2 %>% unlist, list3 %>% unlist)
data.frame(loglkd1, loglkd2, loglkd3)

tmle_params[[1]]$estimates()$psi

initial_likelihood$get_likelihood(tmle_task, "R_1") %>% tail
left_join(tmle_task$data[, 1:4], tmle_params[[1]]$list_all_predicted_lkd$R_1)$output %>% tail %>% as.vector
targeted_likelihood$get_likelihood(tmle_task, "R_1") %>% tail


tmle_params[[1]]$list_D %>% lapply(mean) %>% unlist
abs(unlist(list2)) - abs(unlist(list1))


private$.steps <- self$updater$steps


estimates <- lapply(
  self$tmle_params,
  function(tmle_param) {
    tmle_param$estimates(self$tmle_task, self$updater$update_fold)
  }
)

private$.estimates <- estimates
private$.ED <- ED_from_estimates(estimates)
