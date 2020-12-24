# check part_A

fit_logistic <- glm(data_wide$A_1 ~ data_wide$L1_0 + data_wide$L2_0, family = "binomial")
fit_logistic
predict(fit_logistic, type = "response")

temp_node_names <- tmle_task$npsem %>% names()
cf_task_treatment <- tmle_params[[1]]$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
cf_task_control <- tmle_params[[1]]$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]
part_A <- lapply(2, function(k) tlik$get_likelihoods(cf_task_treatment, temp_node_names[2], fold_number = "full")) %>% pmap_dbl(prod)  # this is the likelihood of being 1

data.frame(part_A, predict(fit_logistic, type = "response"))[1:20, ]

# check delta_Q

fit_logistic <- glm(data_wide$Y_1 ~ ., data = data_wide, family = "binomial")
fit_logistic
Q_L1_A <- predict(fit_logistic, type = "response")
data_wide$Y_1 - Q_L1_A

temp_node_names <- tmle_task$npsem %>% names()
cf_task_treatment <- tmle_params[[1]]$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
cf_task_control <- tmle_params[[1]]$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]
fold_number = "full"
intervention_nodes <- union(names(tmle_params[[1]]$intervention_list_treatment), names(tmle_params[[1]]$intervention_list_control))
obs_data <- tmle_task$data %>% select(-c(id, t))
obs_variable_names <- colnames(obs_data)
intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
intervention_levels_treat <- map_dbl(tmle_params[[1]]$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
intervention_levels_control <- map_dbl(tmle_params[[1]]$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
names(intervention_levels_treat) <- names(tmle_params[[1]]$intervention_list_treatment)
names(intervention_levels_control) <- names(tmle_params[[1]]$intervention_list_control)

# load full_p list first
full_task <- tmle_task
full_node_names <- names(full_task$npsem)
full_data <- full_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))  # exactly the obs data
full_variable_names <- colnames(full_data)
list_all_predicted_lkd <- lapply(1:length(full_node_names), function(loc_node) {
  if (loc_node > 1) {
    # currently only support univariate node for t>0
    current_variable <- full_task$npsem[[loc_node]]$variables
    temp_input <- expand_values(variables = full_variable_names[1:which(full_variable_names == current_variable)])  # all possible inputs
    temp_task <- tmle3_Task$new(temp_input, full_task$npsem[1:loc_node])
    temp_target_node <- intersect(tmle_params[[1]]$update_nodes, full_node_names[loc_node])
    if (length(temp_target_node) == 1) {
      setattr(temp_task, "target_nodes", full_node_names[loc_node])
      temp_output <- tlik$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
    } else {
      setattr(temp_task, "target_nodes", "no_update")
      temp_output <- tlik$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
    }
    data.frame(temp_input, output = temp_output) %>% return
  }
})
names(list_all_predicted_lkd) <- full_node_names

list_H <- get_obs_H(tmle_task, obs_data, current_likelihood = initial_likelihood,
                    cf_task_treatment, cf_task_control,
                    intervention_variables, intervention_levels_treat, intervention_levels_control,
                    fold_number = "full")
list_Q_1 <- get_obs_Q(tmle_task, obs_data, list_H,
                      intervention_variables, intervention_levels_treat, intervention_levels_control,
                      list_all_predicted_lkd,  # val version decided above for fold_number == "validation"
                      lt = 1)
list_Q_0 <- get_obs_Q(tmle_task, obs_data, list_H,
                      intervention_variables, intervention_levels_treat, intervention_levels_control,
                      list_all_predicted_lkd,
                      lt = 0)
list_Q_1 %>% last
list_Q_0 %>% last

list_newH <- list()
for (loc_node in 1:length(list_H)) {
  if(!is.null(list_H[[loc_node]])) {  # for density update, change the sign of some clever covariates
    temp_vec <- list_H[[loc_node]] * (list_Q_1[[loc_node]] - list_Q_0[[loc_node]])
    list_newH[[loc_node]] <- temp_vec
  }
}
names(list_newH) <- temp_node_names


# check epsilon here
H <- list_newH$Y_1
loc_node <- length(temp_node_names)
current_ind <- (obs_data[[tmle_task$npsem[[loc_node]]$variables]] == 1)*1
temp_vec <- list_newH[[loc_node]] %>% as.vector
temp_p <- tlik$get_likelihoods(tmle_task, temp_node_names[loc_node], fold_number)
# if (loc_node %in% loc_Z) temp_p <- self$observed_likelihood$get_likelihoods(cf_task_control, temp_node_names[loc_node], fold_number) else
#   temp_p <- self$observed_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[loc_node], fold_number)
temp_p <- ifelse(current_ind == 1, temp_p, 1 - temp_p)  # transform density to conditional mean
temp_p

# glm(current_ind ~ H - 1 + offset(expit(temp_p)), family = "binomial")
glm(current_ind ~ H - 1 + offset(logit(temp_p)), family = "binomial")
glm(current_ind ~ H - 1, offset = logit(temp_p), family = "binomial")

submodel_data <- updater$generate_submodel_data(
  tlik, tmle_task,
  fold_number, "Y_1",
  drop_censored = TRUE,
  for_fitting = T
)
new_epsilon <- updater$fit_submodel(submodel_data)
new_epsilon

tmle_task$data$A_1 %>% head
tmle_params[[1]]$clever_covariates()[["Y_1"]] %>% head
submodel_data$H %>% head
H %>% head
data.frame(submodel_data$initial,
           ifelse(current_ind == 1, temp_p, 1-temp_p),
           current_ind)

submodel_info <- submodel_spec("logistic")
submodel_info$offset_tranform(temp_p) %>% head
logit(temp_p) %>% head
expit(temp_p) %>% head


new_fit <- fit_tmle3(tmle_task, tlik, tmle_params, updater)
new_fit
updater$epsilons

# for (loc_node in 1:length(list_H)) {
#   if(!is.null(list_H[[loc_node]])) {  # for density update, change the sign of some clever covariates
#     temp_vec <- list_H[[loc_node]] * (list_Q_1[[loc_node]] - list_Q_0[[loc_node]])
#     list_newH[[loc_node]] <- temp_vec
#   }
# }

loc_node <- length(temp_node_names)
current_ind <- (obs_data[[tmle_task$npsem[[loc_node]]$variables]] == 1)*1
temp_vec <- list_newH[[loc_node]] %>% as.vector
temp_p <- tlik$get_likelihoods(tmle_task, temp_node_names[loc_node], fold_number)
# if (loc_node %in% loc_Z) temp_p <- self$observed_likelihood$get_likelihoods(cf_task_control, temp_node_names[loc_node], fold_number) else
#   temp_p <- self$observed_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[loc_node], fold_number)
temp_p <- ifelse(current_ind == 1, temp_p, 1 - temp_p)  # transform density to conditional mean
# list_D[[loc_node]] <- (current_ind - temp_p) * temp_vec
data.frame(current_ind - temp_p, data_wide$Y_1 - Q_L1_A) %>% head(20)





# check ratio part
fit_logistic <- glm(Z_1 ~ L1_0 + L2_0 + A_1 + R_1, data = data_wide, family = "binomial")
fit_logistic
data_wide$A_1

newdata_t <- cf_task_treatment$data %>% select(-c(id, t)) %>% as.matrix %>% as.data.frame
newdata_c <- cf_task_control$data %>% select(-c(id, t)) %>% as.matrix %>% as.data.frame

pZ_t <- predict.glm(fit_logistic, type = "response", newdata = newdata_t)
pZ_c <- predict.glm(fit_logistic, type = "response", newdata = newdata_c)
pZ_t <- ifelse(tmle_task$get_tmle_node("Z_1") == 1, pZ_t, 1-pZ_t)
pZ_c <- ifelse(tmle_task$get_tmle_node("Z_1") == 1, pZ_c, 1-pZ_c)
pZ_t %>% head(20)
pZ_c %>% head(20)

pZ_c/pZ_t

temp_node_names[4]
part_Z <- lapply(4, function(k) {
  tlik$get_likelihoods(cf_task_control, temp_node_names[k], fold_number) /
    tlik$get_likelihoods(cf_task_treatment, temp_node_names[k], fold_number)
}) %>% pmap_dbl(prod)

data.frame(part_Z, pZ_c/pZ_t) %>% head





# check updated tlik
# load full_p list first
full_task <- tmle_task
full_node_names <- names(full_task$npsem)
full_data <- full_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))  # exactly the obs data
full_variable_names <- colnames(full_data)
list_all_predicted_lkd <- lapply(1:length(full_node_names), function(loc_node) {
  if (loc_node > 1) {
    # currently only support univariate node for t>0
    current_variable <- full_task$npsem[[loc_node]]$variables
    temp_input <- expand_values(variables = full_variable_names[1:which(full_variable_names == current_variable)])  # all possible inputs
    temp_task <- tmle3_Task$new(temp_input, full_task$npsem[1:loc_node])
    temp_target_node <- intersect(tmle_params[[1]]$update_nodes, full_node_names[loc_node])
    if (length(temp_target_node) == 1) {
      setattr(temp_task, "target_nodes", full_node_names[loc_node])
      temp_output <- tlik$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
    } else {
      setattr(temp_task, "target_nodes", "no_update")
      temp_output <- tlik$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
    }
    data.frame(temp_input, output = temp_output) %>% return
  }
})
names(list_all_predicted_lkd) <- full_node_names


data.frame(obs_update = tlik$get_likelihood(tmle_task, node = "Z_1"),
temp_update = left_join(tmle_task$data %>% select(-c(id, t)),
          list_all_predicted_lkd[["Z_1"]]
          )$output) %>% head(30)

(left_join(tmle_task$data %>% select(-c(id, t)),
          list_all_predicted_lkd[["Z_1"]]
)$output - tlik$get_likelihood(tmle_task, node = "Z_1")) %>% range
