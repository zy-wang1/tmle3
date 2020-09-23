
#' @export
ipw_middle <- function(task, lik, ipw_args, fold_number){

  cf_likelihood_control = ipw_args$cf_likelihood_control
  cf_likelihood_treatment = ipw_args$cf_likelihood_treatment
  intervention_list_treatment <- ipw_args$intervention_list_treatment
  intervention_list_control <- ipw_args$intervention_list_control
  # todo: extend for stochastic
  cf_task_treatment <- cf_likelihood_treatment$enumerate_cf_tasks(task)[[1]]
  cf_task_control <- cf_likelihood_control$enumerate_cf_tasks(task)[[1]]

  intervention_nodes <- union(names(intervention_list_treatment), names(intervention_list_control))

  temp_node_names <- names(task$npsem)
  loc_A <- grep("A", temp_node_names)
  loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
  if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)

  # get list of all possible predicted lkds
  obs_data <- task$data
  obs_variable_names <- colnames(obs_data)
  # ZW todo: to handle long format and wide format
  temp_to_drop <- which(obs_variable_names %in% c("id", "t"))
  obs_variable_names <- obs_variable_names[-temp_to_drop]
  # ZW todo: see if observed_likelihood needs to change to targeted likelihood

  intervention_variables <- map_chr(task$npsem[intervention_nodes], ~.x$variables)
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
  intervention_levels_treat <- map_dbl(intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
  intervention_levels_control <- map_dbl(intervention_list_control, ~.x$value %>% as.character %>% as.numeric)

  list_H <- get_obs_H(task, obs_data, current_likelihood = lik,
                      cf_task_treatment, cf_task_control,
                      intervention_variables, intervention_levels_treat, intervention_levels_control)
  list_newH <- list()
  for (ind_var in 1:length(list_H)) {
    if(!is.null(list_H[[ind_var]])) {
      current_ind <- (obs_data[[task$npsem[[ind_var]]$variables]] == 1)*1
      if (ind_var %in% loc_Z) temp_p <- lik$get_likelihoods(cf_task_control, temp_node_names[ind_var]) else
        temp_p <- lik$get_likelihoods(cf_task_treatment, temp_node_names[ind_var])
      temp_p <- ifelse(current_ind == 1, temp_p, 1 - temp_p)
      list_newH[[ind_var]] <- ( list_H[[ind_var]] * (current_ind - temp_p)  ) %>% as.matrix
    }
  }
  names(list_newH) <- temp_node_names

  return(list_newH)

  # ZW todo: for categorical variables
}


#' @export
gradient_generator_middle <- function(tmle_task, lik,  node, include_outcome = T, ipw_args = NULL, fold_number){

  task <- tmle_task$get_regression_task(node)
  IC <- ipw_middle(tmle_task, lik,  ipw_args, fold_number)[[node]] %>% as.vector
  cols <- task$add_columns(data.table(IC = IC))
  task <- task$clone()
  nodes <- task$nodes
  nodes$outcome <- "IC"
  nodes$covariates <- c(nodes$covariates, tmle_task$npsem[[node]]$variables)

  task$initialize(
    task$internal_data,
    nodes = nodes,
    folds = task$folds,
    column_names = cols,
    row_index = task$row_index,
    outcome_type = "continuous"
  )
  return(task)
}

