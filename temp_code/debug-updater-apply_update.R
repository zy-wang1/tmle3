self <- tlik$updater
tasks_at_step %>% length
temp_task <- tasks_at_step[[2]]


# get submodel data for all nodes
submodel_data <- self$generate_submodel_data(
  tlik, temp_task,
  fold_number, update_node,
  drop_censored = FALSE
)

tlik$cache$tasks %>% length
tlik$cache$tasks %>% lapply(function(x) x$data$A_1 %>% table)

updated_likelihood <- self$apply_submodel(submodel_data, new_epsilon)

if (any(!is.finite(updated_likelihood))) {
  stop("Likelihood was updated to contain non-finite values.\n
             This is likely a result of unbounded likelihood factors")
}
# un-scale to handle bounded continuous
updated_likelihood <- temp_task$unscale(
  updated_likelihood,
  update_node
)
