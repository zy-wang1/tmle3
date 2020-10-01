self <- tlik$updater

update_nodes <- self$update_nodes

current_step <- self$step_number + 1

# initialize epsilons for this step
na_epsilons <- as.list(rep(NA, length(update_nodes)))
names(na_epsilons) <- update_nodes
epsilons <- self$epsilons
epsilons[[current_step]] <- na_epsilons

tlik$cache$tasks %>% length
tlik$cache$tasks %>% lapply(function(x) x$data$A_1 %>% table)

for (update_node in update_nodes) {
  # get new submodel fit
  submodel_data <- self$generate_submodel_data(
    tlik, tmle_task,
    fold_number = "full", update_node,
    drop_censored = TRUE,
    for_fitting = T
  )

  tlik$cache$tasks %>% length
  tlik$cache$tasks %>% lapply(function(x) x$data$A_1 %>% table)

  new_epsilon <- self$fit_submodel(submodel_data)

  tlik$cache$tasks %>% length
  tlik$cache$tasks %>% lapply(function(x) x$data$A_1 %>% table)


  # update likelihoods
  tlik$update(new_epsilon, current_step, fold_number = "full", update_node)
#
#   tlik$cache$tasks %>% length
#   tlik$cache$tasks %>% lapply(function(x) x$data$A_1 %>% table)
#
#
#   if (fold_number != "full") {
#     # update full fit likelihoods if we haven't already
#     likelihood$update(new_epsilon, self$step_number, "full", update_node)
#   }
#
#   private$.epsilons[[current_step]][[update_node]] <- new_epsilon
# }
#
# # update step number
# private$.step_number <- current_step
