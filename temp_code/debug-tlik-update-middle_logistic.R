self <- targeted_likelihood

tasks_at_step <- self$cache$tasks

task_updates <- lapply(tasks_at_step, self$updater$apply_update, self, fold_number, new_epsilon)  # this returns updated (obs) lkd
task_updates <- lapply(tasks_at_step[5:7], self$updater$apply_update, self, fold_number, new_epsilon)  # this returns updated (obs) lkd

new_epsilon <- new_epsilons  # get it from updater-update_step
updated_values <- self$updater$apply_update(tasks_at_step[[1]], self, fold_number, new_epsilon)
updated_values <- self$updater$apply_update(tasks_at_step[[5]], self, fold_number, new_epsilon)


# then, store all updates
for (task_index in seq_along(tasks_at_step)) {
  task <- tasks_at_step[[task_index]]
  updated_values <- task_updates[[task_index]]
  # updated_values <- updated_values
  for (node in names(updated_values)) {
    likelihood_factor <- self$factor_list[[node]]
    self$cache$set_values(likelihood_factor, task, step_number + 1, fold_number, updated_values[[node]])
  }
}

update_nodes <- self$updater$update_nodes
full_updates <- self$updater$apply_update_full(self$training_task, self, fold_number, new_epsilon)  # this returns updated (full) lkd
for (update_node in update_nodes) private$.list_all_predicted_lkd[[update_node]]$output <- full_updates[[update_node]]
