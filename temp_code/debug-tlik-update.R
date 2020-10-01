self <- tlik

# todo: rethink which tasks need updates here
# tasks_at_step <- self$cache$tasks_at_step(step_number)
tasks_at_step <- self$cache$tasks
# tasks_at_step <- tasks_at_step[map_dbl(tasks_at_step, ~nrow(.x$data)) == nrow(self$training_task$data)]
tasks_at_step %>% length

tlik$cache$tasks %>% length
tlik$cache$tasks %>% lapply(function(x) x$data$A_1 %>% table)

# if (!inherits(self$updater, "tmle3_Update_middle")  ) {
#   # the default version
#   # If task has attr target_nodes then only update these nodes.
  to_update <- sapply(tasks_at_step, function(task) {
    target_nodes <- attr(task, "target_nodes")
    if(is.null(target_nodes)){
      return(T)
    }
    if(update_node %in% c(target_nodes)){
      return(T)
    }
    return(F)
  })

  to_update

  tasks_at_step <- tasks_at_step[to_update]
  tasks_at_step %>% length

  tasks_at_step %>% length
  tasks_at_step %>% lapply(function(x) x$data$A_1 %>% table)

  lapply(tasks_at_step[3:4], self$updater$apply_update, self, fold_number = "full", new_epsilon, update_node)

  # first, calculate all updates
  task_updates <- lapply(tasks_at_step, self$updater$apply_update, self, fold_number = "full", new_epsilon, update_node)
  # tasks_at_step[[6]] %>% self$updater$apply_update(self, fold_number = "full", new_epsilon, update_node)
  tlik$cache$tasks %>% length
  tlik$cache$tasks %>% lapply(function(x) x$data$A_1 %>% table)
#
# fold_number <- "full"
# step_number <- tlik$updater$step_number
#   # then, store all updates
#   for (task_index in seq_along(tasks_at_step)) {
#     task <- tasks_at_step[[task_index]]
#     updated_values <- task_updates[[task_index]]
#
#     likelihood_factor <- self$factor_list[[update_node]]
#     self$cache$set_values(likelihood_factor, task, step_number + 1, fold_number, updated_values, node = update_node)
#   }
#
#   self$cache$tasks %>% length
