self <- updater

update_nodes <- self$update_nodes

# TODO: support not getting observed for case where we're applying
#       updates instead of fitting them
clever_covariates <- lapply(self$tmle_params, function(tmle_param) {
  tmle_param$clever_covariates(task, fold_number
                               # , update = T
  )  # this returns the list of H
})

observed_values <- lapply(update_nodes, task$get_tmle_node)



all_submodels <- lapply(update_nodes, function(update_node) {
  if (update_node %in% names(task$npsem)) {
    node_covariates <- lapply(clever_covariates, `[[`, update_node)
    covariates_dt <- do.call(cbind, node_covariates)  # there might be multiple targets
    # if (self$one_dimensional) {
    #   observed_task <- likelihood$training_task
    #   estimates <- lapply(self$tmle_params, function(tmle_param) {
    #     tmle_param$estimates(observed_task, fold_number)
    #   })
    #   covariates_dt <- self$collapse_covariates(estimates, covariates_dt)
    # }

    observed <- task$get_tmle_node(update_node)  # raw data on this node
    initial <- targeted_likelihood$get_likelihood(
      task, update_node,
      fold_number
    )  # observed (updated or initial) likelihood at this node
    initial <- ifelse(observed == 1, initial, 1 - initial)

    # scale observed and predicted values for bounded continuous
    observed <- task$scale(observed, update_node)
    initial <- task$scale(initial, update_node)

    # protect against qlogis(1)=Inf
    initial <- bound(initial, 0.005)

    submodel_data <- list(
      observed = observed,
      H = covariates_dt,
      initial = initial
    )
  } else {
    return(NULL)
  }
})

names(all_submodels) <- update_nodes

all_submodels
