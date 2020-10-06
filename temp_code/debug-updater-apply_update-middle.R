self <- updater
task <- targeted_likelihood$cache$tasks[[2]]
all_epsilon <- new_epsilons


update_nodes <- self$update_nodes

# all_submodel_data <- all_submodels
# get submodel data for all nodes
all_submodel_data <- self$generate_submodel_data(
  targeted_likelihood, task,
  fold_number
)

# needed to transform between lkd and probs
observed_values <- lapply(update_nodes, task$get_tmle_node)
names(observed_values) <- update_nodes



# apply update to all nodes
updated_likelihoods <- lapply(update_nodes, function(update_node) {
  if (is.null(all_submodel_data[[update_node]])) return(NULL) else {
    submodel_data <- all_submodel_data[[update_node]]
    epsilon <- all_epsilon[[update_node]]
    updated_likelihood <- self$apply_submodel(submodel_data, epsilon)

    # we updated the predicted probabilities, or the non-zero lt likelihoods
    # for lt=0 they are just 1-p or 1- sum of other p
    updated_likelihood <- ifelse(observed_values[[update_node]] == 1, updated_likelihood, 1 - updated_likelihood)

    # un-scale to handle bounded continuous
    updated_likelihood <- task$unscale(
      updated_likelihood,
      update_node
    )
  }
})
names(updated_likelihoods) <- update_nodes
