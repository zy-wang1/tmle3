self <- updater

fold_number <- "full"

all_submodels <- self$generate_submodel_data(
  targeted_likelihood, tmle_task,
  fold_number
)
new_epsilons <- self$fit_submodels(all_submodels)

# update likelihoods
targeted_likelihood$update(new_epsilons, self$step_number, fold_number)  # now full lkd list is updated too

nothing <- suppressWarnings(lapply(self$tmle_params, function(tmle_param) {
  tmle_param$clever_covariates(tmle_task, fold_number
                               , update = T
  )  # this updates the covariates; it does not call full list of lkd
}))
nothing <- suppressWarnings(lapply(self$tmle_params, function(tmle_param) {
  tmle_param$estimates(tmle_task, fold_number
                       , update = T
  )  # this updates the D_list and results (est and ICs)
}))

if (fold_number != "full") {
  # update full fit likelihoods if we haven't already
  targeted_likelihood$update(new_epsilons, self$step_number, "full")  # now full lkd list is updated too
}
# increment step count
private$.step_number <- private$.step_number + 1
