self <- updater

update_fold <- self$update_fold
maxit <- self$maxit

# seed current estimates
current_estimates <- lapply(self$tmle_params, function(tmle_param) {
  tmle_param$estimates(tmle_task, update_fold)
})

if(FALSE) {
  clever_covariates <- lapply(self$tmle_params, function(tmle_param) {
    tmle_param$clever_covariates(tmle_task, update_fold)})
  IC <- lapply(clever_covariates, `[[`, "IC")
  if(!is.null(IC[[1]])){
    n <- length(unique(tmle_task$id))
    IC_vars <- lapply(IC, function(IC) {
      out <- lapply(self$update_nodes, function(node) {
        weights <- tmle_task$get_regression_task(node)$weights
        apply(IC[[node]] * weights,2, function(v) {var(rowSums(matrix(v, nrow = n, byrow = T)))})
      } )
      names(out) <- self$update_nodes
      return(out)
    })
    private$.initial_variances <- IC_vars


  } else {
    n <- length(unique(tmle_task$id))
    IC <- lapply(private$.current_estimates, `[[`, "IC")
    IC_vars <- lapply(IC, function(IC) {
      weights <- tmle_task$get_regression_task(node)$weights
      IC_var <- apply(IC[[node]] * weights,2, function(v) {var(rowSums(matrix(v, nrow = n, byrow = T)))})
      IC_var <- lapply(self$update_nodes, function(node) {IC_var})
      names(IC_var) <- self$update_nodes
      return(IC_var)
    })
    private$.initial_variances <- IC_vars
  }
}


#private$.initial_variances <- lapply(private$.current_estimates, `[[`, "var_comps")

for (steps in seq_len(maxit)) {
  self$update_step(targeted_likelihood, tmle_task, update_fold)

  # update estimates based on updated likelihood
  current_estimates <- lapply(self$tmle_params, function(tmle_param) {
    tmle_param$estimates(tmle_task, update_fold)
  })

  if (self$check_convergence(tmle_task, update_fold)) {
    break
  }

  if (self$use_best) {
    self$update_best(likelihood)
  }
}

if (self$use_best) {
  self$update_best(likelihood)
  likelihood$cache$set_best()
}
