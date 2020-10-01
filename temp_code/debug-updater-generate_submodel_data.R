self <- tlik$updater
for_fitting <- T
fold_number <- "full"
drop_censored <- T

if(!(inherits(tlik, "Targeted_Likelihood"))) {
  submodel_type <- "logistic"
} else {
  submodel_type <- tlik$submodel_type(update_node)
}

submodel_info <- submodel_spec(submodel_type)
# TODO: change clever covariates to allow only calculating some nodes

tmle_param <- self$tmle_params[[1]]

node <- update_node
tmle_params[[1]]$gradient$compute_component(tmle_task, node, fold_number = fold_number)$EIC
tlik$cache$tasks %>% length
tlik$cache$tasks %>% lapply(function(x) x$data$A_1 %>% table)




sl3:::call_with_args(tmle_param$clever_covariates, args)

clever_covariates <- lapply(self$tmle_params, function(tmle_param) {
  # Assert that it supports the submodel type
  tmle_param$supports_submodel_type(submodel_type, update_node)
  #formal_args <- names(formals(tmle_param$clever_covariates))

  # For backwards compatibility:
  # In future, clever covariate functions should accept a "node" and "submodel_type" argument.
  args <- list(for_fitting = for_fitting, submodel_type = submodel_type, fold_number = fold_number, tmle_task = tmle_task
               ,node = update_node)
  return(sl3:::call_with_args(tmle_param$clever_covariates, args))
})



tlik$cache$tasks %>% length
tlik$cache$tasks %>% lapply(function(x) x$data$A_1 %>% table)







node_covariates <- lapply(clever_covariates, `[[`, update_node)
# Get EDs if present. Only for training task
if(self$one_dimensional & for_fitting) {
  IC <- lapply(clever_covariates, `[[`, "IC")
  IC <- do.call(cbind, lapply(IC, `[[`, update_node) )
  if(is.null(IC)) {
    IC <- lapply(private$.current_estimates, `[[`, "IC")
    IC <- do.call(cbind, IC)
  }
}
covariates_dt <- do.call(cbind, node_covariates)

observed <- tmle_task$get_tmle_node(update_node)
initial <- tlik$get_likelihood(
  tmle_task, update_node,
  fold_number)
# scale observed and predicted values for bounded continuous
observed <- tmle_task$scale(observed, update_node)
initial <- tmle_task$scale(initial, update_node)

tlik$cache$tasks %>% length
tlik$cache$tasks %>% lapply(function(x) x$data$A_1 %>% table)


# protect against qlogis(1)=Inf
initial <- bound(initial, 0.005)
weights <- tmle_task$get_regression_task(update_node)$weights
n <- length(unique(tmle_task$id))
if(self$one_dimensional & for_fitting){
  # This computes (possibly weighted) ED and handles long case
  ED <- colSums(IC * weights)/n #apply(IC , 2, function(v) {sum(as.vector(matrix(v, nrow = n, byrow = T)*weights))})/length(weights)
} else {
  ED <- NULL
}





if(length(observed) != length(initial)) {
  ratio <- length(initial) / length(observed)
  if(ratio%%1 == 0){
    warning("Observed and initial length do not match but are multiples of each other. Recycling values...")
    observed <- rep(observed, ratio)
  }
}

if(length(weights) != length(initial)) {
  ratio <- length(initial) / length(weights)
  if(ratio%%1 == 0){
    # This is for likelihood factors that output long_format predictions that dont match nrow of input task
    warning("Weights and initial length do not match but are multiples of each other. Recycling values...")
    weights <- rep(weights, ratio)
  }
}

submodel_data <- list(
  observed = observed,
  H = covariates_dt,
  initial = initial,
  submodel_info = submodel_info,
  ED = ED,
  update_node = update_node,
  weights = weights
)



if (drop_censored) {
  censoring_node <- tmle_task$npsem[[update_node]]$censoring_node$name
  if (!is.null(censoring_node)) {
    observed_node <- tmle_task$get_tmle_node(censoring_node)
    subset <- which(observed_node == 1)
    subset <- intersect(subset, which(tmle_task$get_tmle_node(update_node, compute_risk_set = T)[, at_risk] == 1))
    submodel_data <- list(
      observed = submodel_data$observed[subset],
      H = submodel_data$H[subset, , drop = FALSE],
      initial = submodel_data$initial[subset],
      submodel_info = submodel_info,
      ED = ED,
      update_node = update_node,
      weights = submodel_data$weights[subset]
    )
  }
}




tlik$cache$tasks %>% length
tlik$cache$tasks %>% lapply(function(x) x$data$A_1 %>% table)
