#' Longitudinal Mediation via projection

#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Parameters
#' @keywords data
#'
#' @return \code{Param_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_param(Param_ATT, observed_likelihood, intervention_list, ..., outcome_node)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} corresponding to the observed likelihood
#'     }
#'     \item{\code{intervention_list_treatment}}{A list of objects inheriting from \code{\link{LF_base}}, representing the treatment intervention.
#'     }
#'     \item{\code{intervention_list_control}}{A list of objects inheriting from \code{\link{LF_base}}, representing the control intervention.
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{outcome_node}}{character, the name of the node that should be treated as the outcome
#'     }
#'     }
#'

#' @section Fields:
#' \describe{
#'     \item{\code{cf_likelihood_treatment}}{the counterfactual likelihood for the treatment
#'     }
#'     \item{\code{cf_likelihood_control}}{the counterfactual likelihood for the control
#'     }
#'     \item{\code{intervention_list_treatment}}{A list of objects inheriting from \code{\link{LF_base}}, representing the treatment intervention
#'     }
#'     \item{\code{intervention_list_control}}{A list of objects inheriting from \code{\link{LF_base}}, representing the control intervention
#'     }
#' }
#' @export
Param_middle_projection <- R6Class(
  classname = "Param_middle_projection",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list_treatment, intervention_list_control, outcome_node = "Y", static_likelihood = NULL) {
      temp_node_names <- names(observed_likelihood$training_task$npsem)
      loc_A <- grep("A", temp_node_names)
      loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
      loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
      if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
      tau <- last(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2]))

      all_nodes <- names(observed_likelihood$training_task$npsem)
      A_nodes <- grep("A", all_nodes, value = T)
      Z_nodes <- grep("Z", all_nodes, value = T)
      RLY_nodes <- grep("(R|L|Y).[1-9]$", all_nodes, value = T)

      # node_list <- lapply(observed_likelihood$training_task$npsem, function(x) x$variables)
      #
      # list_gamma_npsem_1 <- lapply(1:tau, function(t) {
      #   gamma_npsem(node_list, type = 1, j = t)
      # })
      # list_gamma_npsem_2 <- lapply(1:tau, function(t) {
      #   gamma_npsem(node_list, type = 2, j = t)
      # })
      #
      # list_gamma_task_1 <- lapply(list_gamma_npsem_1, function(each_npsem) {
      #   raw_data <- observed_likelihood$training_task$data %>% dplyr::select(-c(id, t)) %>% setDT
      #   tmle3_Task$new(raw_data, npsem = each_npsem)
      # })
      # list_gamma_task_2 <- lapply(list_gamma_npsem_2, function(each_npsem) {
      #   raw_data <- observed_likelihood$training_task$data %>% dplyr::select(-c(id, t)) %>% setDT
      #   tmle3_Task$new(raw_data, npsem = each_npsem)
      # })
      #
      # list_gamma_lkd_1 <- lapply(list_gamma_task_1, function(each_task) {
      #   middle_spec$make_initial_likelihood(
      #     each_task,
      #     learner_list
      #   )
      # })
      # list_gamma_lkd_2 <- lapply(list_gamma_task_2, function(each_task) {
      #   middle_spec$make_initial_likelihood(
      #     each_task,
      #     learner_list
      #   )
      # })
      # list_gamma_lkd_1[[1]]$get_likelihood(list_gamma_task_1[[1]], node = "A_1", fold_number = "full")
      #
      # list_gamma_lkd_1[[1]]$factor_list

      private$.static_likelihood <- static_likelihood




      private$.update_nodes <- c(Z_nodes, RLY_nodes)

      super$initialize(observed_likelihood, list(), outcome_node = outcome_node)
      private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
      private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)

      # todo: extend for stochastic
      private$.cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(observed_likelihood$training_task)[[1]]
      private$.cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(observed_likelihood$training_task)[[1]]

      # Train the gradient
      private$.gradient <- Gradient$new(observed_likelihood,
                                        ipw_args = list(cf_likelihood_treatment = self$cf_likelihood_treatment,
                                                        cf_likelihood_control = self$cf_likelihood_control,
                                                        intervention_list_treatment = self$intervention_list_treatment,
                                                        intervention_list_control = self$intervention_list_control,
                                                        cf_task_treatment = self$cf_task_treatment,
                                                        cf_task_control = self$cf_task_control,
                                                        static_likelihood = self$static_likelihood
                                        ),
                                        projection_task_generator = gradient_generator_middle,
                                        target_nodes =  self$update_nodes)


      if(inherits(observed_likelihood, "Targeted_Likelihood")){
        fold_number <- observed_likelihood$updater$update_fold
      } else {
        fold_number <- "full"
      }

      setattr(self$observed_likelihood, "target_nodes", self$update_nodes)
      self$observed_likelihood$get_likelihoods(self$observed_likelihood$training_task)
      for (node in self$update_nodes) {
        temp_long_task <- private$.gradient$expand_task(observed_likelihood$training_task, node)
        self$observed_likelihood$get_likelihood(temp_long_task, node, fold_number)
        self$observed_likelihood$get_likelihood(observed_likelihood$training_task, node, fold_number)
        # private$.gradient$expand_task(private$.cf_task_treatment, node)
        # private$.gradient$expand_task(private$.cf_task_control, node)
      }

      private$.gradient$train_projections(self$observed_likelihood$training_task, fold_number = fold_number)

      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))
      obs_variable_names <- colnames(obs_data)
      list_all_predicted_lkd <- lapply(1:length(temp_node_names), function(loc_node) {
        if (loc_node > 1) {
          # currently only support univariate node for t>0
          current_variable <- tmle_task$npsem[[loc_node]]$variables
          temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)])  # all possible inputs
          temp_task <- tmle3_Task$new(temp_input, tmle_task$npsem[1:loc_node])
          temp_target_node <- intersect(self$update_nodes, temp_node_names[loc_node])
          if (length(temp_target_node) == 1) {
            # for each short task, only the last node (if it is an update_node) needs to be updated
            setattr(temp_task, "target_nodes", temp_target_node)
            for (node in attr(temp_task, "target_nodes")) {
              temp_long_task <- private$.gradient$expand_task(temp_task, node)
              self$observed_likelihood$get_likelihood(temp_long_task, node, fold_number)
            }
            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number)  # corresponding outputs
          } else {
            # A nodes won't get updated
            temp_output <- self$static_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number)  # corresponding outputs
          }
          data.frame(temp_input, output = temp_output) %>% return
        }
      })


    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", node = NULL) {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      update_nodes <- intersect(self$update_nodes, attr(tmle_task, "target_nodes"))
      if(!is.null(node)){
        update_nodes <- c(node)
      }
      islong = F
      if(is.null(update_nodes)){
        update_nodes <- self$update_nodes
      } else {
        islong= T
      }
      EICs <- lapply(update_nodes, function(node){
        return(self$gradient$compute_component(tmle_task, node, fold_number = fold_number)$EIC)
      })

      names(EICs) <- update_nodes

      EICs[[length(EICs) + 1]] <- do.call(cbind, EICs)
      colnames(EICs[[length(EICs)]]) <- update_nodes
      EICs[[length(EICs)]] <- as.data.frame(EICs[[length(EICs)]])

      names(EICs)[length(EICs)] <- "IC"

      return(EICs)
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }

      intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))

      # clever_covariates happen here (for this param) only, but this is repeated computation
      EIC <- (do.call(cbind, self$clever_covariates(tmle_task, fold_number)$IC))

      #TODO need to montecarlo simulate from likleihood to eval parameter.

      # todo: make sure we support updating these params
      # pA <- self$observed_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      # cf_pA_treatment <- self$cf_likelihood_treatment$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      # cf_pA_control <- self$cf_likelihood_control$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      #
      # # todo: extend for stochastic
      # cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
      # cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]
      #
      # Y <- tmle_task$get_tmle_node(self$outcome_node, impute_censoring = TRUE)
      #
      # EY <- self$observed_likelihood$get_likelihood(tmle_task, self$outcome_node, fold_number)
      # EY1 <- self$observed_likelihood$get_likelihood(cf_task_treatment, self$outcome_node, fold_number)
      # EY0 <- self$observed_likelihood$get_likelihood(cf_task_control, self$outcome_node, fold_number)
      #
      # psi <- mean(EY1 - EY0)
      #
      # IC <- EIC + (EY1 - EY0) - psi

      temp_node_names <- names(self$observed_likelihood$training_task$npsem)
      loc_A <- grep("A", temp_node_names)
      loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
      loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
      if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
      tau <- last(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2]))

      obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))
      obs_variable_names <- colnames(obs_data)
      list_all_predicted_lkd <- lapply(1:length(temp_node_names), function(loc_node) {
        if (loc_node > 1) {
          # currently only support univariate node for t>0
          current_variable <- tmle_task$npsem[[loc_node]]$variables
          temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)])  # all possible inputs
          temp_task <- tmle3_Task$new(temp_input, tmle_task$npsem[1:loc_node])
          temp_target_node <- intersect(self$update_nodes, temp_node_names[loc_node])
          if (length(temp_target_node) == 1) {
            # for each short task, only the last node (if it is an update_node) needs to be updated
            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number)  # corresponding outputs
          } else {
            # A nodes won't get updated
            temp_output <- self$static_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number)  # corresponding outputs
          }
          data.frame(temp_input, output = temp_output) %>% return
        }
      })

      intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
      intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
      intervention_levels_treat <- map_dbl(self$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
      intervention_levels_control <- map_dbl(self$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
      # nodes to integrate out in the target identification
      # only support univaraite node for now; assume treatment level is one
      all_possible_RZLY_1 <- expand_values(obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ),
                                           rule_variables = c(intervention_variables,
                                                              last(obs_variable_names)),
                                           rule_values = c(intervention_levels_treat, 1))
      all_possible_RZLY_0 <- expand_values(obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ),
                                           rule_variables = c(intervention_variables,
                                                              last(obs_variable_names)),
                                           rule_values = c(intervention_levels_control, 1))

      # for each observed L_0 vector, generate all needed combinations, one version for A = 1, one version for A = 0
      unique_L0 <- obs_data[, tmle_task$npsem[[1]]$variables] %>% unique
      library_L0 <- data.frame(unique_L0, output =
                                 map_dbl(1:nrow(unique_L0), function(which_row) {
                                   temp_all_comb_0 <- cbind(unique_L0[which_row, ], all_possible_RZLY_0)
                                   temp_all_comb_1 <- cbind(unique_L0[which_row, ], all_possible_RZLY_1)
                                   # for all non-A, non-0 variables, calculate the variable by rule
                                   # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
                                   # note that list_all_predicted_lkd is ordered by node
                                   temp_list_0 <- lapply(loc_Z,
                                                         function(each_t) {
                                                           left_join(temp_all_comb_0, list_all_predicted_lkd[[each_t]])$output
                                                         })
                                   temp_list_1 <- lapply(loc_RLY,
                                                         function(each_t) {
                                                           left_join(temp_all_comb_1, list_all_predicted_lkd[[each_t]])$output
                                                         })
                                   temp_list <- c(temp_list_0, temp_list_1)
                                   pmap_dbl(temp_list, prod) %>% sum %>% return
                                 })
      )
      vec_est <- left_join(obs_data[, tmle_task$npsem[[1]]$variables], library_L0)$output
      psi <- mean(vec_est)

      IC <- rowSums(EIC)
      result <- list(psi =
                       psi
                     # list_all_predicted_lkd
                     ,
                     IC = IC, EIC = colMeans(EIC)
                     # , full_EIC = EIC
      )
      return(result)
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf("ATE[%s_{%s}-%s_{%s}]", self$outcome_node, self$cf_likelihood_treatment$name, self$outcome_node, self$cf_likelihood_control$name)
      return(param_form)
    },
    cf_likelihood_treatment = function() {
      return(private$.cf_likelihood_treatment)
    },
    cf_likelihood_control = function() {
      return(private$.cf_likelihood_control)
    },
    cf_task_treatment = function() {
      return(private$.cf_task_treatment)
    },
    cf_task_control = function() {
      return(private$.cf_task_control)
    },
    intervention_list_treatment = function() {
      return(self$cf_likelihood_treatment$intervention_list)
    },
    intervention_list_control = function() {
      return(self$cf_likelihood_control$intervention_list)
    },
    update_nodes = function() {
      return(c(private$.update_nodes))
    },
    gradient = function(){
      private$.gradient
    },
    static_likelihood = function(){
      private$.static_likelihood
    }
  ),
  private = list(
    .type = "ATE",
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL,
    .cf_task_treatment = NULL,
    .cf_task_control = NULL,
    .supports_outcome_censoring = FALSE,
    .gradient = NULL,
    .submodel_type_supported = c("EIC"),
    .update_nodes = NULL,
    .static_likelihood = NULL
  )
)

