#' Longitudinal Mediation Targets
#'
#' Parameter definition for longitudinal mediation
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
#'   \code{define_param(Param_mediation, observed_likelihood, intervention_list, ..., outcome_node)}
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
Param_medidation <- R6Class(
  classname = "Param_mediation",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list_treatment, intervention_list_control) {
      # outcome_node is used to check self$supports_outcome_censoring; not checked for now
      super$initialize(observed_likelihood, list())
      private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
      private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)
      observed_likelihood$get_likelihoods(observed_likelihood$training_task)
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", update = F, node = NULL, submodel_type = "EIC") {
      if (is.null(tmle_task)) {  # calculate for obs data task if not specified
        tmle_task <- self$observed_likelihood$training_task
      }
      intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))


      if (fold_number == "full") {  # tmle
        list_newH <- private$.list_newH_EIC
      } else if (fold_number == "validation") {  # cvtmle
        list_newH <- private$.list_newH_EIC_val
      }  # load cached obs task clever covariates in case its for convergence check
      if (!is.null(list_newH) & update == F & identical(tmle_task, self$observed_likelihood$training_task)) {  # for faster convergence check
        if (!is.null(node)) {  # return partial list of covariates if requested
          return(list_newH[node])
        } else {
          # list_newH <- c(list_newH, self$clever_covariates(fold_number = fold_number, submodel_type = "EIC"))  # append the IC to clever covariates
          # names(list_newH)[length(list_newH)] <- "IC"
          return(list_newH)
        }
      } else {  # note submodel_type; only calculate when i) no cached newH, ii) forced to update after tlik is updated; or iii) not obs task
        rm(list_newH)

        # load full_p list first
        full_task <- self$observed_likelihood$training_task
        full_node_names <- names(full_task$npsem)
        full_data <- full_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))  # exactly the obs data
        full_variable_names <- colnames(full_data)
        list_all_predicted_lkd <- lapply(1:length(full_node_names), function(loc_node) {
          if (loc_node > 1) {
            # currently only support univariate node for t>0
            current_variable <- full_task$npsem[[loc_node]]$variables
            temp_input <- expand_values(variables = full_variable_names[1:which(full_variable_names == current_variable)])  # all possible inputs
            temp_task <- tmle3_Task$new(temp_input, full_task$npsem[1:loc_node])
            temp_target_node <- intersect(self$update_nodes, full_node_names[loc_node])
            if (length(temp_target_node) == 1) {
              setattr(temp_task, "target_nodes", full_node_names[loc_node])
              temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
            } else {
              setattr(temp_task, "target_nodes", "no_update")
              temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
            }
            data.frame(temp_input, output = temp_output) %>% return
          }
        })
        names(list_all_predicted_lkd) <- full_node_names

        if (all(tmle_task$nrow == self$observed_likelihood$training_task$nrow,
                identical(tmle_task$data[[1]], self$observed_likelihood$training_task$data[[1]])
        )) {  # for cf or obs tasks
          # ZW todo: extend for dynamic treatments
          cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
          cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]

          temp_node_names <- names(tmle_task$npsem)
          loc_A <- grep("A", temp_node_names)
          loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
          loc_RLY <- which(sapply(temp_node_names, function(s) !(strsplit(s, "_")[[1]][1] %in% c("A", "Z")) & strsplit(s, "_")[[1]][2] != 0))

          obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))  # note this is compatible if tmle_task is a cf task
          obs_variable_names <- colnames(obs_data)
          # ZW todo: to handle long format and wide format

          intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
          intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
          intervention_levels_treat <- map_dbl(self$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
          intervention_levels_control <- map_dbl(self$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
          names(intervention_levels_treat) <- names(self$intervention_list_treatment)
          names(intervention_levels_control) <- names(self$intervention_list_control)

          list_H <- get_obs_H(tmle_task, obs_data, current_likelihood = self$observed_likelihood,
                              cf_task_treatment, cf_task_control,
                              intervention_variables, intervention_levels_treat, intervention_levels_control,
                              fold_number)
          list_Q_1 <- get_obs_Q(tmle_task, obs_data, list_H,
                                intervention_variables, intervention_levels_treat, intervention_levels_control,
                                list_all_predicted_lkd,  # val version decided above for fold_number == "validation"
                                lt = 1)
          list_Q_0 <- get_obs_Q(tmle_task, obs_data, list_H,
                                intervention_variables, intervention_levels_treat, intervention_levels_control,
                                list_all_predicted_lkd,
                                lt = 0)

          list_newH <- list()
          for (loc_node in 1:length(list_H)) {
            if(!is.null(list_H[[loc_node]])) {  # for density update, change the sign of some clever covariates
              temp_vec <- list_H[[loc_node]] * (list_Q_1[[loc_node]] - list_Q_0[[loc_node]])
              list_newH[[loc_node]] <- temp_vec
            }
          }
          names(list_newH) <- temp_node_names

          # calculate EIC components
          list_D <- list()
          for (loc_node in 1:length(list_newH)) {
            if(!is.null(list_newH[[loc_node]])) {
              # ZW todo: for discretized variables
              current_ind <- (obs_data[[tmle_task$npsem[[loc_node]]$variables]] == 1)*1
              temp_vec <- list_newH[[loc_node]] %>% as.vector
              temp_p <- self$observed_likelihood$get_likelihoods(tmle_task, temp_node_names[loc_node], fold_number)
              # if (loc_node %in% loc_Z) temp_p <- self$observed_likelihood$get_likelihoods(cf_task_control, temp_node_names[loc_node], fold_number) else
              #   temp_p <- self$observed_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[loc_node], fold_number)
              temp_p <- ifelse(current_ind == 1, temp_p, 1 - temp_p)  # transform density to conditional mean
              list_D[[loc_node]] <- (current_ind - temp_p) * temp_vec
            }
          }
          # list_D[[1]] <- vec_est - psi
          names(list_D) <- names(list_newH)

          if (submodel_type != "logistic") {  # get EIC as clever covariates for EIC submodels
            list_newH <- list_D  # EIC is the clever covariates for EIC submodels
          }

          list_newH[[length(list_newH) + 1]] <- do.call(cbind, list_D)
          names(list_newH)[length(list_newH)] <- "IC"  # to use in by dimension convergence

          if (identical(tmle_task, self$observed_likelihood$training_task)) {  # cache for obs task
            if (fold_number == "full") {
              if (submodel_type == "logistic") private$.list_newH <- list_newH else private$.list_newH_EIC <- list_newH
            } else if (fold_number == "validation") {
              if (submodel_type == "logistic") private$.list_newH_val <- list_newH else private$.list_newH_EIC_val <- list_newH
            }
          }

          if (!is.null(node)) {  # return partial list of covariates if requested
            return(list_newH[node])
          } else return(list_newH)
        } else {  # for library tasks; it's only needed in tlik updates, with single node
          if (is.null(node)) stop("Please specify single update node for library tasks")

          tmle_task_backup <- tmle_task
          tmle_task <- self$observed_likelihood$training_task  # let tmle_task be obs task when calculating for library tasks
          loc_node <- which(names(tmle_task$npsem) == node)
          obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))
          obs_variable_names <- names(obs_data)
          intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
          intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
          intervention_levels_treat <- map_dbl(self$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
          intervention_levels_control <- map_dbl(self$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
          names(intervention_levels_treat) <- names(self$intervention_list_treatment)
          names(intervention_levels_control) <- names(self$intervention_list_control)

          current_newH <- get_current_newH(loc_node,
                                           tmle_task, obs_data,
                                           intervention_variables, intervention_levels_treat, intervention_levels_control,
                                           list_all_predicted_lkd  # this is decided above by fold_number
          )  # this is what we need for logistic submodel

          if (submodel_type != "logistic") {  # for EIC, add needed X - E(X) term
            observed <- tmle_task_backup$get_tmle_node(loc_node)  # get X
            EX <- list_all_predicted_lkd[[node]]$output  # p(X)
            EX <- ifelse(observed == 1, EX, 1-EX)  # transfrom p(X) to E(X) as its needed in D*
            current_newH <- (observed - EX) * current_newH  # it can be raw observed X, since it's a density update
          }

          current_newH <- list(current_newH)
          names(current_newH) <- node

          return(current_newH)
        }

      }
    },
    estimates = function(tmle_task = NULL, fold_number = "full", update = F) {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))

      # todo: extend for stochastic
      cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
      cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]

      temp_node_names <- names(tmle_task$npsem)
      loc_A <- grep("A", temp_node_names)
      loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
      loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
      if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)

      obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))
      obs_variable_names <- colnames(obs_data)
      # ZW todo: to handle long format and wide format

      if (fold_number == "full") {
        list_D <- private$.list_D
        result <- private$.result
      } else if (fold_number == "validation") {
        list_D <- private$.list_D_val
        result <- private$.result_val
      }

      # load full_p list first
      full_task <- self$observed_likelihood$training_task
      full_node_names <- names(full_task$npsem)
      full_data <- full_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))  # exactly the obs data
      full_variable_names <- colnames(full_data)
      list_all_predicted_lkd <- lapply(1:length(full_node_names), function(loc_node) {
        if (loc_node > 1) {
          # currently only support univariate node for t>0
          current_variable <- full_task$npsem[[loc_node]]$variables
          temp_input <- expand_values(variables = full_variable_names[1:which(full_variable_names == current_variable)])  # all possible inputs
          temp_task <- tmle3_Task$new(temp_input, full_task$npsem[1:loc_node])
          temp_target_node <- intersect(self$update_nodes, full_node_names[loc_node])
          if (length(temp_target_node) == 1) {
            setattr(temp_task, "target_nodes", full_node_names[loc_node])
            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
          } else {
            setattr(temp_task, "target_nodes", "no_update")
            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
          }
          data.frame(temp_input, output = temp_output) %>% return
        }
      })
      names(list_all_predicted_lkd) <- full_node_names

      # make sure we force update this after each updating step; this helps speed up updater$check_convergence
      if (!is.null(list_D) & !is.null(result) & update == F) {
        return(result)
      } else {
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
        # substitution estimator
        vec_est <- left_join(obs_data[, tmle_task$npsem[[1]]$variables], library_L0)$output
        psi <- mean(vec_est)

        list_D <- self$clever_covariates(tmle_task, fold_number, submodel_type = "EIC")
        list_D[[1]] <- vec_est - psi
        list_D <- list_D[-which(names(list_D) == "IC")]

        vec_D <- list_D %>% compact %>% pmap_dbl(sum)
        IC <- vec_D

        result <- list(psi = psi, IC = IC
                       # , full_IC = list_D
        )

        # these are cached; unless likelihood is updated, or we force it to update, they shouldn't be changed
        if (fold_number == "full") {
          private$.list_D <- list_D
          private$.result <- result
        } else if (fold_number == "validation") {
          private$.list_D_val <- list_D
          private$.result_val <- result
        }

        return(result)
      }
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
    intervention_list_treatment = function() {
      return(self$cf_likelihood_treatment$intervention_list)
    },
    intervention_list_control = function() {
      return(self$cf_likelihood_control$intervention_list)
    },
    list_D = function() {
      return(private$.list_D)
    },
    list_D_val = function() {
      return(private$.list_D_val)
    },
    update_nodes = function() {
      # if (is.null(tmle_task)) {
      tmle_task <- self$observed_likelihood$training_task
      # }
      temp_node_names <- names(tmle_task$npsem)
      loc_A <- grep("A", temp_node_names)
      if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
      nodes_to_update <- temp_node_names[if_not_0 & !((1:length(temp_node_names)) %in% loc_A)]
      # nodes_to_update <- nodes_to_update[-length(nodes_to_update)]
      return(nodes_to_update)
    },
    list_newH = function() {
      return(private$.list_newH)
    },
    list_newH_val = function() {
      return(private$.list_newH_val)
    },
    list_newH_EIC = function() {
      return(private$.list_newH_EIC)
    },
    list_newH_EIC_val = function() {
      return(private$.list_newH_EIC_val)
    }
  ),
  private = list(
    .type = "mediation",
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL,
    .list_newH = NULL,  # the clever covariates for logistic TMLE
    .list_newH_val = NULL,
    .list_newH_EIC = NULL,  # the clever covariates as the EIC
    .list_newH_EIC_val = NULL,
    .list_D = NULL,
    .list_D_val = NULL,
    .result = NULL,
    .result_val = NULL,
    .submodel_type_supported = c("EIC")
  )
)











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
Param_mediation_survival <- R6Class(
  classname = "Param_mediation",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list_treatment, intervention_list_control, outcome_node = "Y") {
      super$initialize(observed_likelihood, list(), outcome_node = outcome_node)
      private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
      private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)
      observed_likelihood$get_likelihoods(observed_likelihood$training_task)
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", update = F, node = NULL, submodel_type = "logistic") {
      if (is.null(tmle_task)) {  # calculate for obs data task if not specified
        tmle_task <- self$observed_likelihood$training_task
      }
      intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))


      if (fold_number == "full") {  # tmle
        list_newH <- ifelse_vec(submodel_type == "logistic", private$.list_newH, private$.list_newH_EIC)
      } else if (fold_number == "validation") {  # cvtmle
        list_newH <- ifelse_vec(submodel_type == "logistic", private$.list_newH_val, private$.list_newH_EIC_val)
      }  # load cached obs task clever covariates in case its for convergence check
      if (!is.null(list_newH) & update == F & identical(tmle_task, self$observed_likelihood$training_task)) {  # for faster convergence check
        if (!is.null(node)) {  # return partial list of covariates if requested
          return(list_newH[node])
        } else {
          # list_newH <- c(list_newH, self$clever_covariates(fold_number = fold_number, submodel_type = "EIC"))  # append the IC to clever covariates
          # names(list_newH)[length(list_newH)] <- "IC"
          return(list_newH)
        }
      } else {  # note submodel_type; only calculate when i) no cached newH, ii) forced to update after tlik is updated; or iii) not obs task
        rm(list_newH)

        # load full_p list first
        full_task <- self$observed_likelihood$training_task
        full_node_names <- names(full_task$npsem)
        full_data <- full_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))  # exactly the obs data
        full_variable_names <- colnames(full_data)
        list_all_predicted_lkd <- lapply(1:length(full_node_names), function(loc_node) {
          if (loc_node > 1) {
            # currently only support univariate node for t>0
            current_variable <- full_task$npsem[[loc_node]]$variables
            temp_input <- expand_values(variables = full_variable_names[1:which(full_variable_names == current_variable)])  # all possible inputs
            temp_task <- tmle3_Task$new(temp_input, full_task$npsem[1:loc_node])
            temp_target_node <- intersect(self$update_nodes, full_node_names[loc_node])
            if (length(temp_target_node) == 1) {
              setattr(temp_task, "target_nodes", full_node_names[loc_node])
              temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
            } else {
              setattr(temp_task, "target_nodes", "no_update")
              temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
            }
            data.frame(temp_input, output = temp_output) %>% return
          }
        })
        names(list_all_predicted_lkd) <- full_node_names

        if (all(tmle_task$nrow == self$observed_likelihood$training_task$nrow,
                identical(tmle_task$data[[1]], self$observed_likelihood$training_task$data[[1]])
        )) {  # for cf or obs tasks
          # ZW todo: extend for dynamic treatments
          cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
          cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]

          temp_node_names <- names(tmle_task$npsem)
          loc_A <- grep("A", temp_node_names)
          loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
          loc_RLY <- which(sapply(temp_node_names, function(s) !(strsplit(s, "_")[[1]][1] %in% c("A", "Z")) & strsplit(s, "_")[[1]][2] != 0))

          obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))  # note this is compatible if tmle_task is a cf task
          obs_variable_names <- colnames(obs_data)
          # ZW todo: to handle long format and wide format

          intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
          intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
          intervention_levels_treat <- map_dbl(self$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
          intervention_levels_control <- map_dbl(self$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
          names(intervention_levels_treat) <- names(self$intervention_list_treatment)
          names(intervention_levels_control) <- names(self$intervention_list_control)

          list_H <- get_obs_H(tmle_task, obs_data, current_likelihood = self$observed_likelihood,
                              cf_task_treatment, cf_task_control,
                              intervention_variables, intervention_levels_treat, intervention_levels_control,
                              fold_number)
          list_Q_1 <- get_obs_Q(tmle_task, obs_data, list_H,
                                intervention_variables, intervention_levels_treat, intervention_levels_control,
                                list_all_predicted_lkd,  # val version decided above for fold_number == "validation"
                                lt = 1)
          list_Q_0 <- get_obs_Q(tmle_task, obs_data, list_H,
                                intervention_variables, intervention_levels_treat, intervention_levels_control,
                                list_all_predicted_lkd,
                                lt = 0)

          list_newH <- list()
          for (loc_node in 1:length(list_H)) {
            if(!is.null(list_H[[loc_node]])) {  # for density update, change the sign of some clever covariates
              temp_vec <- list_H[[loc_node]] * (list_Q_1[[loc_node]] - list_Q_0[[loc_node]])
              list_newH[[loc_node]] <- temp_vec
            }
          }
          names(list_newH) <- temp_node_names

          # calculate EIC components
          list_D <- list()
          for (loc_node in 1:length(list_newH)) {
            if(!is.null(list_newH[[loc_node]])) {
              # ZW todo: for discretized variables
              current_ind <- (obs_data[[tmle_task$npsem[[loc_node]]$variables]] == 1)*1
              temp_vec <- list_newH[[loc_node]] %>% as.vector
              temp_p <- self$observed_likelihood$get_likelihoods(tmle_task, temp_node_names[loc_node], fold_number)
              # if (loc_node %in% loc_Z) temp_p <- self$observed_likelihood$get_likelihoods(cf_task_control, temp_node_names[loc_node], fold_number) else
              #   temp_p <- self$observed_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[loc_node], fold_number)
              temp_p <- ifelse(current_ind == 1, temp_p, 1 - temp_p)  # transform density to conditional mean
              list_D[[loc_node]] <- (current_ind - temp_p) * temp_vec
            }
          }
          # list_D[[1]] <- vec_est - psi
          names(list_D) <- names(list_newH)

          if (submodel_type != "logistic") {  # get EIC as clever covariates for EIC submodels
            list_newH <- list_D  # EIC is the clever covariates for EIC submodels
          }

          list_newH[[length(list_newH) + 1]] <- do.call(cbind, list_D)
          names(list_newH)[length(list_newH)] <- "IC"  # to use in by dimension convergence

          if (identical(tmle_task, self$observed_likelihood$training_task)) {  # cache for obs task
            if (fold_number == "full") {
              if (submodel_type == "logistic") private$.list_newH <- list_newH else private$.list_newH_EIC <- list_newH
            } else if (fold_number == "validation") {
              if (submodel_type == "logistic") private$.list_newH_val <- list_newH else private$.list_newH_EIC_val <- list_newH
            }
          }

          if (!is.null(node)) {  # return partial list of covariates if requested
            return(list_newH[node])
          } else return(list_newH)
        } else {  # for library tasks; it's only needed in tlik updates, with single node
          if (is.null(node)) stop("Please specify single update node for library tasks")

          tmle_task_backup <- tmle_task
          tmle_task <- self$observed_likelihood$training_task  # let tmle_task be obs task when calculating for library tasks
          loc_node <- which(names(tmle_task$npsem) == node)
          obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))
          obs_variable_names <- names(obs_data)
          intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
          intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
          intervention_levels_treat <- map_dbl(self$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
          intervention_levels_control <- map_dbl(self$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
          names(intervention_levels_treat) <- names(self$intervention_list_treatment)
          names(intervention_levels_control) <- names(self$intervention_list_control)

          current_newH <- get_current_newH(loc_node,
                                           tmle_task, obs_data,
                                           intervention_variables, intervention_levels_treat, intervention_levels_control,
                                           list_all_predicted_lkd  # this is decided above by fold_number
          )  # this is what we need for logistic submodel

          if (submodel_type != "logistic") {  # for EIC, add needed X - E(X) term
            observed <- tmle_task_backup$get_tmle_node(loc_node)  # get X
            EX <- list_all_predicted_lkd[[node]]$output  # p(X)
            EX <- ifelse(observed == 1, EX, 1-EX)  # transfrom p(X) to E(X) as its needed in D*
            current_newH <- (observed - EX) * current_newH  # it can be raw observed X, since it's a density update
          }

          current_newH <- list(current_newH)
          names(current_newH) <- node

          return(current_newH)
        }

      }
    },
    estimates = function(tmle_task = NULL, fold_number = "full", update = F) {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))

      # todo: extend for stochastic
      cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
      cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]

      temp_node_names <- names(tmle_task$npsem)
      loc_A <- grep("A", temp_node_names)
      loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
      loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
      if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)

      obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))
      obs_variable_names <- colnames(obs_data)
      # ZW todo: to handle long format and wide format

      if (fold_number == "full") {
        list_D <- private$.list_D
        result <- private$.result
      } else if (fold_number == "validation") {
        list_D <- private$.list_D_val
        result <- private$.result_val
      }

      # load full_p list first
      full_task <- self$observed_likelihood$training_task
      full_node_names <- names(full_task$npsem)
      full_data <- full_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))  # exactly the obs data
      full_variable_names <- colnames(full_data)
      list_all_predicted_lkd <- lapply(1:length(full_node_names), function(loc_node) {
        if (loc_node > 1) {
          # currently only support univariate node for t>0
          current_variable <- full_task$npsem[[loc_node]]$variables
          temp_input <- expand_values(variables = full_variable_names[1:which(full_variable_names == current_variable)])  # all possible inputs
          temp_task <- tmle3_Task$new(temp_input, full_task$npsem[1:loc_node])
          temp_target_node <- intersect(self$update_nodes, full_node_names[loc_node])
          if (length(temp_target_node) == 1) {
            setattr(temp_task, "target_nodes", full_node_names[loc_node])
            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
          } else {
            setattr(temp_task, "target_nodes", "no_update")
            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = full_node_names[loc_node], fold_number)  # corresponding outputs
          }
          data.frame(temp_input, output = temp_output) %>% return
        }
      })
      names(list_all_predicted_lkd) <- full_node_names

      # make sure we force update this after each updating step; this helps speed up updater$check_convergence
      if (!is.null(list_D) & !is.null(result) & update == F) {
        return(result)
      } else {
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
        # substitution estimator
        vec_est <- left_join(obs_data[, tmle_task$npsem[[1]]$variables], library_L0)$output
        psi <- mean(vec_est)

        list_D <- self$clever_covariates(tmle_task, fold_number, submodel_type = "EIC")
        list_D[[1]] <- vec_est - psi
        list_D <- list_D[-which(names(list_D) == "IC")]

        vec_D <- list_D %>% compact %>% pmap_dbl(sum)
        IC <- vec_D

        result <- list(psi = psi, IC = IC
                       # , full_IC = list_D
        )

        # these are cached; unless likelihood is updated, or we force it to update, they shouldn't be changed
        if (fold_number == "full") {
          private$.list_D <- list_D
          private$.result <- result
        } else if (fold_number == "validation") {
          private$.list_D_val <- list_D
          private$.result_val <- result
        }

        return(result)
      }
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
    intervention_list_treatment = function() {
      return(self$cf_likelihood_treatment$intervention_list)
    },
    intervention_list_control = function() {
      return(self$cf_likelihood_control$intervention_list)
    },
    list_D = function() {
      return(private$.list_D)
    },
    list_D_val = function() {
      return(private$.list_D_val)
    },
    update_nodes = function() {
      # if (is.null(tmle_task)) {
      tmle_task <- self$observed_likelihood$training_task
      # }
      temp_node_names <- names(tmle_task$npsem)
      loc_A <- grep("A", temp_node_names)
      if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
      nodes_to_update <- temp_node_names[if_not_0 & !((1:length(temp_node_names)) %in% loc_A)]
      # nodes_to_update <- nodes_to_update[-length(nodes_to_update)]
      return(nodes_to_update)
    },
    list_newH = function() {
      return(private$.list_newH)
    },
    list_newH_val = function() {
      return(private$.list_newH_val)
    },
    list_newH_EIC = function() {
      return(private$.list_newH_EIC)
    },
    list_newH_EIC_val = function() {
      return(private$.list_newH_EIC_val)
    }
  ),
  private = list(
    .type = "mediation",
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL,
    .list_newH = NULL,  # the clever covariates for logistic TMLE
    .list_newH_val = NULL,
    .list_newH_EIC = NULL,  # the clever covariates as the EIC
    .list_newH_EIC_val = NULL,
    .list_D = NULL,
    .list_D_val = NULL,
    .result = NULL,
    .result_val = NULL,
    .submodel_type_supported = c("EIC")
  )
)




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
Param_mediation_projection <- R6Class(
  classname = "Param_mediation_projection",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list_treatment, intervention_list_control, outcome_node = "Y", static_likelihood = NULL, n_resampling = NULL) {
      if(inherits(observed_likelihood, "Targeted_Likelihood")){
        fold_number <- observed_likelihood$updater$update_fold
      } else {
        fold_number <- "full"
      }

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

      private$.static_likelihood <- static_likelihood
      private$.update_nodes <- c(Z_nodes, RLY_nodes)

      super$initialize(observed_likelihood, list(), outcome_node = outcome_node)
      private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
      private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)

      tmle_task <- observed_likelihood$training_task
      obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t))
      obs_variable_names <- colnames(obs_data)

      if (!is.null(n_resampling)) {  # use expanded Monte Carlo samples to train the HAL projection
        # temp_input <- generate_Zheng_data(B = n_resampling, tau = 1, if_LY_misspec = T) %>% data.frame()
        # names(temp_input)[grep("L1_1", names(temp_input))] <- "L_1"

        temp_input <- tmle_task$get_tmle_node(temp_node_names[1])
        temp_input <- rbind(
          temp_input[c(
            sample(nrow(temp_input), abs(round(n_resampling)), replace = T)
            # 1:n_resampling
          ), ]
        )
        for (i in 2:length(temp_node_names)) {
          temp_input <- cbind(temp_input, 1) %>% as.data.frame()
          names(temp_input)[ncol(temp_input)] <- temp_node_names[i]
          temp_task <- tmle3_Task$new(temp_input, tmle_task$npsem[1:i])
          temp_input[,ncol(temp_input)] <- rbinom(temp_task$nrow, 1, static_likelihood$get_likelihood(temp_task, node = temp_node_names[i]))
        }


        # temp_input <- expand_values(variables = obs_variable_names)  # all possible inputs
        # temp_input <- obs_data
        # temp_input <- rbind(
        #   # obs_data,
        #                     temp_input[c(
        #                            # sample(nrow(temp_input), abs(round(n_resampling)), replace = T)
        #                       1:n_resampling
        #                            ), ]
        # )
        temp_input <- data.frame(temp_input,
                                 id = 1:nrow(temp_input),
                                 t = 0)
        # temp_input <- rbind(temp_input, obs_data)
        temp_task <- tmle3_Task$new(temp_input, tmle_task$npsem)

        # private$.cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(temp_task)[[1]]
        # private$.cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(temp_task)[[1]]
        # Train the gradient
        private$.gradient <- Gradient$new(observed_likelihood,
                                          ipw_args = list(cf_likelihood_treatment = self$cf_likelihood_treatment,
                                                          cf_likelihood_control = self$cf_likelihood_control,
                                                          intervention_list_treatment = self$intervention_list_treatment,
                                                          intervention_list_control = self$intervention_list_control,
                                                          # cf_task_treatment = self$cf_task_treatment,
                                                          # cf_task_control = self$cf_task_control,
                                                          static_likelihood = self$static_likelihood
                                          ),
                                          projection_task_generator = gradient_generator_middle,
                                          target_nodes =  self$update_nodes)

        private$.gradient$train_projections(temp_task, fold_number = fold_number)
      } else {
        # todo: extend for stochastic
        # private$.cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(observed_likelihood$training_task)[[1]]
        # private$.cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(observed_likelihood$training_task)[[1]]
        # Train the gradient
        private$.gradient <- Gradient$new(observed_likelihood,
                                          ipw_args = list(cf_likelihood_treatment = self$cf_likelihood_treatment,
                                                          cf_likelihood_control = self$cf_likelihood_control,
                                                          intervention_list_treatment = self$intervention_list_treatment,
                                                          intervention_list_control = self$intervention_list_control,
                                                          # cf_task_treatment = self$cf_task_treatment,
                                                          # cf_task_control = self$cf_task_control,
                                                          static_likelihood = self$static_likelihood
                                          ),
                                          projection_task_generator = gradient_generator_middle,
                                          target_nodes =  self$update_nodes)

        private$.gradient$train_projections(self$observed_likelihood$training_task, fold_number = fold_number)
      }









      setattr(self$observed_likelihood, "target_nodes", self$update_nodes)
      self$observed_likelihood$get_likelihoods(self$observed_likelihood$training_task, fold_number = fold_number)
      for (node in self$update_nodes) {
        temp_long_task <- private$.gradient$expand_task(observed_likelihood$training_task, node)
        self$observed_likelihood$get_likelihood(temp_long_task, node, fold_number)
        self$observed_likelihood$get_likelihood(observed_likelihood$training_task, node, fold_number)
        # private$.gradient$expand_task(private$.cf_task_treatment, node)
        # private$.gradient$expand_task(private$.cf_task_control, node)
      }

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
            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number)  # corresponding outputs
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
            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number)  # corresponding outputs
          } else {
            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number)  # corresponding outputs
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

      EIC <- cbind(EIC, vec_est - psi)

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
Param_mediation_projection_survival <- R6Class(
  classname = "Param_mediation_projection",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list_treatment, intervention_list_control, outcome_node = "Y", static_likelihood = NULL, n_resampling = NULL) {
      if(inherits(observed_likelihood, "Targeted_Likelihood")){
        fold_number <- observed_likelihood$updater$update_fold
      } else {
        fold_number <- "full"
      }

      temp_node_names <- names(observed_likelihood$training_task$npsem)
      temp_node_names <- temp_node_names[-grep("delta_", temp_node_names)]  # remove delta nodes for wide format fitting
      loc_A_E <- grep("A_E", temp_node_names)
      loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
      loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
      if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
      tau <- as.numeric(last(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2])))

      A_E_nodes <- grep("A_E", temp_node_names, value = T)
      Z_nodes <- grep("Z", temp_node_names, value = T)
      RLY_nodes <- grep("(R|L|Y).[1-9]$", temp_node_names, value = T)

      private$.static_likelihood <- static_likelihood
      private$.update_nodes <- c(Z_nodes, RLY_nodes)

      private$.supports_outcome_censoring <- T
      super$initialize(observed_likelihood, list(), outcome_node = outcome_node)

      private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
      private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)

      tmle_task <- observed_likelihood$training_task
      obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t)) %>% dplyr::select(!starts_with("delta"))
      obs_variable_names <- colnames(obs_data)
      loc_A_C_or_Y <- grep("A_C_|Y_", temp_node_names)

      if (!is.null(n_resampling)) {  # use expanded Monte Carlo samples to train the HAL projection
        # temp_input <- generate_Zheng_data(B = n_resampling, tau = 1, if_LY_misspec = T) %>% data.frame()
        # names(temp_input)[grep("L1_1", names(temp_input))] <- "L_1"
        temp_input <- tmle_task$get_tmle_node(temp_node_names[1])
        temp_input <- rbind(
          temp_input[c(
            sample(nrow(temp_input), abs(round(n_resampling)), replace = T)
          ), ]
        )
        for (i in 2:length(temp_node_names)) {
          temp_input <- cbind(temp_input, 1) %>% as.data.frame()
          names(temp_input)[ncol(temp_input)] <- temp_node_names[i]
          delta_vars <- names(sapply(paste0("delta_", temp_node_names[1:i]), function(x) grep(x, names(tmle_task$npsem))) %>% compact %>% unlist)
          if (length(delta_vars) > 0) {
            full_input <- cbind(temp_input, matrix(T, 1, length(delta_vars)))
            colnames(full_input)[(ncol(temp_input) + 1):ncol(full_input)] <- delta_vars
          } else {
            full_input <- temp_input
          }
          temp_task <- tmle3_Task_drop_censored$new(full_input, tmle_task$npsem[c(1:i,
                                                                                  sapply(paste0("delta_", temp_node_names[1:i]), function(x) grep(x, names(tmle_task$npsem))) %>% compact %>% unlist
          )])
          temp_input[,ncol(temp_input)] <- rbinom(temp_task$nrow, 1, static_likelihood$get_likelihood(temp_task, node = temp_node_names[i]))
        }
        loc_A_C_or_Y_needed <- loc_A_C_or_Y[loc_A_C_or_Y < i]
        loc_current_var <- which(colnames(temp_input) == tmle_task$npsem[[i]]$variables)
        if (length(loc_A_C_or_Y_needed) > 0) for (s in loc_A_C_or_Y_needed) {
          loc_var <- which(colnames(temp_input) == tmle_task$npsem[[s]]$variables)
          if_censored <- which(temp_input[, loc_var] == 0)
          if_censored <- if_censored[!is.na(if_censored)]
          temp_input[if_censored, (loc_var+1):loc_current_var] <- NA
        }
        for (name in colnames(temp_input)) {
          if(any(is.na(temp_input[, name]))) {
            temp_input <- cbind(temp_input, !is.na(temp_input[, name]))
            colnames(temp_input)[ncol(temp_input)] <- paste0("delta_", name)
          }
        }
        temp_input <- data.frame(temp_input,
                                 id = 1:nrow(temp_input),
                                 t = 0)
        # temp_input <- rbind(temp_input, obs_data)
        temp_task <- tmle3_Task_drop_censored$new(temp_input, tmle_task$npsem)

        # private$.cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(temp_task)[[1]]
        # private$.cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(temp_task)[[1]]
        # Train the gradient
        private$.gradient <- Gradient_wide$new(observed_likelihood,
                                               ipw_args = list(cf_likelihood_treatment = self$cf_likelihood_treatment,
                                                               cf_likelihood_control = self$cf_likelihood_control,
                                                               intervention_list_treatment = self$intervention_list_treatment,
                                                               intervention_list_control = self$intervention_list_control,
                                                               # cf_task_treatment = self$cf_task_treatment,
                                                               # cf_task_control = self$cf_task_control,
                                                               static_likelihood = self$static_likelihood
                                               ),
                                               projection_task_generator = gradient_generator_middle_survival,
                                               target_nodes =  self$update_nodes)

        private$.gradient$train_projections(temp_task, fold_number = fold_number)
      } else {
        # todo: extend for stochastic
        # private$.cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(observed_likelihood$training_task)[[1]]
        # private$.cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(observed_likelihood$training_task)[[1]]
        # Train the gradient
        private$.gradient <- Gradient_wide$new(observed_likelihood,
                                               ipw_args = list(cf_likelihood_treatment = self$cf_likelihood_treatment,
                                                               cf_likelihood_control = self$cf_likelihood_control,
                                                               intervention_list_treatment = self$intervention_list_treatment,
                                                               intervention_list_control = self$intervention_list_control,
                                                               # cf_task_treatment = self$cf_task_treatment,
                                                               # cf_task_control = self$cf_task_control,
                                                               static_likelihood = self$static_likelihood
                                               ),
                                               projection_task_generator = gradient_generator_middle_survival,
                                               target_nodes =  self$update_nodes)

        private$.gradient$train_projections(self$observed_likelihood$training_task, fold_number = fold_number)
      }

      setattr(self$observed_likelihood, "target_nodes", self$update_nodes)
      for (node in temp_node_names) {
        self$observed_likelihood$get_likelihood(self$observed_likelihood$training_task, node = node, fold_number = fold_number)
      }
      for (node in self$update_nodes) {
        temp_long_task <- private$.gradient$expand_task(observed_likelihood$training_task, node)
        self$observed_likelihood$get_likelihood(temp_long_task, node, fold_number)
        self$observed_likelihood$get_likelihood(observed_likelihood$training_task, node, fold_number)
        # private$.gradient$expand_task(private$.cf_task_treatment, node)
        # private$.gradient$expand_task(private$.cf_task_control, node)
      }

      list_all_predicted_lkd <- lapply(1:length(temp_node_names), function(loc_node) {
        if (loc_node > 1) {
          # currently only support univariate node for t>0
          current_variable <- tmle_task$npsem[[loc_node]]$variables
          loc_impute <- grep("Y_|A_C_", temp_node_names)  # remain alive and uncensored before current variable
          loc_impute <- loc_impute[loc_impute < loc_node]
          if (length(loc_impute) == 0) {
            temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)])  # all possible inputs
          } else {
            temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)],
                                        rule_variables = sapply(loc_impute, function(s) tmle_task$npsem[[s]]$variables),
                                        rule_values = rep(1, length(loc_impute))
            )  # all possible inputs
          }
          delta_vars <- names(sapply(paste0("delta_", temp_node_names[1:loc_node]), function(x) grep(x, names(tmle_task$npsem))) %>% compact %>% unlist)
          if (length(delta_vars) > 0) {
            temp_input <- cbind(temp_input, matrix(T, 1, length(delta_vars)))
            colnames(temp_input)[(ncol(temp_input) - length(delta_vars) + 1):ncol(temp_input)] <- delta_vars
          }

          temp_task <- tmle3_Task$new(temp_input, tmle_task$npsem[c(1:loc_node,
                                                                    sapply(paste0("delta_", temp_node_names[1:loc_node]), function(x) grep(x, names(tmle_task$npsem))) %>% compact %>% unlist
          )])
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
            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number)  # corresponding outputs
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

      EICs_impute <- lapply(update_nodes, function(node){
        temp_vec <- self$gradient$compute_component(tmle_task, node, fold_number = fold_number)$EIC
        if (!is.null(tmle_task$npsem[[node]]$censoring_node)) {
          full_vec <- if_observed <- tmle_task$get_tmle_node(tmle_task$npsem[[node]]$censoring_node$variables)
          full_vec[if_observed] <- temp_vec
          full_vec[!if_observed] <- 0
          temp_vec <- full_vec
        }
        return(temp_vec)
      })
      EICs[[length(EICs) + 1]] <- do.call(cbind, EICs_impute)
      EICs[[length(EICs)]]
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

      temp_node_names <- names(self$observed_likelihood$training_task$npsem)
      temp_node_names <- temp_node_names[-grep("delta_", temp_node_names)]  # remove delta nodes for wide format fitting
      loc_A_E <- grep("A_E", temp_node_names)
      loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
      loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
      if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
      tau <- last(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2]))

      obs_data <- tmle_task$data %>% as.data.frame %>% dplyr::select(-c(id, t)) %>% dplyr::select(!starts_with("delta"))
      obs_variable_names <- colnames(obs_data)
      list_all_predicted_lkd <- lapply(1:length(temp_node_names), function(loc_node) {
        if (loc_node > 1) {
          # currently only support univariate node for t>0
          current_variable <- tmle_task$npsem[[loc_node]]$variables
          loc_impute <- grep("Y_|A_C_", temp_node_names)  # remain alive and uncensored before current variable
          loc_impute <- loc_impute[loc_impute < loc_node]
          if (length(loc_impute) == 0) {
            temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)])  # all possible inputs
          } else {
            temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)],
                                        rule_variables = sapply(loc_impute, function(s) tmle_task$npsem[[s]]$variables),
                                        rule_values = rep(1, length(loc_impute))
            )  # all possible inputs
          }
          delta_vars <- names(sapply(paste0("delta_", temp_node_names[1:loc_node]), function(x) grep(x, names(tmle_task$npsem))) %>% compact %>% unlist)
          if (length(delta_vars) > 0) {
            temp_input <- cbind(temp_input, matrix(T, 1, length(delta_vars)))
            colnames(temp_input)[(ncol(temp_input) - length(delta_vars) + 1):ncol(temp_input)] <- delta_vars
          }

          temp_task <- tmle3_Task$new(temp_input, tmle_task$npsem[c(1:loc_node,
                                                                    sapply(paste0("delta_", temp_node_names[1:loc_node]), function(x) grep(x, names(tmle_task$npsem))) %>% compact %>% unlist
          )])
          temp_target_node <- intersect(self$update_nodes, temp_node_names[loc_node])
          if (length(temp_target_node) == 1) {
            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number)  # corresponding outputs
          } else {
            # A nodes won't get updated
            temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number)  # corresponding outputs
          }
          data.frame(temp_input[1:which(obs_variable_names == current_variable)], output = temp_output) %>% return
        }
      })

      intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
      intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
      intervention_levels_treat <- map_dbl(self$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
      intervention_levels_control <- map_dbl(self$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
      # nodes to integrate out in the target identification
      # only support univaraite node for now; assume treatment level is one
      loc_impute <- grep("Y_|A_C_", temp_node_names)  # in integral, Y and A_C always 1
      all_possible_RZLY_1 <- expand_values(variables = obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ),
                                           rule_variables = c(intervention_variables, sapply(loc_impute, function(s) tmle_task$npsem[[s]]$variables)),
                                           rule_values = c(intervention_levels_treat, rep(1, length(loc_impute)))
      )
      all_possible_RZLY_0 <- expand_values(variables = obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ),
                                           rule_variables = c(intervention_variables, sapply(loc_impute, function(s) tmle_task$npsem[[s]]$variables)),
                                           rule_values = c(intervention_levels_control, rep(1, length(loc_impute)))
      )
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

      EIC <- cbind(EIC, vec_est - psi)

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
    .type = "mediation_survival",
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
