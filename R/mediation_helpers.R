#' Helper Functions for Longitudinal Mediation
#'
#' Handles the ARZL(Y) structure
#'
#' @param data a \code{data.frame}, or \code{data.table} containing data for use in estimation
#' @param node_list a list of character vectors, listing the variables that comprise each node
#' @param variable_types a list of variable types, one for each node. If missing, variable types will be guessed
#' @param tmle_task a \code{\link{tmle3_Task}} as constructed via \code{point_tx_task}
#' @param learner_list a list of sl3 learners, one for A and one for Y to be used for likelihood estimation
#' @param ... extra arguments.
#' @export
#' @rdname point_tx
mediation_npsem <- function(node_list, variable_types = NULL) {
  # make tmle_task; NULL variable_types will be learned by sl3 function
  npsem <- c(define_node("L_0", node_list$L_0, variable_type = variable_types$L_0),
             lapply(2:length(node_list), function(k) {
                 define_node(names(node_list)[k],
                             node_list[[k]],
                             names(node_list)[1:(k-1)],
                             variable_type = variable_types[[ names(node_list)[k] ]])
             })
  )

  return(npsem)
}



#' @export
#' @rdname point_tx
mediation_task <- function(data, node_list, variable_types = NULL, ...) {
  setDT(data)

  npsem <- mediation_npsem(node_list, variable_types)

  if (!is.null(node_list$id)) {
    tmle_task <- tmle3_Task$new(data, npsem = npsem, id = node_list$id, ...)
  } else {
    tmle_task <- tmle3_Task$new(data, npsem = npsem, ...)
  }

  return(tmle_task)
}

#' @export
#' @rdname point_tx
mediation_task_drop_censored <- function(data, node_list, variable_types = NULL, ...) {
  setDT(data)

  npsem <- mediation_npsem(node_list, variable_types)

  if (!is.null(node_list$id)) {
    tmle_task <- tmle3_Task_drop_censored$new(data, npsem = npsem, id = node_list$id, ...)
  } else {
    tmle_task <- tmle3_Task_drop_censored$new(data, npsem = npsem, ...)
  }

  return(tmle_task)
}

#' @export
#' @rdname point_tx
mediation_likelihood <- function(tmle_task, learner_list) {
  factor_list <- list()

  # covariates
  L_0_factor <- define_lf(LF_emp, "L_0")
  factor_list[[1]] <- L_0_factor

  # treatment (bound likelihood away from 0 (and 1 if binary))
  # A_type <- tmle_task$npsem[["A"]]$variable_type
  A_type <- tmle_task$npsem[[ grep("A", names(tmle_task$npsem))[1] ]]$variable_type  # evaluate the first A node for now
  if (A_type$type == "continous") {
    A_bound <- c(1 / tmle_task$nrow, Inf)
  } else if (A_type$type %in% c("binomial", "categorical")) {
    A_bound <- 0.025
  } else {
    A_bound <- NULL
  }

  temp_names <- names(tmle_task$npsem)
  loc_A <- grep("A", temp_names)
  factor_list[loc_A] <- lapply(loc_A, function(k) {
    LF_fit$new(temp_names[k], learner = learner_list[[ temp_names[k] ]], bound = A_bound, type = "density")
  })

  # others
  loc_others <- (1:length(temp_names))[-c(grep("A", temp_names), 1)]
  factor_list[loc_others] <- lapply(loc_others, function(k) {
    LF_fit$new(temp_names[k], learner = learner_list[[ temp_names[k] ]], type = "density")
  })

  likelihood_def <- Likelihood$new(factor_list)
  likelihood <- likelihood_def$train(tmle_task)
  return(likelihood)
}


#' @export
get_current_Q <- function(loc_node, which_Q,
                          tmle_task, obs_variable_names,
                          intervention_variables, intervention_levels_treat, intervention_levels_control,
                          list_all_predicted_lkd
) {
  temp_node_names <- names(tmle_task$npsem)  # this has to be the obs task, but only the node names are needed
  loc_A <- grep("A", temp_node_names)
  loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) !(strsplit(s, "_")[[1]][1] %in% c("A", "Z")) & strsplit(s, "_")[[1]][2] != 0))
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))

  # ind_var is the order among variables
  ind_var <- loc_current_var <- which(obs_variable_names == tmle_task$npsem[[loc_node]]$variables)
  # decide which Q to get
  ind_var <- ind_var + which_Q

  intervention_variables_loc_needed <- intervention_variables_loc[intervention_variables_loc < ind_var]

  # only update t!=0 and non-A nodes
  if (loc_node %in% c(loc_Z, loc_RLY)) {
    df_to_update <- temp_current <- list_all_predicted_lkd[[loc_node]] # decide input matrix; the last column might not be fully used

    # get Q current
    if (ind_var <= length(obs_variable_names)) {  # Q_X, except for Q_R_tao+1
      data_temp <- df_to_update[1:(ind_var-1)]
      # generate all possible 0, 1 valued rlz; set A to 0 first; set y_tau to always 1;
      all_possible_rlz_1 <- expand_values(obs_variable_names, to_drop = c(1:(ind_var-1) ),
                                          A = 1,
                                          rule_variables = c(last(obs_variable_names)), rule_values = c(1))
      all_possible_rlz_0 <- expand_values(obs_variable_names, to_drop = c(1:(ind_var-1) ),
                                          A = 0,
                                          rule_variables = c(last(obs_variable_names)), rule_values = c(1))
      # for each observed L_0 vector, generate all needed combinations, one version for A = 1, one version for A = 0
      unique_input <- data_temp[1:(ind_var-1)] %>% unique
      library_output <- data.frame(unique_input, output =
                                     map_dbl(1:nrow(unique_input), function(which_row) {
                                       # probs in the integrals, A=1 or A=0 is inserted
                                       temp_all_comb_0 <- temp_all_comb_1 <- unique_input[which_row, ]
                                       if (length(all_possible_rlz_0) > 0) {
                                         temp_all_comb_0 <- suppressWarnings(cbind(temp_all_comb_0, all_possible_rlz_0))
                                         temp_all_comb_1 <- suppressWarnings(cbind(temp_all_comb_1, all_possible_rlz_1))
                                       }
                                       if (length(intervention_variables_loc_needed) > 0) {
                                         for (i in 1:length(intervention_variables_loc_needed)) {
                                           temp_all_comb_0[, intervention_variables_loc[i]] <- intervention_levels_control[i]
                                           temp_all_comb_1[, intervention_variables_loc[i]] <- intervention_levels_treat[i]
                                         }
                                       }
                                       # for all non-A, non-0 variables, calculate the variable by rule
                                       # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
                                       loc_Z_need <- loc_Z[loc_Z >= loc_node + which_Q]  # only integrate out future variables
                                       temp_list_0 <- lapply(loc_Z_need,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_0, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       loc_other_needed <- loc_RLY[loc_RLY >= loc_node + which_Q]
                                       temp_list_1 <- lapply(loc_other_needed,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_1, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       temp_list <- c(temp_list_0, temp_list_1)
                                       pmap_dbl(temp_list, prod) %>% sum %>% return
                                     })
      )
      Q_current <- left_join(data_temp, library_output)$output
    } else {
      Q_current <- df_to_update[, ind_var-1]  # Q_R_tao+1 is just the output
    }

    return(Q_current)
  } else {
    return(NULL)
  }
}


#' @export
get_current_H <- function(loc_node,
                          tmle_task, obs_variable_names,
                          intervention_variables, intervention_levels_treat, intervention_levels_control,
                          list_all_predicted_lkd
) {
  temp_node_names <- names(tmle_task$npsem)  # this has to be the obs task, but only the node names are needed
  loc_A <- grep("A", temp_node_names)
  loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) !(strsplit(s, "_")[[1]][1] %in% c("A", "Z")) & strsplit(s, "_")[[1]][2] != 0))
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))

  # ind_var is the order among variables
  ind_var <- loc_current_var <- which(obs_variable_names == tmle_task$npsem[[loc_node]]$variables)
  intervention_variables_loc_needed <- intervention_variables_loc[intervention_variables_loc < ind_var]

  temp_current <- list_all_predicted_lkd[[loc_node]] # decide input matrix; the last column might not be fully used

  # get H current
  {
    data_temp <- temp_current[1:ind_var]  # take =1 probs
    temp_ind <- ind_var  # the variable loc

    all_observed_1 <- all_observed_0 <- data_temp
    if (length(intervention_variables_loc_needed) > 0) {
      for (i in 1:length(intervention_variables_loc_needed)) {
        all_observed_0[, intervention_variables_loc[i]] <- intervention_levels_control[i]
        all_observed_1[, intervention_variables_loc[i]] <- intervention_levels_treat[i]
      }
    }

    # R and L (including Y before t = tau) (or not Z not A not t=0, not Y_tau) can be calculated in the same way
    if (loc_node %in% loc_RLY) {
      loc_A_needed <- loc_A[loc_A < loc_node]  # all needed A nodes
      loc_Z_needed <- loc_Z[loc_Z < loc_node]  # all needed Z nodes

      # all predicted probs at these locs are needed
      # needed inputs are in either all_observed_1 or all_observed_0

      # these A probs will be taken as product
      part_A <- lapply(loc_A_needed, function(k) {
        left_join(all_observed_1, list_all_predicted_lkd[[k]])$output
      }) %>% pmap_dbl(prod)

      # these ratios Z probs will be taken as product
      part_Z <- ifelse_vec(length(loc_Z_needed) == 0, rep(1, length(part_A)),
                           lapply(loc_Z_needed, function(k) {
                             left_join(all_observed_0, list_all_predicted_lkd[[k]])$output /
                               left_join(all_observed_1, list_all_predicted_lkd[[k]])$output
                           }) %>% pmap_dbl(prod))

      H_current <- ifelse(
        apply(data_temp[intervention_variables_loc_needed] == intervention_levels_treat, 1, prod) == 1  # this is the A_ind decided by the library tasks
        , 1/part_A*part_Z, 0) %>% as.vector
    }
    # Z nodes
    if (loc_node %in% loc_Z) {
      loc_A_needed <- loc_A[loc_A < loc_node]  # all needed A nodes
      loc_RLY_needed <- loc_RLY[loc_RLY < loc_node]

      # these A probs will be taken as product
      part_A <- lapply(loc_A_needed, function(k) {
        left_join(all_observed_0, list_all_predicted_lkd[[k]])$output
      }) %>% pmap_dbl(prod)

      # these ratios Z probs will be taken as product
      part_LR <- lapply(loc_RLY_needed, function(k) {
        left_join(all_observed_1, list_all_predicted_lkd[[k]])$output /
          left_join(all_observed_0, list_all_predicted_lkd[[k]])$output
      }) %>% pmap_dbl(prod)

      H_current <- ifelse(
        # ZW todo: stochastic intervention
        apply(data_temp[intervention_variables_loc_needed] == intervention_levels_control, 1, prod) == 1
        , 1/part_A*part_LR, 0)
    }
  }
  return(H_current)
}


















#' @export
get_obs_Q_list <- function(tmle_task, obs_data,
                      intervention_variables, intervention_levels_treat, intervention_levels_control,
                      list_all_predicted_lkd) {
  # decide validation or full by list_all_predicted_lkd
  obs_variable_names <- colnames(obs_data)
  temp_node_names <- names(tmle_task$npsem)
  loc_A <- grep("A", temp_node_names)
  loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) !(strsplit(s, "_")[[1]][1] %in% c("A", "Z")) & strsplit(s, "_")[[1]][2] != 0))
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
  list_Q <- list()
  for (loc_node in 1:length(temp_node_names)) {
    if (loc_node %in% c(loc_Z, loc_RLY)) {
      # generate all possible 0, 1 valued rlz; set A to 0 first; set y_tau to always 1;
      loc_current_var <- which(obs_variable_names == tmle_task$npsem[[loc_node]]$variables)
      all_possible_RZLY_1 <- expand_values(obs_variable_names, to_drop = c(1:(loc_current_var-1) ),
                                           rule_variables = c(last(obs_variable_names), intervention_variables),
                                           rule_values = c(1, intervention_levels_treat))
      all_possible_RZLY_0 <- expand_values(obs_variable_names, to_drop = c(1:(loc_current_var-1) ),
                                           rule_variables = c(last(obs_variable_names), intervention_variables),
                                           rule_values = c(1, intervention_levels_control))

      # for each observed L_0 vector, generate all needed combinations, one version for A = 1, one version for A = 0
      unique_input <- obs_data[, 1:(loc_current_var-1)] %>% unique
      library_output <- data.frame(unique_input, output =
                                     map_dbl(1:nrow(unique_input), function(which_row) {
                                       # probs in the integrals, A=1 or A=0 is inserted
                                       temp_all_comb_0 <- data.frame(unique_input[which_row, ], all_possible_RZLY_0)
                                       temp_all_comb_1 <- data.frame(unique_input[which_row, ], all_possible_RZLY_1)
                                       for (i in 1:length(intervention_variables_loc)) {
                                         temp_all_comb_0[, intervention_variables_loc[i]] <- intervention_levels_control[i]
                                         temp_all_comb_1[, intervention_variables_loc[i]] <- intervention_levels_treat[i]
                                       }
                                       # for all non-A, non-0 variables, calculate the variable by rule
                                       # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
                                       loc_Z_needed <- loc_Z[loc_Z >= loc_node]  # only product children variables
                                       # temp_list_0 <- lapply(loc_Z_needed,
                                       temp_list_0 <- lapply(loc_Z_needed,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_0, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       loc_RLY_needed <- loc_RLY[loc_RLY >= loc_node]
                                       temp_list_1 <- lapply(loc_RLY_needed,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_1, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       temp_list <- c(temp_list_0, temp_list_1)
                                       pmap_dbl(temp_list, prod) %>% sum %>% return
                                     })
      )
      list_Q[[loc_node]] <- left_join(obs_data[, 1:(loc_current_var-1)], library_output)$output
    }
  }
  return(list_Q)
}





#' @export
get_obs_H_list <- function(tmle_task, obs_data, current_likelihood,
                      cf_task_treatment, cf_task_control,
                      intervention_variables, intervention_levels_treat, intervention_levels_control,
                      fold_number = "full"
) {
  obs_variable_names <- colnames(obs_data)
  temp_node_names <- names(tmle_task$npsem)
  loc_A <- grep("A", temp_node_names)
  loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))

  list_H <- list()  # get a list of corresponding H covariates; ordered by nodes, not variables
  for (temp_ind in loc_RLY) {  # calculate RLY nodes
    loc_A_needed <- loc_A[loc_A < temp_ind]  # all needed A nodes
    loc_Z_needed <- loc_Z[loc_Z < temp_ind]  # all needed Z nodes
    A_ind <- apply(sapply(loc_A_needed, function(k) {
      obs_data[[tmle_task$npsem[[k]]$variables]] == intervention_levels_treat[tmle_task$npsem[[k]]$name]
    }), 1, all)
    part_A <- lapply(loc_A_needed, function(k) current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k], fold_number)) %>% pmap_dbl(prod)  # this is the likelihood of being 1
    part_Z <- lapply(loc_Z_needed, function(k) {
      current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k], fold_number) /
        current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k], fold_number)
    }) %>% pmap_dbl(prod)
    if(length(part_Z) == 0) part_Z <- 1
    list_H[[temp_ind]] <- ifelse(A_ind, 1/part_A*part_Z, 0) %>% as.vector
  }
  for (temp_ind in loc_Z) {    # calculate Z nodes
    loc_A_needed <- loc_A[loc_A < temp_ind]  # all needed A nodes
    loc_RLY_needed <- loc_RLY[loc_RLY < temp_ind]
    A_ind <- apply(sapply(loc_A_needed, function(k) {
      obs_data[[tmle_task$npsem[[k]]$variables]] == intervention_levels_control[tmle_task$npsem[[k]]$name]
    }), 1, all)
    part_A <- lapply(loc_A_needed, function(k) current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k], fold_number)) %>% pmap_dbl(prod)
    part_RLY <- lapply(loc_RLY_needed, function(k) {
      current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k], fold_number) /
        current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k], fold_number)
    }) %>% pmap_dbl(prod)
    if(length(part_RLY) == 0) part_RLY <- 1
    list_H[[temp_ind]] <- ifelse(A_ind, 1/part_A*part_RLY, 0) %>% as.vector
  }
  return(list_H)
}
