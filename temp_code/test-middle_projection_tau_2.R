library(R6)  # R6class
library(data.table)  # setDT
library(sl3)
library(digest)
library(uuid)  # UUIDgenerate
library(delayed)  # bundle_delayed
library(assertthat)  # assert_that
library(speedglm)  # speedglm
# library(methods)  # is

library(dplyr)  # dplyr::select, to mask other packages

code_list <- list.files("./R", full.names = T)
for (code in code_list) source(code)
source("./temp_code/generate_data.R")

# record <- record_cor <- list()

set.seed(1234)

data_sim <- generate_Zheng_data(B = 1000, tau = 2, if_LY_misspec = F)
data_wide <- data.frame(data_sim)

node_list <- list(L_0 = c("L1_0", "L2_0"),
                  A_1 = "A_1",
                  R_1 = "R_1",
                  Z_1 = "Z_1",
                  L_1 = "L1_1",
                  Y_1 = "Y_1"
                  ,
                  A_2 = "A_2",
                  R_2 = "R_2",
                  Z_2 = "Z_2",
                  L_2 = "L1_2",
                  Y_2 = "Y_2"
)
node_list$L_1 <- "L_1"
names(data_sim[[2]])[grep("L1_1", names(data_sim[[2]]))] <- "L_1"
names(data_wide)[grep("L1_1", names(data_wide))] <- "L_1"
node_list$L_2 <- "L_2"
names(data_sim[[3]])[grep("L1_2", names(data_sim[[3]]))] <- "L_2"
names(data_wide)[grep("L1_2", names(data_wide))] <- "L_2"



middle_spec <- tmle_middle(
  treatment_level = 1,
  control_level = 0
)
tmle_task <- middle_spec$make_tmle_task(data_wide, node_list)
# choose base learners
lrnr_glm_fast <- Lrnr_glm_fast$new(outcome_type = "binomial")
learner_list <- lapply(1:length(tmle_task$npsem), function(s) Lrnr_sl$new(
  learners = list(
    lrnr_glm_fast
  )
))
names(learner_list) <- names(tmle_task$npsem)  # the first will be ignored; empirical dist. will be used for covariates
initial_likelihood <- middle_spec$make_initial_likelihood(
  tmle_task,
  learner_list
)

tlik <- Targeted_Likelihood$new(initial_likelihood,
                                submodel_type_by_node = "EIC" ,
                                updater = list(convergence_type = "scaled_var",
                                               constrain_step = T,
                                               optim_delta_epsilon = F,
                                               one_dimensional=TRUE,
                                               delta_epsilon = 0.01,
                                               maxit = 3))
tmle_params <- middle_spec$make_params(tmle_task, tlik, if_projection = T, initial_likelihood)
tmle_params[[1]]$estimates()$psi

aaa <- tlik$cache$tasks %>% names
tlik$cache$tasks %>% length
tlik$cache$tasks %>% lapply(function(x) x$data %>% dim)
tlik$updater$epsilons

tlik$updater$update(tlik, tmle_task)

tlik$updater$update_step(tlik, tmle_task)



self <- tlik

# todo: rethink which tasks need updates here
# tasks_at_step <- self$cache$tasks_at_step(step_number)
tasks_at_step <- self$cache$tasks
if (inherits(self$updater, "tmle3_Update")) {} else tasks_at_step <- tasks_at_step[map_dbl(tasks_at_step, ~nrow(.x$data)) == nrow(self$training_task$data)]


  # the default version
  # If task has attr target_nodes then only update these nodes.
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

  tasks_at_step <- tasks_at_step[to_update]
  task_updates <- tasks_at_step %>% lapply(function(x) if (update_node %in% names(x$npsem)) {
    self$updater$apply_update(x, self, fold_number, new_epsilon, update_node)
  })

  tasks_at_step[3] %>% lapply(function(x) if (update_node %in% names(x$npsem)) {
    self$updater$apply_update(x, self, fold_number, new_epsilon, update_node)
  })




  # the short task can't be updated
