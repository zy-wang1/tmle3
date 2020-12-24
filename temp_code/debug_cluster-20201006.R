{
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

  timepoint <- 1
  if_misspec <- T

set.seed(1234)
data_sim <- generate_Zheng_data(B = 1000, tau = timepoint, if_LY_misspec = if_misspec)
data_wide <- data.frame(data_sim)
node_list <- list(L_0 = c("L1_0", "L2_0"),
                  A_1 = "A_1",
                  R_1 = "R_1",
                  Z_1 = "Z_1",
                  L_1 = "L1_1",
                  Y_1 = "Y_1"
                  # ,
                  # A_2 = "A_2",
                  # R_2 = "R_2",
                  # Z_2 = "Z_2",
                  # L_2 = "L1_2",
                  # Y_2 = "Y_2"
)
node_list$L_1 <- "L_1"
names(data_sim[[2]])[grep("L1_1", names(data_sim[[2]]))] <- "L_1"
names(data_wide)[grep("L1_1", names(data_wide))] <- "L_1"
# node_list$L_2 <- "L_2"
# names(data_sim[[3]])[grep("L1_2", names(data_sim[[3]]))] <- "L_2"
# names(data_wide)[grep("L1_2", names(data_wide))] <- "L_2"
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
tmle_params_ini <- middle_spec$make_params(tmle_task, initial_likelihood, static_likelihood = NULL)
suppressMessages(
  nontargeting <- tmle_params_ini[[1]]$estimates(tmle_task)
)
# tmle_params_ini[[1]]$estimates()$IC %>% mean
# (tmle_params_ini[[1]]$estimates()$full_IC %>% lapply(mean) %>% compact() %>% unlist)[-c(1, 2)]


tlik <- Targeted_Likelihood$new(initial_likelihood,
                                submodel_type_by_node = "EIC" ,
                                updater = list(convergence_type = "scaled_var",
                                               constrain_step = T,
                                               optim_delta_epsilon = T,
                                               one_dimensional=T,
                                               delta_epsilon=function(x) {
                                                 ifelse(abs(mean(x %>% as.vector)) < sqrt(var(x %>% as.vector)/1000)/log(1000),
                                                        0.00000001,
                                                        ifelse(mean(x %>% as.vector) > 0, 0.01, -0.01)
                                                 )
                                                 # ifelse(mean(x %>% as.vector) > 0, 0.0001, -0.0001)
                                               },
                                               maxit=100
                                               ,
                                               cvtmle=F
                                ))
tmle_params <- middle_spec$make_params(tmle_task, tlik, if_projection = T, initial_likelihood)



self <- tmle_params[[1]]
observed_likelihood <- tlik



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

static_likelihood <- initial_likelihood

update_nodes <- c(Z_nodes, RLY_nodes)
intervention_list_treatment <- self$intervention_list_treatment
intervention_list_control <- self$intervention_list_control

}
cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)
cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(observed_likelihood$training_task)[[1]]
cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(observed_likelihood$training_task)[[1]]

gradient <- Gradient$new(observed_likelihood,
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
  temp_long_task <- gradient$expand_task(observed_likelihood$training_task, node)
  self$observed_likelihood$get_likelihood(temp_long_task, node, fold_number)
  self$observed_likelihood$get_likelihood(observed_likelihood$training_task, node, fold_number)
  # private$.gradient$expand_task(private$.cf_task_treatment, node)
  # private$.gradient$expand_task(private$.cf_task_control, node)
}


gradient$train_projections(self$observed_likelihood$training_task, fold_number = fold_number)


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
        temp_long_task <- gradient$expand_task(temp_task, node)
        self$observed_likelihood$get_likelihood(temp_long_task, node, fold_number)
      }
      temp_output <- self$observed_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number = "full")  # corresponding outputs
    } else {
      # A nodes won't get updated
      temp_output <- self$static_likelihood$get_likelihood(temp_task, node = temp_node_names[loc_node], fold_number = "full")  # corresponding outputs
    }
    data.frame(temp_input, output = temp_output) %>% return
  }
})

list_all_predicted_lkd %>% lapply(head)


fold_number <- "full"
node <- NULL
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



self <- gradient
node <- "Y_1"

print("Gradient")
print(node)
print(fold_number)
time <- tmle_task$npsem[[node]]$time

self$assert_trained()
#Converts squashed basis to R functions of tmle3_tasks

fit_obj <- self$component_fits[[node]]

long_task <- self$expand_task(tmle_task, node)

IC_task <- self$generate_task(tmle_task, node, include_outcome = F, fold_number = fold_number)


col_index <- which(colnames(IC_task$X) == node )

long_preds <- NULL


tryCatch({
  value_step1 <- self$likelihood$cache$get_update_step(self$likelihood$factor_list[[node]], tmle_task, fold_number, node = node)
  value_step2 <- self$likelihood$cache$get_update_step(self$likelihood$factor_list[[node]], long_task, fold_number, node = node)
  if(is.null(value_step1)){
    value_step1 <- 0
  }
  if(is.null(value_step2)){
    value_step2 <- 0
  }
  if(value_step1!=value_step2) {
    stop("Long_task and tmle_task out of sync.")
  }
  long_preds <- self$likelihood$get_likelihood(long_task, node, fold_number = fold_number, drop_id = T, drop_time = T, drop = T  )},
  error = function(e){
    #long_task is probably out of sync with tmle_task
    #Update it to the same level
    value_step <- self$likelihood$cache$get_update_step(self$likelihood$factor_list[[node]], tmle_task, fold_number, node = node)

    self$likelihood$sync_task(long_task, fold_number = fold_number, check = F, max_step = value_step)
    value_step1 <- self$likelihood$cache$get_update_step(self$likelihood$factor_list[[node]], tmle_task, fold_number, node = node)
    value_step2 <- self$likelihood$cache$get_update_step(self$likelihood$factor_list[[node]], long_task, fold_number, node = node)
    if(is.null(value_step1)){
      value_step1 <- 0
    }
    if(is.null(value_step2)){
      value_step2 <- 0
    }
    if(value_step1!=value_step2) {
      stop("Long_task and tmle_task out of sync.")
    }
    long_preds <<- self$likelihood$get_likelihood(long_task, node, fold_number = fold_number, drop_id = T, drop_time = T, drop = T, check_sync = F  )

  })


data <- IC_task$data

variables <- node
#TODO check id order

data <- data.table(cbind(merge(long_task$data[,c("id", "trueid", "t")], cbind(long_task$get_tmle_node(node, format = T, include_id = T, include_time = T), long_preds), by = c("id", "t"))
))


long_task$data[,c("id", "trueid", "t")] %>% head
long_task$data %>% head

long_task$data %>% filter(trueid == 100)
identical(long_task$data$id, long_task$get_tmle_node(node, format = T, include_id = T, include_time = T)$id)

cbind(long_task$get_tmle_node(node, format = T, include_id = T, include_time = T), long_preds)

cbind(long_task$get_tmle_node(node, format = T, include_id = T, include_time = T), long_preds) %>% filter(long_preds < 0.1)
sum(long_preds<0.1)
long_preds %>% unique


which(long_task$data$trueid == 999)
long_task$data %>% filter(trueid == 999)
self$likelihood$get_likelihood(long_task, node, fold_number = fold_number, drop_id = T, drop_time = T, drop = T)[which(long_task$data$trueid == 999)]



likelihood_factor <- self$likelihood$factor_list[[node]]
likelihood_values <- self$likelihood$cache$get_values(likelihood_factor, long_task, fold_number, node = paste0(node, collapse = "%"))
likelihood_values[which(long_task$data$trueid == 999)]




self <- gradient

key <- paste0(tmle_task$uuid, node, sep = "%")

cached_task <- get0(key, self$cache, inherits = FALSE)
variables <- tmle_task$npsem[[node]]$variables

# if(length(variables) >1) stop("Multivariate nodes not supported")
data <- tmle_task$data
data$trueid <- data$id
time <- tmle_task$npsem[[node]]$time
levels <- sort(unique(unlist(tmle_task$get_tmle_node(node))))  #data[, variables, with = F])))

long_data <- rbindlist(lapply(levels, function(level) {
  data <- copy(data)
  set(data , which(data$t == time), variables, level)
  data$levelcopy <- level
  return(data)
}))
long_data

long_data$id <-  paste(long_data$trueid,long_data$levelcopy, sep = "_")
long_data$id

# suppressWarnings(long_task <- tmle3_Task$new(long_data, tmle_task$npsem, id = "id", time = "t", force_at_risk = tmle_task$force_at_risk, summary_measure_columns = c(tmle_task$summary_measure_columns, "trueid")))

data <- long_data
npsem <- tmle_task$npsem
id = "id"
time = "t"
force_at_risk = tmle_task$force_at_risk
summary_measure_columns = c(tmle_task$summary_measure_columns, "trueid")
long_format = NULL
folds_for_ids = NULL

!inherits(data, "Shared_Data")

data[, id := as.factor(id)]
data$id %>% as.factor
data <- setkey(data, id, t)
shared_data <- data



long_data
long_task$data

long_data, tmle_task$npsem, id = "id", time = "t", force_at_risk = tmle_task$force_at_risk, summary_measure_columns = c(tmle_task$summary_measure_columns, "trueid")

setattr(long_task, "target_nodes", node)
if(is.null(attr(long_task, "target_nodes"))) {
  print(node)
  stop("wrong")
}
assign(key, long_task, self$cache)

private$.uuid_expanded_history[[long_task$uuid]] <- node

