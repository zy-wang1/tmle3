self <- tmle_params[[1]]$gradient

node <- update_node
tmle_params[[1]]$gradient$likelihood$cache$get_update_step(self$likelihood$factor_list[[node]], long_task, fold_number, node = node)

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


tlik$cache$tasks %>% length
tlik$cache$tasks %>% lapply(function(x) x$data$A_1 %>% table)

self$likelihood$get_likelihood(long_task, node, fold_number = fold_number, drop_id = T, drop_time = T, drop = T  )
tlik$cache$tasks %>% length
tlik$cache$tasks %>% lapply(function(x) x$data$A_1 %>% table)

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


tlik$cache$tasks %>% length
tlik$cache$tasks %>% lapply(function(x) x$data$A_1 %>% table)


data <- IC_task$data

# TODO

variables <- node
#TODO check id order

data <- data.table(cbind(merge(long_task$data[,c("id", "trueid", "t")], cbind(long_task$get_tmle_node(node, format = T, include_id = T, include_time = T), long_preds), by = c("id", "t"))
))

idkey <- data$trueid
data$t <- NULL
data$id <- NULL

setnames(data, c("id", node, "pred"))

setkeyv(data, cols = c("id", node))

#TODO handle updating of expanded task
#This is done by stacking copies of the cdf

data <- dcast(data, as.formula(paste0("id ~ ", node)), value.var = "pred")
id <- data$id
data$id <- NULL
levels <- as.numeric(colnames(data))

cdf <- as.data.table(t(apply(data, 1, cumsum)))
setnames(cdf, as.character(levels))
#print(cdf)

if(long_task$uuid == tmle_task$uuid){
  #if expanded task is tmle_task then obtain then expand cdf to match
  #This ensures we dont have any recursion errors by expanding an expanded task

  match_index <- match(idkey, id)
  cdf <- cdf[match_index]
}


fit_obj <- private$.component_fits[[node]]
basis_list <- fit_obj$basis_list
coefs <- fit_obj$coefs
col_index <- which(colnames(IC_task$X) == node )


keep <- sapply(basis_list, function(b){
  col_index %in% b$cols
}) %>% unlist

basis_list <- basis_list[keep]
coefs <- coefs[c(T, keep)]

#Should already be sorted
X <- as.matrix(IC_task$X)
design <- as.data.table(as.matrix(hal9001::make_design_matrix(X, basis_list)))

diff_map <- sapply(seq_along(basis_list), function(i) {
  basis <- basis_list[[i]]
  result <- (list(which(levels == basis$cutoffs[which(basis$cols == col_index)])))

  return(result)
})


center_basis <- lapply(seq_along(diff_map), function(i){
  col_index <- diff_map[[i]]
  diff <- design[[as.integer(i)]] - 1 + cdf[[col_index]]
  set(design, , as.integer(i), diff)
})
min_val <- min(IC_task$X[[node]]) - 5
clean_basis <- function(basis){
  index = which(basis$cols == col_index)
  basis$cutoffs[index] <- min_val
  return(basis)
}
clean_list = lapply(basis_list, clean_basis)
#print(table(unlist( lapply(basis_list, `[[`, "cols"))))
clean_design <- hal9001::make_design_matrix(X, clean_list)
clean_design <- data.table(as.matrix(clean_design))

#print(as.data.table(clean_design))
#print(as.data.table(coefs))

mid_result <- as.matrix(design * clean_design)
result =  mid_result %*% coefs[-1]
out = list(col_index = col_index,Y = IC_task$Y, cdf = cdf,design = design,  mid_result = mid_result, coefs = coefs[-1], EIC = result)
return(out)
