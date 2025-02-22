
context("Fitting survival hazards via risk sets and pooled over time nodes")
library(simcausal)
D <- DAG.empty()
D <- D +
  node("W", distr = "runif", min = 0, max = 1.5) +
  node("A", distr = "rbinom", size = 1, prob = .15 + .5 * as.numeric(W > .75)) +
  node("Trexp", distr = "rweibull", shape = 1 + .5 * W, scale = 6) +
  node("Cweib", distr = "rweibull", shape = 1 + .5 * W, scale = 15) +
  node("T", distr = "rconst", const = min(ceiling(Trexp ), 12)) +
  node("C", distr = "rconst", const = ceiling(Cweib )) +
  # Observed random variable (follow-up time):
  node("T.tilde", distr = "rconst", const = ifelse(T <= C, T, C)) +
  # Observed random variable (censoring indicator, 1 - failure event, 0 - censored):
  node("Delta", distr = "rconst", const = ifelse(T <= C, 1, 0))
setD <- set.DAG(D)
dat <- sim(setD, n = 1000)
# only grab ID, W's, A, T.tilde, Delta
Wname <- grep("W", colnames(dat), value = TRUE)
dat <- dat[, c("ID", Wname, "A", "T.tilde", "Delta")]
dat
dat$id <- dat$ID
dat$ID <- NULL
dat$processN <- dat$Delta
dat$processA <- 1 - dat$Delta
dat$t <- dat$T.tilde
dat$T.tilde <- NULL

init_dat = copy(dat)
init_dat$t = 0
init_dat$processN <- 0
init_dat$processA <- 0
long_data <- data.table(rbind(init_dat, dat))
setkey(long_data, id, t)
long_data



# at_risk indicators for degeneracy
at_risk_A <- function(X,time, ...){
  index_strict <- X$t < time
  as.numeric(all(X$processN==0) & all(X$processA[index_strict] ==0))
}
at_risk_N <- function(X, time, ...){
  index_strict <- X$t < time
  as.numeric(all(X$processN[index_strict]==0) & all(X$processA[index_strict] ==0))
}

long_data[t <= 10, at_risk_N(.SD, 10), by = id]

at_risk_A <- Summary_measure$new(c("processA", "processN", "t"), at_risk_A)
at_risk_N <- Summary_measure$new(c("processA", "processN", "t"), at_risk_N)

npsem = list(define_node("processN", "processN", c("W", "A"), time = sort(unique(setdiff(long_data$t, 0))), risk_set_map = at_risk_N, missing_row_implies_not_at_risk = F),
             define_node("processA", "processA", c("W", "A"), time = sort(unique(setdiff(long_data$t, 0))), risk_set_map = at_risk_A, missing_row_implies_not_at_risk = F),
             define_node("A", "A", c("W"), time = 0),
             define_node("W", "W", c(), time = 0))
task <- tmle3_Task$new(long_data, npsem, t= "t", id = "id")

factor_list = list(
  LF_emp$new("W"),
  LF_fit$new("A", make_learner(Lrnr_glm)),
  LF_fit$new("processA", make_learner(Lrnr_glm), type = "mean"),
  LF_fit$new("processN", make_learner(Lrnr_glm), type = "mean")
)



task <- tmle3_Task$new(long_data, npsem, t= "t", id = "id")

lik <- Likelihood$new(factor_list)

lik_trained <- lik$train(task)



cf_task <- task$clone()
cf_task$force_at_risk <- T


cf_task$get_regression_task("processN", expand = T)$get_data()
task$get_regression_task("processN", expand = T)$get_data()
lik_trained$get_likelihood(cf_task, "processN",  to_wide = T)
lik_trained$get_likelihood(task, "processN",  to_wide =T, drop=F)
