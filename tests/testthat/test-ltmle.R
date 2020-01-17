library(simcausal)
library(tmle3)
options(simcausal.verbose=FALSE)
D <- DAG.empty()
D <- D +
  node("L2", t = 0, distr = "rbern",
       prob = 0.05) +
  node("L1", t = 0, distr = "rbern",
       prob = ifelse(L2[0] == 1, 0.5, 0.1)) +
  node("A1", t = 0, distr = "rbern",
       prob = ifelse(L1[0] == 1 & L2[0] == 0, 0.5,
                     ifelse(L1[0] == 0 & L2[0] == 0, 0.1,
                            ifelse(L1[0] == 1 & L2[0] == 1, 0.9, 0.5))))

t.end <- 2
D <- D +
  node("L2", t = 1:t.end, distr = "rbern",
       prob =
         ifelse(A1[t-1] == 1, 0.1,
                ifelse(L2[t-1] == 1, 0.9, min(1, 0.1 + t / 16)))) +
  node("A1", t = 1:t.end, distr = "rbern",
       prob = ifelse(A1[t-1] == 1, 0.5,
                     ifelse(L1[0] == 1 & L2[t] == 0, 0.3,
                            ifelse(L1[0] == 0 & L2[t] == 0, 0.1,
                                   ifelse(L1[0] == 1 & L2[t] == 1, 0.7, 0.5)))))+
  node("Y", t = 0:t.end, distr = "rbern",
      prob =
        plogis(-6.5 + 5 * A1[t] + L1[0] + 4 * L2[t] + 0.05 * sum(I(L2[0:t] == rep(0, t+1)))),
      EFU = FALSE)
lDAG <- set.DAG(D)
plotDAG(lDAG, xjitter = 0.3, yjitter = 0)
Odat <- sim(DAG = lDAG, n = 100, rndseed = 123)
ldat <- sim(DAG = lDAG, n = 100, rndseed = 123, wide=FALSE)
node_list <- list(W=c("L1"),
                  A=c("A1"),
                  Y="Y",
                  time="t",
                  id="ID")

npsem <- list(
  define_node("W", node_list$W),
  define_node("A", node_list$A, list("W")),
  define_node("Y", node_list$Y, list("A", "W", list(name="Y", lag=-1)))
)

# parent format, either a character = just get that node
# or a list(name="", lag=0) = get that node lagged
tmle_task <- tmle3_Task$new(ldat, npsem=npsem, nodes=list(id = "ID", time="t"))
Q_task <- tmle_task$get_regression_task("Y")

# specify npsem
# npsem <- list(
#   define_node("W", node_list$W, variable_type = variable_types$W),
#   define_node("A", node_list$A, c("W"), variable_type = variable_types$A),
#   define_node("Y", node_list$Y, c("A", "W"), variable_type = variable_types$Y, scale = TRUE)
# )

# function of nodes at previous times
node_list <- setdiff(names(Odat),"ID")
names(node_list) <- node_list
i <- 1

# assume time ordering in variables
nodes <- lapply(seq_along(node_list),function(i){
  node <- node_list[i]
  if(i==1){
    parents <- c()
  } else {
    parents <- unlist(node_list[1:(i-1)])
  }
  
  node_obj <- define_node(as.vector(node), names(node), parents)
  return(node_obj)
})

nodes[[10]]$group

# get groupped parent nodes
# make regression task with repeated measures, appropriately lagged cols
# fit task


lnpsem <- list(
  define_node("L2_0", "L2_0"),
  define_node("L1_0", "L1_0", "L2_0"),
  define_node("A1_0", "A1_0", c("L2_0","L1_0")),
  define_node("Y_1", "Y_1", c("L2_0","L1_0","A1_0")),
  define_node("L2_1", "L2_1", c("L2_0","L1_0","A1_0")),
  define_node("A", node_list$A, c("W"), variable_type = variable_types$A),
  define_node("Y", node_list$Y, c("A", "W"), variable_type = variable_types$Y, scale = TRUE)
)


