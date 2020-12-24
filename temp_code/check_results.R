library(dplyr)
library(xtable)

if_misspec <- T

if (!if_misspec) {  # correct models
  allofthem <- list.files(paste0("/home/leo42k/Documents/projects/brc/20201008_output/", "correct", "/temp_output"), full.names = T)
} else {
  allofthem <- list.files(paste0("/home/leo42k/Documents/projects/brc/20201008_output/", "misspec", "/temp_output"), full.names = T)
}

# allofthem <- list.files("/home/leo42k/Documents/projects/20201006_output/misspec/temp_output/", full.names = T)
# allofthem <- list.files("/home/leo42k/Documents/projects/20201006_output/correct/temp_output/", full.names = T)

n_total_finished <- allofthem %>% length
XX <- do.call(rbind, lapply(allofthem, function(x) readRDS(x) %>% unlist))

XX[, 1]

library(R6)
code_list <- list.files("./R", full.names = T)
for (code in code_list) source(code)
source("./temp_code/generate_data.R")


timepoint <- 2
data_truth <- generate_Zheng_data(B = 100000, tau = timepoint, seed = 202008, setAM = c(1, 0), if_LY_misspec = if_misspec)
truth <- data_truth[[timepoint + 1]]$Y %>% mean
truth

list_mse <- lapply(1:6, function(i) {
  mean((XX[, i*3 - 2] - truth)^2)
})
names(list_mse) <- colnames(XX)[(1:6)*3 - 2]
list_mse

list_bias <- lapply(1:6, function(i) {
  abs(XX[, i*3 - 2] - truth) %>% mean
})

list_sd <- lapply(1:6, function(i) {
  sd(XX[, i*3 - 2])
})

list_coverage <- lapply(1:6, function(i) {
  (truth < XX[, i*3] & truth > XX[, i*3 - 1]) %>% mean
})

report <- data.frame(MSE = list_mse %>% unlist,
                     Bias = list_bias %>% unlist,
                     SD = list_sd %>% unlist,
                     Coverage = list_coverage %>% unlist)

rownames(report) <- c("No targeting",
                      "First-step logistic",
                      "Iterative logistic, maxit 10",
                      "One-step EIC, maxit 100, fixed direction",
                      "One-step EIC, maxit 100, diff direction",
                      "Projection-based, one-step EIC, maxit 100")
report

report %>%
  xtable(type = "latex",
         caption = paste0("Sample size ", 1000, "; iteration: ", n_total_finished),
         digits = 6) %>%
  print(caption.placement = "top",
        file = paste0("./temp_output/", 1000, "_LY_", ifelse(if_misspec, "misspec", "correct"), "_glm_sample_size20201008_.tex")
        )
