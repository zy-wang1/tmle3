# specify task name
identifier <- "cvtmle"
new_folder <- paste0("./parallel_", identifier)

# specify the template code
raw_code <- readLines("./temp_code/simulation_cvtmle-20201027.R")

# specify wd on cluster (each date should have a different folder) in the first line
raw_code <- c("setwd(\"~/20201027/tmle3\")",
              raw_code)

# Specify total number of tasks to run
set.seed(123)
K <- 100
vec_seeds <- sample(K*K, K)

# generated task codes are saved locally too
if (!dir.exists(new_folder)) dir.create(new_folder)

# set.seed according to generated vec_seeds
# let output be saved as newfolder/temp_output/m.RDS where newfolder = ./parallel_task_name
for (m in 1:K) {
  temp_code <- raw_code
  temp_code[grep("set.seed", raw_code)] <- paste0("set.seed(", vec_seeds[m], ")")
  temp_code[grep("temp_output", raw_code)] <- c(paste0("if (!dir.exists(\"", new_folder, "/temp_output\")) dir.create(\"", new_folder, "/temp_output\")"),
                                                paste0("saveRDS(temp_result, \"", new_folder, "/temp_output/", m, ".RDS\")"))
  writeLines(temp_code, paste0(new_folder, "/", m, ".R"))
}

# generate the needed study.txt
writeLines(as.character(1:K), paste0(new_folder, "/study.txt"))

