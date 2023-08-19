library(parallelly)
library(parallel)
library(future)
library(future.apply)
library(progressr)
library(tidyverse)

options(future.wait.interval=0L)

### Initialise progress bar
handlers(handler_progress(format="[:bar] :percent :eta :message"))


### Initialise computing cluster
worker_names <- NULL # put your worker names here
working_directory <- NULL # put your working directory here

### Wake workers
for (worker in worker_names){
  system(paste0("/alt/bin/wake ", worker))
}
Sys.sleep(60)

nodes_per_worker <- 1

my_cluster <- parallelly::makeClusterPSOCK(
  rep(worker_names, nodes_per_worker),
  outfile="", # this option ensures that error messages and status messages are printed to the console of the host (i.e. the computer the running the script)
  homogeneous=FALSE) # homogeneous = FALSE is crucial if the operating system of the host (i.e. the computer running the script) differs from the operating system of the workers


devtools::load_all()

clusterEvalQ(my_cluster, library(devtools))
clusterExport(my_cluster, c("working_directory"), envir = .GlobalEnv)
clusterEvalQ(my_cluster, setwd(working_directory))
clusterEvalQ(my_cluster, devtools::load_all())

plan(cluster, workers=my_cluster)


### Define simulations
reps <- 1000

ex_settings <- c("normal", "mixture2", "mixture3", "logistic", "t4", "401k")
f_settings <- c("plm", "additive", "interaction")

for (ex_setting in ex_settings){
  for (f_setting in f_settings){

    sim_df <- expand.grid(reps=1:reps)

    with_progress( {
      prog_bar <- progressor(along=1:(nrow(sim_df)))
      sim_res <- future_apply(sim_df, MARGIN=1, future.seed=TRUE, future.scheduling = Inf, simplify=FALSE, FUN = function(x) {
        prog_bar()

        out <- simulate_data(n = 9915, ex_setting = ex_setting, f_setting = f_setting)
        return(c(out$sample_theta, out$sample_std_dev))
      })})

    sim_res_mat <- Reduce(rbind, sim_res)

    theta <- round(mean(sim_res_mat[,1]), 4)
    if (ex_setting == "401k"){output <- list("theta" = theta)}
    else{
      std_dev <- round(mean(sim_res_mat[,2]), 4)
      output <- list("theta" = theta,
                     "std_dev" = std_dev)
    }

    jsonData <- rjson::toJSON(output)
    write(jsonData, paste0("target_", f_setting, "_", ex_setting, ".json"))

  }
}

stopCluster(my_cluster)
