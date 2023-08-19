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

ns <- 9915 # 1000
ex_settings <- "401k" # c("normal", "mixture2", "mixture3", "logistic", "t4")
f_settings <- c("plm", "additive", "interaction")

param_grid <- expand.grid(n=ns,
                          ex=ex_settings,
                          f=f_settings)
sim_df <- dplyr::slice(param_grid, rep(1:n(), each=reps))
sim_df$rep <- rep(1:reps, nrow(param_grid))


with_progress( {
  prog_bar <- progressor(along=1:(nrow(sim_df)))
  sim_res <- future_apply(sim_df, MARGIN=1, future.seed=TRUE, future.scheduling = Inf, simplify=FALSE, FUN = function(x) {
    prog_bar()
    n <- as.numeric(x["n"])
    ex <- paste(x["ex"])
    f <- paste(x["f"])

    out <- suppressWarnings(compare_rothenhausler(n=n,
                                                     ex_setting=ex,
                                                     f_setting=f))
    return(c("sample_theta" = out$sample_theta, "est" = out$est, "se" = out$se))
  })
})


stopCluster(my_cluster)

sim_res_df <- cbind(sim_df, data.frame(Reduce(rbind, sim_res), row.names = NULL))
write.csv(sim_res_df, "simulation_rothenhausler_results.csv", row.names=FALSE)
