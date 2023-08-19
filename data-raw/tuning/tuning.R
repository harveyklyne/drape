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

#####

ex_settings <- "401k" # c("normal", "mixture2", "mixture3", "logistic", "t4")
f_settings <- c("plm", "additive", "interaction")

n_tr <- 7932 # 800
n_te <- 1983 # 1000
reps <- 100

for (ex_setting in ex_settings){

  basis <- tune_basis(ex_setting = ex_setting,
                      reps = reps,
                      n_tr = n_tr,
                      n_te = n_te)

  xgb_m <- tune_xgb(f_setting = "m",
                    ex_setting = ex_setting,
                    reps = reps,
                    n_tr = n_tr,
                    n_te = n_te)

  spl <- tune_spline(ex_setting = ex_setting,
                        xgb_params = xgb_m,
                        reps = reps,
                        n_tr = n_tr,
                        n_te = n_te)

  stddev <- tune_stddev(ex_setting = ex_setting,
                        reps = reps,
                        n = n_tr)

  difference <- stddev / 4
  bw <- exp(seq(-3, 2, 0.1))/(2*sqrt(3)) * stddev

  rothenhausler_m <- tune_rothenhausler(ex_setting = ex_setting,
                                      f_setting = "m",
                                      reps = reps,
                                      n_tr = n_tr,
                                      n_te = n_te)

  for (f_setting in f_settings){

    xgb_f <- tune_xgb(f_setting= f_setting,
                      ex_setting = ex_setting,
                      reps = reps,
                      n_tr = n_tr,
                      n_te = n_te,
                      partially_linear = FALSE)

    xgb_f_pl <- tune_xgb(f_setting= f_setting,
                      ex_setting = ex_setting,
                      reps = reps,
                      n_tr = n_tr,
                      n_te = n_te,
                      partially_linear = TRUE)

    rothenhausler_f <- tune_rothenhausler(ex_setting = ex_setting,
                                          f_setting = f_setting,
                                          reps = reps,
                                          n_tr = n_tr,
                                          n_te = n_te)

    output <- list("basis" = basis,
                   "rothenhausler_m" = rothenhausler_m,
                   "xgb_m" = xgb_m,
                   "spl" = spl,
                   "stddev" = stddev,
                   "difference" = difference,
                   "rothenhausler_f" = rothenhausler_f,
                   "xgb_f" = xgb_f,
                   "xgb_f_pl" = xgb_f_pl,
                   "reps" = reps,
                   "n_tr" = n_tr,
                   "n_te" = n_te)

    jsonData <- toJSON(output)
    write(jsonData, paste0("tune_", f_setting, "_", ex_setting, "_results.json"))

  }

}

stopCluster(my_cluster)



