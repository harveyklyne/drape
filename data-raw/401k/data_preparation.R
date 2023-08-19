df <- DoubleML::fetch_401k(return_type = "data.frame",
                           polynomial_features = FALSE,
                           instrument = FALSE)
data <- list("x" = df$inc,
             "z" = as.matrix(df[ , !(names(df) %in% c("net_tfa", "inc"))]))

saveRDS(data, "401k_data.rds")
