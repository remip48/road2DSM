##--------------------------------------------------------------------------------------------------------
## Authors : Remi Pigeault, Matthieu Authier
## Last update : 2025-03-30
## R version 4.4.0 (2024)
## Copyright (C) 2024 The R Foundation for Statistical Computing
## Platform: x86_64-w64-mingw32/x64 (64-bit)
##--------------------------------------------------------------------------------------------------------

beta <- function(fitted_model, n_sim = 200, seed_id = 19790428) {
  if(!is.null(seed_id)) {
    set.seed(seed_id)
  } else {
    set.seed(sample.int(1e6, size = 1))
  }
  theta <- mvtnorm::rmvnorm(n_sim, mean = fitted_model$coefficients, sigma = fitted_model$Vp)
  return(theta)
}

n_cores <- floor(nc_cores / 4)
cat("Parallel processing with", n_cores, "for BIAS estimate and correction.\n")

cl <- makeCluster(n_cores)
registerDoParallel(cl)

theta <- beta(fitted_model = gam_mod, n_sim = 1e3)

allX <- foreach(X = to_pred,
                .packages = c("mgcv"),
                .noexport = ls()[!(ls() %in% c("gam_mod"))]) %dopar% {
                  return(predict.gam(gam_mod, newdata = X, type = "lpmatrix"))
                }

stopCluster(cl)
gc()

cl <- makeCluster(floor(detectCores() / 4))
registerDoParallel(cl)

allpred <- foreach(l = allX,
                   .noexport = ls()[!(ls() %in% c("theta"))]) %dopar% { return(exp(l %*% t(theta))) }


stopCluster(cl)
gc()
