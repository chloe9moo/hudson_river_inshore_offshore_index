## functions for use across scripts ##

#progress function for for loop
print_progress <- function(iter_num, total_iter) {
  percent_complete <- round((iter_num / total_iter) * 100, 2)
  cat(sprintf("\rloop progress: %0.2f%% complete", percent_complete))
}


#pull out abundance index
get_index <- function(model.save.location, model.save.name, vast.mod = NULL, index_name = "Index_ctl") {
  
  if(is.null(vast.mod)) {
    #read in model
    vast.mod <- readRDS(paste(model.save.location, model.save.name, sep = "/"))
    # vast.mod <- reload_model(vast.mod)
  }
  
  sdreport <- vast.mod$parameter_estimates$SD
  
  par_hat <- TMB:::as.list.sdreport( sdreport, what="Estimate", report=TRUE )
  par_SE = TMB:::as.list.sdreport( sdreport, what="Std. Error", report=TRUE )
  
  #just in case, will need to double check later if adding in seasonal time steps!!
  time.dat <- data.table(o = vast.mod$data_list$t_i,
                         time = vast.mod$data_frame$t_i) |>
    unique()
  time.dat <- time.dat[order(o)]
  
  s.names <- vast.mod$settings$strata.limits$STRATA
  
  if(length(s.names) == 1) {
    dt <- data.table(estimate = as.vector(par_hat[[index_name]]),
                     std_error_for_estimate = as.vector(par_SE[[index_name]]),
                     time = as.vector(time.dat$time))
    
    dt[, `:=` (norm_estimate = estimate / mean(estimate),
               norm_std_error = std_error_for_estimate / mean(estimate))]
    
    dt[, `:=` (upper = estimate + 1.96 * std_error_for_estimate,
               lower = estimate - 1.96 * std_error_for_estimate,
               norm_upper = norm_estimate + 1.96 * norm_std_error,
               norm_lower = norm_estimate - 1.96 * norm_std_error)]
    
  } else {
    dt <- data.table(strata = c(rep(s.names[[1]], length(time.dat$time)), rep(s.names[[2]], length(time.dat$time)), rep(s.names[[3]], length(time.dat$time))),
                     estimate = as.vector(par_hat[[index_name]]),
                     std_error_for_estimate = as.vector(par_SE[[index_name]]),
                     time = rep(as.vector(time.dat$time), 3))
    
    dt[, mn_estimate := mean(estimate), by = strata]
    
    dt[, `:=` (norm_estimate = estimate / mn_estimate,
               norm_std_error = std_error_for_estimate / mn_estimate)]
    
    dt[, `:=` (upper = estimate + 1.96 * std_error_for_estimate,
               lower = estimate - 1.96 * std_error_for_estimate,
               norm_upper = norm_estimate + 1.96 * norm_std_error,
               norm_lower = norm_estimate - 1.96 * norm_std_error)]
  }
  
  return(dt)
}