### COMPARE VAST MODELS ###

library(data.table); library(VAST)

PATH <- getwd()

##### SET VARS #####
multi_spp.l = c("Alewife", "American Shad", "Blueback Herring", "Striped Bass")
model.run.type = "cov comp"#"cov comp" #"knot comp"
model.run.date = "20250919"
model.run.location = "amarel_cluster/" #""
survey2check = c("all", "fjs", "bss")
#-###################-#

#pull in model directories ----
for(spp in multi_spp.l) {
  
  cat("Starting comparison for ", spp, "\n")
  
  spp.wd = paste0(PATH, "/inshore_offshore/", model.run.location, gsub(" ", "_", spp))
  
  #load in directory outputs ----
  model.dir <- lapply(survey2check, function(x) {
    
    if(grepl("cov", model.run.type)) {
      tmp <- list.files(path = spp.wd, pattern = paste0("model_directory_", x, "_(", model.run.date, ")"), full.names = T)
    } else {
      if(grepl("knot", model.run.type)) {
        tmp <- list.files(path = paste0(spp.wd, "/knot_dir"), 
                          pattern = paste0(x, "_(", model.run.date, ")"), full.names = T)
      } else {
        stop("Model run type not recognized (you probably typed it wrong :( )")
      }
    }
    
    if(length(tmp) < 1) { stop("Something is wrong pulling the model directory file, double check the string and what is available.") }
    
    tmp.l <- vector("list", length(tmp))
    for(i in 1:length(tmp)) {
      tmp.l[[i]] <- fread(tmp[i])
      tmp.l[[i]][, survey_set := x]
    }
    
    return(rbindlist(tmp.l, fill = TRUE))
  })
  
  model.dir <- rbindlist(model.dir, use.names = TRUE, fill = TRUE)
  
  #get model diagnostics + settings ----
  model.list <- model.dir[is.na(error_message) | error_message == "", save_name]
  res <- vector("list", length = length(model.list))
  for(i in 1:length(model.list)) {
    
    if(FALSE) { mod.name <- model.list[1] } #testing
    
    mod.name <- model.list[i]
    
    ##read in model
    vast.mod <- readRDS(paste0(spp.wd, "/vast_models/", mod.name))
    # vast.mod <- reload_model(vast.mod)
    
    ##get diagnostics
    ###the knot run basically already does this, but can double check with a reload
    tmp <- data.table()
    
    if(length(vast.mod$Report) == 1) {
      
      tmp[, `:=` (model_name = mod.name,
                  convergence_check = vast.mod$Report)]
      
    } else {
      
      ##RMSE prep
      d_pred <- vast.mod$Report$D_i
      d_obs <- strip_units(vast.mod$data_frame$b_i / vast.mod$data_frame$a_i)
      
      #diagnostic table output
      tmp[, `:=` (
        model_name = mod.name,
        #convergence check
        convergence_check = vast.mod$parameter_estimates$Convergence_check,
        PDH = vast.mod$parameter_estimates$SD$pdHess,
        max_grad = max(abs(vast.mod$parameter_estimates$diagnostics$final_gradient)),
        max_grad_check = max(abs(vast.mod$parameter_estimates$diagnostics$final_gradient)) < 10e-4,
        param_issue = check_fit(vast.mod$parameter_estimates, check_gradients = TRUE),
        #covariate usefulness
        deviance = vast.mod$Report$deviance,
        #model selection metrics
        aic = vast.mod$parameter_estimates$AIC,
        rrmse = sqrt(mean((d_obs - d_pred)^2)) / mean(d_obs),
        rmse = sqrt(mean((d_obs - d_pred)^2)),
        #settings
        obs_mod = paste(vast.mod$data_list$ObsModel_ez, collapse = " "),
        catch_cov = paste0(vast.mod$Q1_formula)[2],
        dens_cov = paste0(vast.mod$X1_formula)[2],
        knots = vast.mod$settings$n_x
      )]
      
    }
    
    res[[i]] <- tmp
    
    cat(sprintf("\r%-20s", sprintf("%.0f%% complete...", 100 * i / length(model.list))))
    
    if(i == length(model.list)) { cat("\nDone!\n")}
    
  }
  
  res <- rbindlist(res, fill = TRUE)
  
  if(grepl("knot", model.run.type)) {
    model.output <- res[model.dir[, .(species, save_name, error_message, survey_set, knots)], on = c("model_name" = "save_name", "knots")]
  } else {
    model.output <- res[model.dir, on = c("model_name" = "save_name")]
  }
  
  
  ##calc comparison values
  # nocov.dev <- model.output[catch_cov == "~0" & dens_cov == "~0", deviance, by = survey_set]
  model.output[, `:=` (#pct_dev_exp = round((1 - (deviance/nocov.dev))*100, 3), 
    delta_aic = aic - min(aic, na.rm = T)),
    by = survey_set] 
  
  #format output results + save ----
  model.output <- model.output[order(survey_set, delta_aic, rmse)]
  setcolorder(model.output, c("species", "catch_cov", "dens_cov", "delta_aic", "aic", "rmse", "deviance", "max_grad", "PDH", "convergence_check", "param_issue"))
  
  fwrite(model.output, paste0(spp.wd, "/model_directory_w_results_combined_", 
                              gsub(" ", "_", model.run.type), "_", gsub("\\|", "_", model.run.date), ".csv"))
  
  rm(spp.wd, model.dir, model.output, res, tmp, vast.mod, d_obs, d_pred, i, mod.name, model.list)
}


#multi species summary ----
if(length(multi_spp.l) > 1) {
  
  ##load in all model results across species
  model.output <- lapply(multi_spp.l, function(spp) {
    spp.wd <- paste0(PATH, "/inshore_offshore/", model.run.location, gsub(" ", "_", spp))
    
    #load in results
    model.output <- list.files(spp.wd, 
                               pattern = paste0("model_directory_w_results_combined_", 
                                                gsub(" ", "_", model.run.type), "_", gsub("\\|", "_", model.run.date)), 
                               full.names = T)
    if(length(model.output) != 1) { stop("Rethink pulling in the model results for ", spp, ".") }
    model.output <- fread(model.output)
    
    return(model.output)
  })
  
  model.output <- rbindlist(model.output)
  
  fwrite(model.output, paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_", gsub(" ", "_", model.run.type), "_results.csv"))
  
  #select the models to use for next steps
  ##first diagnostic check
  mod.top <- model.output[max_grad_check & PDH & !param_issue & convergence_check == "There is no evidence that the model is not converged"]
  
  ##redo delta aic (sometimes the lowest aic has diagnostic issues)
  mod.top[, delta_aic := aic - min(aic, na.rm = T), by = .(species, survey_set)]
  
  ##if models are basically indistinguishable by aic, use rmse
  mod.top <- mod.top[
    delta_aic <= 4
    ,
    if (all(.N == 1)) .SD
    else .SD[.N == 1 | rmse == min(rmse)],
    by = .(species, survey_set)
  ]
  
  if(nrow(mod.top) != length(multi_spp.l) * length(survey2check)) { stop("One model for every species x survey combo was not returned!!!") }
  
  fwrite(mod.top, paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_", gsub(" ", "_", model.run.type), "_results_filtered.csv"))
  
}





