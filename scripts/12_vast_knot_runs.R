## VAST MODEL RUNS :: COVARIATE COMPARISON ##
###notes for some model decisions can be found in io_vast_examples.R, VAST notes on google docs
###aim: run through model selection options we are testing for specific species and survey
library(data.table); library(VAST)
library(splines)

PATH <- getwd()
run_start_date <- gsub("-", "", Sys.Date()) ##for organizing saves

#set variables ----
### SELECT SURVEY INDEX ###
#all, fjs, or bss
survey2run = "all"

### SELECT A SPECIES TO RUN ###
spp = "Alewife"
spp.wd = paste0(PATH, "/inshore_offshore/", gsub(" ", "_", spp))

if(!dir.exists(spp.wd)) { dir.create(spp.wd) }
if(!dir.exists(paste0(spp.wd, "/vast_models"))) { dir.create(paste0(spp.wd, "/vast_models")) }
if(!dir.exists(paste0(spp.wd, "/knot_dir"))) { dir.create(paste0(spp.wd, "/knot_dir")) }

### MODEL SETTINGS/OPTIONS ###
num_knots = 200 
calc_range_in_models = FALSE
calc_eff_area_in_models = FALSE
run_model_opt = TRUE

##############start################
#load in extrap grid ----
ext.g <- fread(paste0(PATH, "/inshore_offshore/hudson_river_VAST_extrap_grid_250cs.csv"))
ext.g <- as.data.frame(ext.g)
ext.g$Depth <- ext.g$depth ##necessary for strata limits

#pull in data and settings from selected model ----
##covariates ----
#directory of 'top' selected models
mod.top <- fread(paste0(PATH, "/inshore_offshore/multi-species_cov_comp_results_filtered.csv"))
mod.top <- mod.top[species == spp & survey_set == survey2run]

top.v.m <- readRDS(paste0(spp.wd, "/vast_models/", mod.top$model_name))

##catch dat
io.dat.spp <- top.v.m$data_frame

##cov dat
catchability_data <- top.v.m$catchability_data

if(survey2run != "bss") {
  catchability_data$gear_def <- relevel(as.factor(catchability_data$gear_def), ref = "1m Tucker Trawl")
} else {
  catchability_data$gear_def <- relevel(as.factor(catchability_data$gear_def), ref = "100' x 10' Beach Seine")
}

covariate_data <- top.v.m$covariate_data

Q1_formula = top.v.m$Q1_formula
Q2_formula = top.v.m$Q2_formula
X1_formula = top.v.m$X1_formula
X2_formula = top.v.m$X2_formula

## settings ----
settings = top.v.m$settings

##can turn these off to save time on runs (doesn't have to be the same from cov run)
settings$Options[["Calculate_Range"]] <- calc_range_in_models
settings$Options[["Calculate_effective_area"]] <- calc_eff_area_in_models

#set knots for this run: ----
settings$n_x = num_knots

#run model(s) ----
#copy reduced extrap grid if it exists into species folder (creating it takes a long time)
if(file.exists(paste0(PATH, "/inshore_offshore/Kmeans_extrapolation-2000.RData")) &
   !file.exists(paste0(spp.wd, "/Kmeans_extrapolation-2000.RData"))) {
  
  file.copy(from = paste0(PATH, "/inshore_offshore/Kmeans_extrapolation-2000.RData"),
            to = paste0(spp.wd, "/Kmeans_extrapolation-2000.RData"))
}

model.save.name <- paste0(gsub(" ", "_", spp), "_", survey2run, "_knots_", num_knots, "_", run_start_date, ".rds")
    
fit <- tryCatch(
  fit_model( settings = settings,
                   Lat_i = io.dat.spp$Lat_i,
                   Lon_i = io.dat.spp$Lon_i,
                   t_i = io.dat.spp$t_i,
                   b_i = io.dat.spp$b_i,
                   a_i = io.dat.spp$a_i,
                   input_grid = ext.g,
                   #density:
                   X1_formula = X1_formula,
                   X2_formula = X2_formula,
                   covariate_data = covariate_data,
                   #catchability:
                   Q1_formula = Q1_formula,
                   Q2_formula = Q2_formula,
                   catchability_data = catchability_data,
             # bias.correct = F, #i think it already does this with index2, warnings suggests so when = T
             run_model = run_model_opt,
             working_dir = spp.wd
             ),
  error = function(e) e
)
    
if(inherits(fit, "error")) {
  error.message <- fit$message
} else {
  error.message <- NA
  
  #save model
  saveRDS(fit, paste0(spp.wd, "/vast_models/", model.save.name))
}

tmp <- data.table(species = spp, 
                  knots = num_knots, 
                  save_name = model.save.name,
                  error_message = error.message)

if(!is.null(fit$parameter_estimates$AIC)) {
  
  ##RMSE prep
  d_pred <- fit$Report$D_i
  d_obs <- strip_units(fit$data_frame$b_i / fit$data_frame$a_i)
  
  #diagnostic table output
  tmp[, `:=` (
    #convergence check
    convergence_check = fit$parameter_estimates$Convergence_check,
    PDH = fit$parameter_estimates$SD$pdHess,
    max_grad = max(abs(fit$parameter_estimates$diagnostics$final_gradient)),
    max_grad_check = max(abs(fit$parameter_estimates$diagnostics$final_gradient)) < 10e-4,
    param_issue = check_fit(fit$parameter_estimates, check_gradients = TRUE),
    #model selection metrics
    aic = fit$parameter_estimates$AIC,
    rrmse = sqrt(mean((d_obs - d_pred)^2)) / mean(d_obs),
    rmse = sqrt(mean((d_obs - d_pred)^2))
  )]
  
}

fwrite(tmp, paste0(spp.wd, "/knot_dir/model_dir_knots_", num_knots, "_", survey2run, "_", run_start_date, ".csv"))

