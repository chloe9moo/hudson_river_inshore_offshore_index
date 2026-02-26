## VAST MODEL RUNS :: POS CATCH DISTRIBUTION COMPARISON ##
##vast author recommends exploring multiple distributions, so this is that

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

### MODEL SETTINGS/OPTIONS ###
# num_knots = 200
do.bias.correct = FALSE
calc_range_in_models = FALSE
calc_eff_area_in_models = FALSE
run_model_opt = TRUE
compare_models = TRUE

#distribution options:
##author rec. trying multiple distributions, using conventional delta model link function (do not use default (1) link, that is for biomass only)
obs.model.m <- matrix(c(2, 4, 5, 7, 11,
                        0, 0, 0, 0, 0), nrow = 5, ncol = 2)  

##############start################
#load in species dat ----
io.dat.spp <- fread(paste0(spp.wd, "/", gsub(" ", "_", spp), "_vast_prep_dat_", survey2run, ".csv"))

#load in extrap grid ----
ext.g <- fread(paste0(PATH, "/inshore_offshore/hudson_river_VAST_extrap_grid_250cs.csv"))
ext.g <- as.data.frame(ext.g)
ext.g$Depth <- ext.g$depth ##necessary for strata limits

# model selection: settings determ. ----
settings = make_settings( n_x = 200, #number of vertices in SPDE mesh (knots)
                          Region = "user",
                          knot_method = "grid", #determine location of GMRF vertices, other option is samples
                          purpose = "index2", 
                          use_anisotropy = TRUE,
                          bias.correct = do.bias.correct,
                          fine_scale = FALSE )

## catchability ----
#gear is minimum cov
# set catchability formula(s)
Q1_formula = ~ gear_def
Q2_formula = ~ gear_def

# relevel gear as factor, 1m Tucker Trawl is base level unless bss survey
if(survey2run != "bss") {
  io.dat.spp$gear_def <- relevel(as.factor(io.dat.spp$gear_def), ref = "1m Tucker Trawl")
} else {
  io.dat.spp$gear_def <- relevel(as.factor(io.dat.spp$gear_def), ref = "100' x 10' Beach Seine")
}

catch_data <- io.dat.spp[, .(gear_def)]

## specify derived quantities ----
##can turn these off to save time on runs
settings$Options[["Calculate_Range"]] <- calc_range_in_models
settings$Options[["Calculate_effective_area"]] <- calc_eff_area_in_models

#other settings default:
#specify spatial / spatiotemporal variation params ---
##default:
# settings$FieldConfig <- matrix(rep("IID", 6), ncol = 2, nrow = 3, dimnames = list(c("Omega", "Epsilon", "Beta"), c("Component_1", "Component_2")))

#specify temporal correlation ---
#default:
# settings$RhoConfig <- c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0)

model_opt <- data.table()
starttime <- Sys.time()
for(i in seq_len(nrow(obs.model.m))) {
  
  # if(i %in% 1:4) { next }
  
  model.save.name <- paste0(gsub(" ", "_", spp), "_settingscomp_", survey2run, "_", i, "_", gsub("-", "", run_start_date), ".rds")
  
  settings$ObsModel <- obs.model.m[i,]
  
  fit <- tryCatch(
    fit_model( settings = settings,
               Lat_i = io.dat.spp$lat,
               Lon_i = io.dat.spp$long,
               t_i = io.dat.spp$year,
               b_i = as_units(io.dat.spp$ct_yoy, "count"),
               a_i = as_units(io.dat.spp$area_swept, "m^2"), #is this true? : For VAST, you typically log-transform effort or area before including it as an offset in the model.
               Q1_formula = Q1_formula,
               Q2_formula = Q2_formula,
               catchability_data = catch_data,
               input_grid = ext.g,
               run_model = run_model_opt,
               working_dir = spp.wd),
    error = function(e) e
  )
  
  if(inherits(fit, "error")) {
    error.message <- fit$message
  } else {
    error.message <- NA
    
    #save model
    saveRDS(fit, paste0(spp.wd, "/vast_models/", model.save.name))
  }
  
  tmp1 <- data.table(species = spp, 
                     obs_mod = paste(obs.model.m[i,], collapse = " "), 
                     save_name = model.save.name,
                     error_message = error.message)
  model_opt <- rbindlist(list(model_opt, tmp1))
  
  fwrite(model_opt, paste0(spp.wd, "/model_directory_settingscomp_", survey2run, "_", gsub("-", "", run_start_date), ".csv"))
  
  rm(fit, tmp1)
}
endtime <- Sys.time()
endtime-starttime

## read in models, compare ----
model_opt <- fread(paste0(spp.wd, "/model_directory_settingscomp_", survey2run, "_", gsub("-", "", run_start_date), ".csv"))

for(i in  seq_len(nrow(model_opt))) {
  
  mod.n <- model_opt[i, save_name]
  
  mod.n <- list.files(paste0(spp.wd, "/vast_models"), mod.n, full.names = T)
  
  if(length(mod.n) == 0) { next }
  
  #read in model
  fit <- readRDS(mod.n)
  # fit <- reload_model(fit)
  
  #get aic
  if(!is.null(fit$parameter_estimates$AIC)) {
    
    ##RMSE prep
    d_pred <- fit$Report$D_i
    d_obs <- strip_units(fit$data_frame$b_i / fit$data_frame$a_i)
    
    #get aic
    model_opt[i, `:=` (aic = fit$parameter_estimates$AIC,
                       #get rmse
                       rmse = sqrt(mean((d_obs - d_pred)^2)),
                       #pos def hess
                       PDH = fit$parameter_estimates$SD$pdHess,
                       #convergence check
                       convergence_check = fit$parameter_estimates$Convergence_check,
                       #final max abs gradient
                       max_grad = max(abs(fit$parameter_estimates$diagnostics$final_gradient)),
                       #parameter issue check
                       param_issue = check_fit(fit$parameter_estimates, check_gradients = TRUE)
    )]

  }
  
}

fwrite(model_opt, paste0(spp.wd, "/model_directory_settingscomp_", survey2run, "_", gsub("-", "", run_start_date), ".csv"))


#pull in and compare ALL results ----
if(FALSE) { #BUT DON'T RUN IN FULL RUN!
  
  spp.l <- c("Alewife", "Striped Bass", "American Shad", "Blueback Herring")
  dates2check <- c("20250918")
  run_location <- "amarel_cluster/" #"amarel_cluster/" #or ""
  
mod.out.all <- lapply(spp.l, function(spp) {
    
    dir.list <- list.files(paste0(PATH, "/inshore_offshore/", run_location, gsub(" ", "_", spp)), "model_directory_settingscomp", full.names = T)
    dir.list <- dir.list[grepl(paste(dates2check, collapse = "|"), dir.list)]
    
    model_opt.l <- vector("list", length = length(dir.list))
    for(i in seq_along(dir.list)) {
      
      DT.str <- dir.list[[i]]
      
      DT <- fread(DT.str)
      DT[, survey := fifelse(grepl("fjs", DT.str), "fjs",
                             fifelse(grepl("bss", DT.str), "bss",
                                     "all"))]
      
      model_opt.l[[i]] <- DT
    }
    
    model_opt.l <- rbindlist(model_opt.l, fill = TRUE)
    
    #for run mess ups that saved the model but didn't save the results
    for(i in seq_len(nrow(model_opt.l))) {
      
      if(is.na(model_opt.l[i, aic]) && (model_opt.l[i, error_message] == "" | is.na(model_opt.l[i, error_message]))) {
        
        #read in model
        fit <- readRDS(paste0(PATH, "/inshore_offshore/", run_location, gsub(" ", "_", spp), "/vast_models/", model_opt.l[i, save_name]))
        # fit <- reload_model(fit)
        
        if(!is.null(fit$parameter_estimates$AIC)) {
          
          ##RMSE prep
          d_pred <- fit$Report$D_i
          d_obs <- strip_units(fit$data_frame$b_i / fit$data_frame$a_i)
         
          #get aic
          model_opt.l[i, `:=` (aic = fit$parameter_estimates$AIC,
                               #get rmse
                               rmse = sqrt(mean((d_obs - d_pred)^2)),
                               #pos def hess
                               PDH = fit$parameter_estimates$SD$pdHess,
                               #convergence check
                               convergence_check = fit$parameter_estimates$Convergence_check,
                               #final max abs gradient
                               max_grad = max(abs(fit$parameter_estimates$diagnostics$final_gradient)),
                               #parameter issue check
                               param_issue = check_fit(fit$parameter_estimates, check_gradients = TRUE)
                               )]
           
        } else {
          
          model_opt.l[i, error_message := "Likely returned max absolute gradient of -Inf"]
          
        }
      } else {
        next ##already pulled results
      }
    } ##end run mess ups (re pull results)
    
    return(model_opt.l)
    
  })
  
mod.out.all <- rbindlist(mod.out.all, use.names = T)

mod.out.all[, `:=` (delta_aic = aic - min(aic, na.rm = T),
                    max_grad_check = max_grad < 10e-4), by = c("species", "survey")]

fwrite(mod.out.all, paste0(PATH, "/inshore_offshore/", run_location, "multi-species_omsens_results_", dates2check, ".csv"))

#select the models to use for next steps
##first diagnostic check
mod.top <- mod.out.all[max_grad_check & PDH & !param_issue & convergence_check == "There is no evidence that the model is not converged"]

##redo delta aic
mod.top[, delta_aic := aic - min(aic, na.rm = T), by = .(species, survey)]

##if models are basically indistinguishable by aic, use rmse
mod.top <- mod.top[
  delta_aic <= 4
  ,
  if (all(.N == 1)) .SD
  else .SD[.N == 1 | rmse == min(rmse)],
  by = .(species, survey)
]

if(nrow(mod.top) != length(spp.l) * 3) { stop("One model for every species x survey combo was not returned!!!") }


fwrite(mod.top, paste0(PATH, "/inshore_offshore/", run_location, "multi-species_omsens_results_filtered_", dates2check, ".csv"))

}
