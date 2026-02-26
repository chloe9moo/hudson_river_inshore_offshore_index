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

### SET COV OPTIONS ###
##catchability options: gear, time (diff from solar noon), sample depth
catch.opt <- c(~0,
               ~ gear_def, 
               ~ gear_def + bs( log(sam_dpth + 1), df = 3),
               ~ gear_def + ns(solar_noon_diff, df = 5),
               ~ gear_def + ns(solar_noon_diff, df = 5) + bs( log(sam_dpth + 1), df = 3)
)

##density options: river depth, bottom type
cov.opt <- c(~ 0,
             ~ bs( log(riv_dpth + 1), df = 3),
             ~ factor(benth_cat),
             ~ bs( log(riv_dpth + 1), df = 3) + factor(benth_cat)
)

#skip any covariate sets?
skip.catch = FALSE
skip.dens = FALSE

### MODEL SETTINGS/OPTIONS ###
include_strata = FALSE
# num_knots = 200
do.bias.correct = TRUE 
calc_range_in_models = FALSE
calc_eff_area_in_models = FALSE
run_model_opt = TRUE

##############start################
# species + env data prep ----
##load in cleaned/filtered data for target species x survey (see 02_prep_bio_dat.R)
io.dat.spp <- fread(paste0(spp.wd, "/", gsub(" ", "_", spp), "_vast_prep_dat_", survey2run, ".csv"))

##prepare covariate columns
###relevel gear as factor, 1m Tucker Trawl is base level
if(survey2run != "bss") {
  io.dat.spp$gear_def <- relevel(as.factor(io.dat.spp$gear_def), ref = "1m Tucker Trawl")
} else {
  io.dat.spp$gear_def <- relevel(as.factor(io.dat.spp$gear_def), ref = "100' x 10' Beach Seine")
}

###pull covariate col names (some need to be standardized)
get_cov_col_names <- function(x) {
  cn <- strsplit(gsub("~|\\+|log|\\(|\\)|poly|ns|bs|factor|1|3|5|2|4|\\=|\\,|df|degree", "", paste(x, collapse = " ")), split = " ")[[1]]
  cn <- unique(cn[!cn %in% c("", "0")])
}

catch.names <- get_cov_col_names(catch.opt)
cov.names <- get_cov_col_names(cov.opt)

##double check for missing data
if(any(io.dat.spp[, is.na(.SD), .SDcols = intersect(names(io.dat.spp), c(catch.names, cov.names))])) { stop("CHECK SPECIES DATA, WHY ARE THERE NAs IN THE COV COLS?\n") }

###from VAST wiki: rescale covariates being used to have an SD > 0.1 and < 10 (for numerical stability)
stand_cols <- names(which(sapply(io.dat.spp[, .SD, .SDcols = c(catch.names, cov.names)], is.numeric)))
io.dat.spp[, names(.SD) := lapply(.SD, function(x) x / 10), .SDcols = stand_cols]

###VAST cov matrix requires specific names
setnames(io.dat.spp, c("lat", "long"), c("Lat", "Lon"))

#load in extrap grid ----
ext.g <- fread(paste0(PATH, "/inshore_offshore/hudson_river_VAST_extrap_grid_250cs.csv"))
ext.g <- as.data.frame(ext.g)
ext.g$Depth <- ext.g$depth ##necessary for strata limits

#set settings ----
settings = make_settings( n_x = 200, #number of vertices in SPDE mesh (knots)
                          Region = "user",
                          knot_method = "grid", #determine location of GMRF vertices, other option is samples
                          purpose = "index2", 
                          use_anisotropy = TRUE,
                          bias.correct = do.bias.correct,
                          fine_scale = FALSE)

#using :
top.mod <- fread(paste0(PATH, "/inshore_offshore/multi-species_omsens_results_filtered_20250915.csv"))
om <- top.mod[species == spp & survey == survey2run]

if(nrow(om) == 0) { stop("no obs model selected for ", spp, " and ", survey2run) }

om <- om$obs_mod
om <- strsplit(om, split = " ")[[1]]

settings$ObsModel[1] <- as.numeric(om[1])
settings$ObsModel[2] <- 0 ##always 0

#set strata ---
if(include_strata & survey2run == "all") {
  # range(ext.g$Depth)
  # ifelse(region_grid.df$depth < 3.048, "inshore", "offshore")
  settings$strata.limits <- data.frame('STRATA' = as.factor(c("All", "inshore", "offshore")),
                                       'shallow_border' = c(0, 0, 3.05),
                                       'deep_border' = c(55, 3.05, 55))
}

##can turn these off to save time on runs
settings$Options[["Calculate_Range"]] <- calc_range_in_models
settings$Options[["Calculate_effective_area"]] <- calc_eff_area_in_models

#run model(s) ----
#copy reduced extrap grid if it exists into species folder (creating it takes a long time)
if(file.exists(paste0(PATH, "/inshore_offshore/Kmeans_extrapolation-2000.RData")) &
   !file.exists(paste0(spp.wd, "/Kmeans_extrapolation-2000.RData"))) {
  
  file.copy(from = paste0(PATH, "/inshore_offshore/Kmeans_extrapolation-2000.RData"),
            to = paste0(spp.wd, "/Kmeans_extrapolation-2000.RData"))
}

##make dataframe for saving model options used in set up
model_opt <- data.table()
starttime <- Sys.time()
i <- 4; y <- 2 ##testing
for(i in seq_along(catch.opt)) {
  
  #skip some set ups for reruns
  if(skip.catch) {
    if(i %in% c(10)) { next } 
  }
  
  # set catchability formula(s)
  Q1_formula = catch.opt[[i]]
  Q2_formula = catch.opt[[i]]
  
  if(i == 1) {
    catchability_data <- NULL
  } else {
    col_sub <- catch.names[sapply(catch.names, function(x) grepl(x, paste(catch.opt[[i]], collapse = "")))]
    catchability_data <- io.dat.spp[, ..col_sub]
  }
  
 for(y in seq_along(cov.opt)) {
    
    # if(survey2run == "bss" & y %in% c(2, 4)) { next } ##do not run models with river depth for bss
    if(skip.dens) {
      if(y %in% c(10)) { next }
    }
    
    model.save.name <- paste0(gsub(" ", "_", spp), "_", survey2run, "_", i, "_", y, "_", run_start_date, ".rds")
    
    # set density formula(s)
    X1_formula = cov.opt[[y]]
    X2_formula = cov.opt[[y]]
    
    if(y == 1) {
      covariate_data <- NULL
    } else {
      col_sub <- cov.names[sapply(cov.names, function(x) grepl(x, paste(cov.opt[[y]], collapse = "")))]
      col_sub <- c("Lon", "Lat", col_sub)
      covariate_data <- io.dat.spp[, ..col_sub]
      covariate_data[, Year := NA]
    }
    
    fit <- tryCatch(
      fit_model( settings = settings,
                       Lat_i = io.dat.spp$Lat,
                       Lon_i = io.dat.spp$Lon,
                       t_i = io.dat.spp$year,
                       b_i = as_units(io.dat.spp$ct_yoy, "count"),
                       a_i = as_units(io.dat.spp$area_swept, "m^2"), #is this true? : For VAST, you typically log-transform effort or area before including it as an offset in the model.
                       input_grid = ext.g,
                       #density:
                       X1_formula = X1_formula,
                       X2_formula = X2_formula,
                       covariate_data = covariate_data,
                       #catchability:
                       # e_i = as.numeric(io.dat.spp$gear_def)-1, #-1 so factors start at 0 not 1
                       Q1_formula = Q1_formula,
                       Q2_formula = Q2_formula,
                       catchability_data = catchability_data,
                 bias.correct = do.bias.correct, #i think it already does this with index2, warnings suggests so when = T
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

    tmp1 <- data.table(species = spp, 
                       catch_cov = paste(catch.opt[[i]], collapse = ""), dens_cov = paste(cov.opt[[y]], collapse = ""), 
                       save_name = model.save.name,
                       error_message = error.message)
    model_opt <- rbindlist(list(model_opt, tmp1))
    
    fwrite(model_opt, paste0(spp.wd, "/model_directory_", survey2run, "_", run_start_date, ".csv"))
    
    rm(list = ls()[ls() %in% c("fit", "tmp1", "covariate_data", "X1_formula", "X2_formula")])
  }
  
  rm(list = ls()[ls() %in% c("catchability_data", "Q1_formula", "Q2_formula")])
}

endtime <- Sys.time()

(endtime - starttime)/20 ## ~ 30 min per model (on average)

# #to reload where necessary
# fit <- readRDS(paste0(PATH, "/inshore_offshore/Alewife/vast_model.rds"))
# fit <- reload_model(fit)
# 
# #get generic output plots ---
# plot(fit, working_dir = paste0(spp.wd, "/vast_base_plots/plots_20250529_stratatest"))
#get index estimates ---
# Sdreport <- fit$parameter_estimates$SD
# par_hat <- TMB:::as.list.sdreport( Sdreport, what="Estimate", report=TRUE )
# est <- as.vector(par_hat[["Index_ctl"]])
# ##quick check plot
# plot(seq(1, length(est), 1), est)
# lines(seq(1, length(est), 1), est)
