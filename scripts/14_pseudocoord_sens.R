## psuedo coord sensitivity analysis ##
#see how much creating substitute coordinates effects the results

library(data.table); library(VAST)
library(splines)
library(ggplot2)

PATH <- getwd()
target.species <- c("Alewife", "American Shad", "Blueback Herring", "Striped Bass")
model.run.location = "/amarel_cluster"

#sections to run:
prep_data = FALSE
iter = FALSE
compare_results = TRUE

#prep catch data ----
if(prep_data) {

  cov.col.names <- c("gear_def", "sam_dpth", "solar_noon_diff", "riv_dpth", "benth_cat")
  
  ##(see 02_prep_bio_dat.R, unfortunately I didn't save a unique dataset for with 
  ##and without pseudocoords used, so some of the steps need to be repeated)
  ##load in cleaned/filtered data
  io.dat <- fread(paste0(PATH, "/inshore_offshore/combined_surveys_fjs_bss.csv"))
  ##subset years to study temporal extent
  io.dat <- io.dat[year <= 2014 & year >= 1980] ##most consistent overlap gears x years
  ##set FJS tows with no assigned river depth (coordinate errors) to be pseudocoords
  io.dat <- io.dat[(survey == "fjs" & flag_coord_issue == TRUE), `:=` (lat = NA, long = NA)]
  ##subset individual survey types
  io.dat.l <- split(io.dat, io.dat$survey)
  
  io.dat.l[[3]] <- io.dat
  names(io.dat.l)[3] <- "all"
  
  target.surveys <- names(io.dat.l)
  
  #apply data filters to each species for each survey set
  summ_datasets <- data.table()
  
  for(s in target.species) {
    
    summ_datasets.tmp <- lapply(io.dat.l, function(DT, spp = s) {
      
      if(FALSE) {
        DT = copy(io.dat.l[[3]])
        spp = "American Shad"
      }
      
      spp.wd <- paste0(PATH, "/inshore_offshore/", gsub(" ", "_", spp))
      
      ##get all sample efforts ----
      all.samples <- unique(DT$uniq_id)
      
      ##subset data to only include target species (and desired columns) ----
      DT.spp <- DT[common_name == spp,
                   .(uniq_id, common_name, ct_yoy, riv_mile, survey, year, date, time, area_swept,
                     gear_def, lat, long, strata_new, volume, riv_dpth, sam_dpth, 
                     benth_cat, benth_cat_psc, solar_noon_diff,
                     month, jul_day, jul_week, lat_psc, long_psc)]
      
      if(nrow(DT.spp) == 0) { stop("Something went wrong setting the target species, check spelling?") }
      
      ## add in zero catch counts ----
      ##identify tows/hauls without target species
      zero.samples <- all.samples[!all.samples %in% DT.spp$uniq_id]
      
      zero.samples <- DT[uniq_id %in% zero.samples, unique(.SD), .SDcols = names(DT.spp)[!names(DT.spp) %in% c("common_name", "ct_yoy")]]
      
      zero.samples[, `:=` (common_name = spp, ct_yoy = 0)]
      
      ## combine zeros and species dat
      DT.spp <- rbindlist(list(DT.spp, zero.samples), use.names = TRUE)
      
      ##fix tows with more than 1 row per species ----
      if(length(unique(DT.spp$uniq_id)) != nrow(DT.spp)) {
        
        keep_cols <- setdiff(names(DT.spp), c("uniq_id", "ct_yoy"))  # everything else
        
        dup_rows <- DT.spp[DT.spp[, .I[.N > 1], by = uniq_id]$V1]
        
        dup_rows <- dup_rows[, c(
          list(ct_yoy = sum(ct_yoy)), 
          lapply(.SD, first)            
        ), by = uniq_id, .SDcols = keep_cols]
        
        DT.spp <- rbindlist(list(DT.spp[!uniq_id %in% dup_rows$uniq_id], dup_rows), use.names = TRUE)
        
      }
      
      ##remove obs if missing data ----
      missing_dat_rows <- DT.spp[!complete.cases(DT.spp[, ..cov.col.names])] # returns rows with at least one NA in cov columns
      DT.spp <- na.omit(DT.spp, cols = cov.col.names)
      
      ## remove pseudo coords for missing lat long // extracted vars ----
      missing.coords <- DT.spp[is.na(lat)]
      DT.spp <- DT.spp[!is.na(lat)]
      
      ## update benthic categories ----
      ## benthic category causes errors in the model, need to combine rare levels
      #just doing three cats
      if(any(grepl("benth_cat", cov.col.names))) {
        DT.spp[, benth_cat := fifelse(benth_cat %in% c("sandy mud", "gravelly mud", "mud"), "mud", 
                                      fifelse(benth_cat %in% c("gravelly sand", "muddy sand", "sand"), "sand",
                                              fifelse(benth_cat %in% c("sandy gravel", "muddy gravel", "gravel"), "gravel",
                                                      "other")))] ##none should be in this category otherwise something is missing
        DT.spp[, benth_cat_psc := fifelse(benth_cat %in% c("sandy mud", "gravelly mud", "mud"), "mud", 
                                          fifelse(benth_cat %in% c("gravelly sand", "muddy sand", "sand"), "sand",
                                                  fifelse(benth_cat %in% c("sandy gravel", "muddy gravel", "gravel"), "gravel",
                                                          "other")))] ##none should be in this category otherwise something is missing
        
        
        #check for too few per group again
        min_n <- 10
        n_group <- DT.spp[, .N, by = benth_cat][order(-N)]
        
        #get small groups
        small_groups <- n_group[N < min_n]$benth_cat
        
        benth_check <- length(small_groups) > 0
        
        while(length(small_groups) != 0) {
          if(length(small_groups) == 1) {
            
            next_grp <- n_group$benth_cat[nrow(n_group)-1]
            
            warning("Making new benthic category: ", paste(c(next_grp, small_groups), collapse = "-"))
            
            DT.spp[benth_cat %in% c(small_groups, next_grp), benth_cat := paste(c(next_grp, small_groups), collapse = "-")]
            
          } else {
            
            warning("Making new benthic category: ", paste(small_groups, collapse = "-"))
            DT.spp[benth_cat %in% small_groups, benth_cat := paste(small_groups, collapse = "-")]
            
          } 
          
          #check again
          n_group <- DT.spp[, .N, by = benth_cat][order(-N)]
          
          #get small groups
          small_groups <- n_group[N < min_n]$benth_cat
        }
        
        ##now do psc benth
        n_group <- DT.spp[, .N, by = benth_cat_psc][order(-N)]
        
        #get small groups
        small_groups <- n_group[N < min_n]$benth_cat_psc
        
        benth_check <- length(small_groups) > 0
        
        while(length(small_groups) != 0) {
          if(length(small_groups) == 1) {
            
            next_grp <- n_group$benth_cat_psc[nrow(n_group)-1]
            
            warning("Making new benthic category: ", paste(c(next_grp, small_groups), collapse = "-"))
            
            DT.spp[benth_cat_psc %in% c(small_groups, next_grp), benth_cat_psc := paste(c(next_grp, small_groups), collapse = "-")]
            
          } else {
            
            warning("Making new benthic category: ", paste(small_groups, collapse = "-"))
            DT.spp[benth_cat_psc %in% small_groups, benth_cat_psc := paste(small_groups, collapse = "-")]
            
          } 
          
          #check again
          n_group <- DT.spp[, .N, by = benth_cat_psc][order(-N)]
          
          #get small groups
          small_groups <- n_group[N < min_n]$benth_cat_psc
        }
      }
      
      surv <- unique(DT.spp$survey)
      
      if(length(surv) == 2) {
        surv = "all"
      }
      
      if(length(surv) != 1) { stop("Too many surveys returned or not enough .... double check code") }
      
      fwrite(DT.spp, paste0(spp.wd, "/", gsub(" ", "_", spp), "_pc_sens_vast_prep_dat_", surv, ".csv"))
      
      #summarize final counts
      tmp <- data.table(species = spp,
                        survey = surv,
                        total_tows = nrow(DT.spp),
                        mn_pos_catch = mean(DT.spp[ct_yoy > 0]$ct_yoy),
                        total_zero_catch = nrow(DT.spp[ct_yoy == 0]),
                        tows_per_fjs = nrow(DT.spp[survey == "fjs"]),
                        tows_per_bss = nrow(DT.spp[survey == "bss"]),
                        rows_removed_for_na = nrow(missing_dat_rows),
                        rows_removed_for_missing_coords = nrow(missing.coords),
                        collapsed_benth_cat = benth_check)
      
      #return summarized
      return(tmp)
      
    })
    
    summ_datasets.tmp <- rbindlist(summ_datasets.tmp)
    
    summ_datasets <- rbindlist(list(summ_datasets, summ_datasets.tmp))
    
  }
  
  fwrite(summ_datasets, paste0(PATH, "/inshore_offshore/species_datasets_final_counts_for_model_input_psc_sens.csv"))
  
  rm(io.dat, io.dat.l, summ_datasets, summ_datasets.tmp, cov.col.names, s)
  
}

#run sens analysis of pseudo coords ----
if(iter) {
  run_start_date <- gsub("-", "", Sys.Date()) ##for organizing saves
  
  ##set variables ----
  run.location = "/amarel_cluster" # "/amarel_cluster"
  do.bias.correct = TRUE 
  calc_range_in_models = FALSE
  calc_eff_area_in_models = FALSE
  run_model_opt = TRUE
  
  
  ### SELECT SURVEY INDEX ###
  #all, fjs, or bss
  survey2run = "all"
  
  ### SELECT A SPECIES TO RUN ###
  spp = "Blueback Herring"
  spp.wd = paste0(PATH, "/inshore_offshore", run.location, "/", gsub(" ", "_", spp))
  
  if(!dir.exists(spp.wd)) { dir.create(spp.wd) }
  if(!dir.exists(paste0(spp.wd, "/vast_models"))) { dir.create(paste0(spp.wd, "/vast_models")) }
  if(!dir.exists(paste0(spp.wd, "/psuedo_coord_iter"))) { dir.create(paste0(spp.wd, "/psuedo_coord_iter")) }
  
  ##load in prepped catch dat ----
  io.dat.spp <- fread(paste0(spp.wd, "/", gsub(" ", "_", spp), "_pc_sens_vast_prep_dat_", survey2run, ".csv"))
  setnames(io.dat.spp, c("lat", "long"), c("Lat", "Lon"))
  
  ##fjs coordinates only started in 2000, so only use years starting there (regardless of which survey set) ----
  io.dat.spp <- io.dat.spp[year >= 2000]
  
  ##load in extrap grid ----
  ext.g <- fread(paste0(PATH, "/inshore_offshore/hudson_river_VAST_extrap_grid_250cs.csv"))
  ext.g <- as.data.frame(ext.g)
  ext.g$Depth <- ext.g$depth ##necessary for strata limits
  
  #copy reduced extrap grid if it exists into species folder (creating it takes a long time)
  if(file.exists(paste0(PATH, "/inshore_offshore/Kmeans_extrapolation-2000.RData")) &
     !file.exists(paste0(spp.wd, "/Kmeans_extrapolation-2000.RData"))) {
    
    file.copy(from = paste0(PATH, "/inshore_offshore/Kmeans_extrapolation-2000.RData"),
              to = paste0(spp.wd, "/Kmeans_extrapolation-2000.RData"))
  }
  
  ##pull in settings from selected model ----
  #directory of 'top' selected models
  mod.top <- fread(paste0(PATH, "/inshore_offshore", run.location, "/multi-species_cov_comp_results_filtered.csv"))
  mod.top <- mod.top[species == spp & survey_set == survey2run]
  
  top.v.m <- readRDS(paste0(spp.wd, "/vast_models/", mod.top$model_name))
  
  settings = top.v.m$settings
  
  settings$n_x <- 200 #reset knots to 200, otherwise some will take FOREVER to run through it all
  settings$ObsModel <- c(4, 0) #optim does not like 5 0 in this format for some reason
  
  ###covariates ----
  catchability_data <- names(top.v.m$catchability_data)
  catchability_data <- io.dat.spp[, .SD, .SDcols = catchability_data]
  
  if(survey2run != "bss") {
    catchability_data$gear_def <- relevel(as.factor(catchability_data$gear_def), ref = "1m Tucker Trawl")
  } else {
    catchability_data$gear_def <- relevel(as.factor(catchability_data$gear_def), ref = "100' x 10' Beach Seine")
  }
  
  covariate_data <- names(top.v.m$covariate_data)
  
  if(!is.null(covariate_data)) {
    covariate_data <- io.dat.spp[, .SD, .SDcols = covariate_data[covariate_data != "Year"]]
    covariate_data[, Year := NA]
  }
  
  Q1_formula = top.v.m$Q1_formula
  Q2_formula = top.v.m$Q2_formula
  # X1_formula = top.v.m$X1_formula
  # X2_formula = top.v.m$X2_formula
  X1_formula <- ~ bs(log(riv_dpth + 1), df = 3)
  X2_formula <- ~ bs(log(riv_dpth + 1), df = 3)
  
  ##can turn these off to save time on runs (doesn't have to be the same from earlier runs)
  settings$Options[["Calculate_Range"]] <- calc_range_in_models
  settings$Options[["Calculate_effective_area"]] <- calc_eff_area_in_models
  
  ## vast model run w/ missing coord dat removed: ----
  model.save.name <- paste0(gsub(" ", "_", spp), "_", survey2run, "_pscs_coords_removed_only_", run_start_date, ".rds")
  
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
    saveRDS(fit, paste0(spp.wd, "/psuedo_coord_iter/", model.save.name))
  }
  
  tmp1 <- data.table(species = spp, 
                     sens_type = "missing coords removed, no pseudocoords at all",
                     save_name = model.save.name,
                     error_message = error.message)
  
  fwrite(tmp1, paste0(spp.wd, "/psuedo_coord_iter/model_directory_", survey2run, "_pscs_coords_removed_only_", run_start_date, ".csv"))
  
  ## vast model run w/ missing coord dat removed, iterate in fake coordinates: ----
  #get proportion of data that was pseudocoords
  o.counts <- fread(paste0(PATH, "/inshore_offshore/species_datasets_final_counts_for_model_input.csv"))
  psc.counts <- fread(paste0(PATH, "/inshore_offshore/species_datasets_final_counts_for_model_input_psc_sens.csv"))
  ## num. removed rows (missing coords) / total tows or rows for that survey
  fjs.prop <- psc.counts[species == spp & survey == "fjs", rows_removed_for_missing_coords] / o.counts[species == spp & survey == "fjs", total_tows]
  bss.prop <- psc.counts[species == spp & survey == "bss", rows_removed_for_missing_coords] / o.counts[species == spp & survey == "bss", total_tows]
  
  ##each iteration should have how many pseudocoords?
  fjs.psc.cts <- round(nrow(io.dat.spp[survey == "fjs"]) * fjs.prop)
  bss.psc.cts <- round(nrow(io.dat.spp[survey == "bss"]) * bss.prop)
  
  tmp <- data.table(ind = 1:nrow(io.dat.spp))
  
  for(i in 1:10) { #repeat randomization 10 times
    
    model.save.name <- paste0(gsub(" ", "_", spp), "_", survey2run, "_pscs_randomized_iter", i, "_", run_start_date, ".rds")
    
    ###set fake coords ----
    #first, randomly pick tows that will have the pscs, based on the previous calculated proportion
    io.dat.spp$assign_psc <- numeric(nrow(io.dat.spp))
    
    if(survey2run != "all") {
      
      if(survey2run == "fjs") { 
        io.dat.spp[sample(nrow(io.dat.spp), fjs.psc.cts), assign_psc := 1] 
      } else {
        if(survey2run == "bss") { 
          io.dat.spp[sample(nrow(io.dat.spp), bss.psc.cts), assign_psc := 1] 
        } else {
          stop("Something went wrong, did you specify the survey incorrectly?")
        }
      }
      
      io.dat.spp[assign_psc == 1, `:=` (Lat = lat_psc, Lon = long_psc, benth_cat = benth_cat_psc)]
      
    } else {
      
      tmp.io.dat <- split(io.dat.spp, io.dat.spp$survey)
      
      tmp.io.dat[["fjs"]][sample(nrow(tmp.io.dat[["fjs"]]), fjs.psc.cts), assign_psc := 1]
      tmp.io.dat[["bss"]][sample(nrow(tmp.io.dat[["bss"]]), bss.psc.cts), assign_psc := 1]
      
      io.dat.spp <- rbindlist(tmp.io.dat)
      
      io.dat.spp[assign_psc == 1, `:=` (Lat = lat_psc, Lon = long_psc, benth_cat = benth_cat_psc)]
      
      ##have to remake covariate dfs because it's not in the same order anymore
      catchability_data <- names(top.v.m$catchability_data)
      catchability_data <- io.dat.spp[, .SD, .SDcols = catchability_data]
      
      if(survey2run != "bss") {
        catchability_data$gear_def <- relevel(as.factor(catchability_data$gear_def), ref = "1m Tucker Trawl")
      } else {
        catchability_data$gear_def <- relevel(as.factor(catchability_data$gear_def), ref = "100' x 10' Beach Seine")
      }
      
      covariate_data <- names(top.v.m$covariate_data)
      covariate_data <- io.dat.spp[, .SD, .SDcols = covariate_data[covariate_data != "Year"]]
      covariate_data[, Year := NA]
      
    }
    
    #run model
    fit <- tryCatch(
      fit_model( settings = settings,
                 Lat_i = io.dat.spp$Lat,
                 Lon_i = io.dat.spp$Lon,
                 t_i = io.dat.spp$year,
                 b_i = as_units(io.dat.spp$ct_yoy, "count"),
                 a_i = as_units(io.dat.spp$area_swept, "m^2"),
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
      saveRDS(fit, paste0(spp.wd, "/psuedo_coord_iter/", model.save.name))
    }
    
    tmp1 <- data.table(species = spp,
                       sens_type = "random psc",
                       iter = i,
                       num_psc_fjs = nrow(io.dat.spp[survey == "fjs" & assign_psc == 1]),
                       num_psc_bss = nrow(io.dat.spp[survey == "bss" & assign_psc == 1]),
                       save_name = model.save.name,
                       error_message = error.message)
    
    fwrite(tmp1, paste0(spp.wd, "/psuedo_coord_iter/model_directory_", survey2run, "_pscs_randomized_iter", i, "_", run_start_date, ".csv"))
    
  }
}


#compare results of psc runs ----
if(compare_results) {
  
  model.run.date = "20251031|20260126"
  
  #combine 'directories'
  mod.dir <- lapply(target.species, function(s) {
    
    spp.wd <- paste0(PATH, "/inshore_offshore", model.run.location, "/", gsub(" ", "_", s), "/psuedo_coord_iter")
    
    f <- list.files(spp.wd, pattern = "model_directory", full.names = T)
    f <- grep(model.run.date, f, value = T)
    
    f <- lapply(f, fread)
    
    f <- rbindlist(f, use.names = T, fill = T)
    
    f[, survey_set := fifelse(grepl("all", save_name), "inshore+offshore",
                              fifelse(grepl("fjs", save_name), "offshore",
                                      "inshore"))]
    
    return(f)
  })
  
  mod.dir <- rbindlist(mod.dir)
  
  mod.top <- fread(paste0(PATH, "/inshore_offshore", model.run.location, "/multi-species_cov_comp_results_filtered.csv"))
  
  #for every model that exists, pull abundance estimate and se
  source("scripts/inshore_offshore/00_functions.R")
  
  ai.io.spp.l <- lapply(target.species, function(s) {
    
    spp.wd <- paste0(PATH, "/inshore_offshore", model.run.location, "/", gsub(" ", "_", s))
    
    mod.name <- mod.dir[species == s & is.na(error_message), save_name]
    mod.name <- c(mod.name, mod.top[species == s, model_name])
    
    ai.io <- vector("list", length(mod.name))
    for(i in 1:length(mod.name)) {
      
      if(!file.exists(file.path(spp.wd, "psuedo_coord_iter", mod.name[i])) &
         !file.exists(file.path(spp.wd, "vast_models", mod.name[i]))) { next }
      
      if(grepl("pscs", mod.name[i])) {
        ai.io[[i]] <- get_index(model.save.location = paste0(spp.wd, "/psuedo_coord_iter"), model.save.name = mod.name[i])
      } else {
        ai.io[[i]] <- get_index(model.save.location = paste0(spp.wd, "/vast_models"), model.save.name = mod.name[i])
      }
      
      ai.io[[i]][, `:=` (cv_yr = std_error_for_estimate / estimate,
                         species = s)]

      for(surv in c("fjs", "bss", "all")) {
        if(grepl(surv, mod.name[i])) { ai.io[[i]]$survey_set <- surv }
      }
      
      if(grepl("randomized", mod.name[i])) {
        ai.io[[i]]$run_type <- stringr::str_extract(mod.name[i], "randomized_iter[0-9]+")
      } else {
        if(grepl("coords_removed", mod.name[i])) {
          ai.io[[i]]$run_type <- "pscs coords removed"
        } else {
          ai.io[[i]] <- ai.io[[i]][time > 1999]
          ai.io[[i]]$run_type <- "original"
          ##recalc normalized estimate based only on years with all coord available
          ai.io[[i]][, `:=` (norm_estimate = estimate / mean(estimate),
                             norm_std_error = std_error_for_estimate / mean(estimate))]
        }
      }
      
    }
    
    ai.io <- rbindlist(ai.io)
    
    return(ai.io)
    
  })
  
  ai.io.spp.l <- rbindlist(ai.io.spp.l)
  
  ai.io.spp.l[survey_set == "all", survey_set := "inshore+offshore"]
  ai.io.spp.l[survey_set == "bss", survey_set := "inshore"]
  ai.io.spp.l[survey_set == "fjs", survey_set := "offshore"]
  
  io_theme <- theme(panel.border = element_rect(color = "black", fill = NA),
                    panel.background = element_rect(fill = "white"),
                    panel.grid.major = element_line(color = "lightgray"),
                    strip.background = element_rect(color = "black", fill = "gray"),
                    axis.text.x = element_text(angle = 47, vjust = 1, hjust = 1),
                    axis.text = element_text(size = 12, color = "black"),
                    strip.text = element_text(size = 12, face = "bold"),
                    axis.title = element_text(size = 12),
                    legend.text = element_text(size = 12, color = "black"))
  
  iter.colors <- viridis::viridis(10)
  names(iter.colors) <- paste0("randomized_iter", seq(1, 10, 1))
  iter.colors["pscs coords removed"] <- "blue"
  iter.colors["original"] <- "red"
  ai.io.spp.l[, run_type := factor(run_type, levels = c("original", "pscs coords removed", 
                                                        paste0("randomized_iter", seq(1, 10, 1))))]
  
  ggplot() +
    geom_path(data = ai.io.spp.l, aes(x = time, y = norm_estimate, color = run_type), linewidth = 1) +
    # geom_point(data = ai.io.spp.l, aes(x = time, y = norm_estimate, color = run_type), size = 2) +
    scale_color_manual(values = iter.colors) +
    facet_wrap(species ~ survey_set, scales = "free_y", ncol = 3) +
    scale_y_continuous(name = "normalized abundance index") +
    scale_x_continuous(n.breaks = 5, name = "year") +
    io_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.key.size = unit(1, 'cm'),
          legend.title = element_text(hjust = 0, size = 12),
          legend.text = element_text(hjust = 0, margin = margin(r = 7)))
  ggsave(paste0(PATH, "/inshore_offshore", model.run.location, "/figures/pscs_sensitivity_abundance_indices.png"), width = 11, height = 7)
  
  ai.io.spp.l[, random := gsub("_iter[0-9]+", "", run_type)]
  ggplot(data = ai.io.spp.l, aes(x = random, y = cv_yr)) +
    geom_boxplot(outliers = FALSE) +
    geom_point(shape = 21, position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
               show.legend = FALSE, alpha = 0.6, size = 2) +
    facet_wrap(species ~ survey_set, scales = "free_y", ncol = 3) +
    scale_y_continuous(transform = "log10", n.breaks = 5, 
                       name = bquote(annual~CV~(SE~"/"~abundance~estimate))) +
    io_theme +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1), legend.position = "right") +
    labs(x = "")
  ggsave(paste0(PATH, "/inshore_offshore", model.run.location, "/figures/pscs_sensitivity_cv.png"), width = 10, height = 11)
  
  #test for differences
  test <- ai.io.spp.l[run_type != "original"]
  wt.res <- data.table()
  for(spp in target.species) {
    for(surv in unique(test$survey_set)) {
      tmp <- test[species == spp & survey_set == surv]
      cv_rem <- tmp[run_type == "pscs coords removed", cv_yr]
      cv_rand <- tmp[random == "randomized", cv_yr]
      wt <- wilcox.test(cv_rem, cv_rand)
      
      df <- data.table(species = spp, survey_set = surv, statistic = wt$statistic, p_val = round(wt$p.value, 4),
                       removed_larger_cv = median(cv_rem) > median(cv_rand)) 
      wt.res <- rbindlist(list(wt.res, df))
      
    }
  }
  
  fwrite(wt.res, file.path(PATH, "inshore_offshore", model.run.location, "pseudo_coord_sens_wilcox_res.csv"))
  
  #get sensitivity
  ai.cv <- ai.io.spp.l[, 
                       .(cv = sd(norm_estimate) / mean(norm_estimate)),
                       by = c("species", "survey_set", "time")]
  fwrite(ai.cv, file.path(PATH, "inshore_offshore", model.run.location, "pseudo_coord_sens_ai_cv.csv"))
  
  ai.env <- ai.io.spp.l[, {
    orig <- norm_estimate[run_type == "original"]
    rem_val <- norm_estimate[run_type == "pscs coords removed"]
    rand_vals <- norm_estimate[grepl("random", run_type)]
    
    .(
      orig = orig,
      rem_val = rem_val,
      q10  = quantile(rand_vals, 0.05),
      q90  = quantile(rand_vals, 0.95)
    )
  }, by = c("species", "survey_set", "time")]
  fwrite(ai.env, file.path(PATH, "inshore_offshore", model.run.location, "pseudo_coord_sens_ai_env.csv"))
  
  cv.cv <- ai.io.spp.l[, 
                       .(cv = sd(cv_yr) / mean(cv_yr)),
                       by = c("species", "survey_set", "time")]
  fwrite(cv.cv, file.path(PATH, "inshore_offshore", model.run.location, "pseudo_coord_sens_cv_cv.csv"))
  
  cv.env <- ai.io.spp.l[, {
    orig <- cv_yr[run_type == "original"]
    rem_val <- cv_yr[run_type == "pscs coords removed"]
    rand_vals <- cv_yr[grepl("random", run_type)]
    
    .(
      orig = orig,
      rem_val = rem_val,
      q10  = quantile(rand_vals, 0.05),
      q90  = quantile(rand_vals, 0.95)
    )
  }, by = c("species", "survey_set", "time")]
  fwrite(cv.env, file.path(PATH, "inshore_offshore", model.run.location, "pseudo_coord_sens_cv_env.csv"))
  
  
  ##summarize
  ai.cv <- fread(file.path(PATH, "inshore_offshore", model.run.location, "pseudo_coord_sens_ai_cv.csv"))
  summ.ai <- ai.cv[, .(mn_cv_ai = round(mean(cv), 3)), by = c("species", "survey_set")]
  
  cv.cv <- fread(file.path(PATH, "inshore_offshore", model.run.location, "pseudo_coord_sens_cv_cv.csv"))
  summ.cv <- cv.cv[, .(mn_cv_cv = round(mean(cv), 3)), by = c("species", "survey_set")]
  cv.cv[, .(mn_cv = round(mean(cv), 3))]
  
  cv.env <- fread(file.path(PATH, "inshore_offshore", model.run.location, "pseudo_coord_sens_cv_env.csv"))
  cv.env[, env_location := fifelse(orig <= q90 & orig >= q10, "inside", 
                                   fifelse(orig > q90, "up", "down"))]
  table(cv.env[, .(survey_set, env_location, species)])
  
  ##supplementary table
  summ.psc <- summ.ai[summ.cv, on = c("species", "survey_set")]
  fwrite(summ.psc, file.path(PATH, "inshore_offshore", model.run.location, "pseudo_coord_sens_summary_res.csv"))
  
}

