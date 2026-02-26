### COMPARE INDICES AND MODEL OUTPUTS ###
#this script includes:
##comparison of vast derived abundance indices
##comparison of predicted density and observed density estimates

library(data.table); library(VAST)

PATH <- getwd()

source(paste0(PATH, "/scripts/inshore_offshore/XX_colors.R"))

##### SET VARS #####
spp.l = c("Alewife", "American Shad", "Blueback Herring", "Striped Bass")
model.run.location = "amarel_cluster/" #""
model.run.type = "knot_comp"

#-###############-#

#directory of 'top' selected models
mod.top <- fread(paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_", gsub(" ", "_", model.run.type), "_results_filtered.csv"))

# compare survey cpue ----
load_spp_dat <- function(spp) {
  
  spp.wd <- paste0(PATH, "/inshore_offshore/", gsub(" ", "_", spp))
  
  sp.dat <- list.files(spp.wd, pattern = "catch_dat", full.names = T)
  
  sp.dat <- lapply(sp.dat, function(DT) {
    
    DT <- fread(DT)
    if(all(DT$survey == "fjs")) { DT$surv_label <- "fjs" }
    if(all(DT$survey == "bss")) { DT$surv_label <- "bss" }
    if(length(unique(DT$survey)) == 2) { DT$surv_label <- "all" }
    
    return(DT)
    
  })
  
  sp.dat <- rbindlist(sp.dat)
  return(sp.dat)
}

se <- function(x) sd(x)/sqrt(length(x))

spp.summ.all <- vector("list", length = length(spp.l))
catch.summ <- vector("list", length = length(spp.l))
sp.dat.all <- vector("list", length = length(spp.l))

for(i in seq_along(spp.l)) {
  
  sp.dat <- load_spp_dat(spp = spp.l[[i]])
  sp.dat <- sp.dat[surv_label != "all"]
  
  #get cpue (catch / m2)
  # sp.dat[, km2 := area_swept / 1e6]
  sp.dat[, cpue := ct_yoy / area_swept]
  
  sp.dat.all[[i]] <- sp.dat
  
  #get survey averages
  t1 <- sp.dat[ct_yoy > 0, .(mn_pos_catch_per_haul = round(mean(ct_yoy), 3), 
                             mn_pos_catch_cpue = round(mean(cpue), 3),
                             min_cpue = round(min(cpue), 3)), 
               by = c("common_name", "survey")]
  t2 <- sp.dat[, .(prop_zero_catch = round(mean(ct_yoy == 0), 3),
                   mn_cpue = round(mean(cpue), 3), 
                   min_cpue = round(min(cpue), 3),
                   max_cpue = round(max(cpue), 3)), by = c("common_name", "survey")]
  
  catch.summ[[i]] <- t1[t2, on = c("common_name", "survey")]
  
  #summarize for time series
  spp.summ.all[[i]] <-  sp.dat[, .(mn_cpue = mean(cpue),
                                   sd_cpue = sd(cpue),
                                   se_cpue = se(cpue)), 
                               by = c("common_name", "survey", "year")]
}

spp.summ.all <- rbindlist(spp.summ.all)

fwrite(spp.summ.all, paste0(PATH, "/inshore_offshore/multi-species_summarized_cpue.csv")) 

catch.summ <- rbindlist(catch.summ)

fwrite(catch.summ, paste0(PATH, "/inshore_offshore/multi-species_mean_cpue_per_survey.csv")) 

sp.dat.all <- rbindlist(sp.dat.all)

##survey summary
surv.summ <- sp.dat.all[, .(survey, year, uniq_id, km2)] |> unique()
surv.summ <- surv.summ[,  .(ttl_tows = .N, ttl_km = sum(km2)), by = c("survey", "year")]
surv.summ[, .(ttl_tows = sum(ttl_tows), mn_tows = mean(ttl_tows), 
              min_tows = min(ttl_tows), max_tows = max(ttl_tows)), c("survey")]

ggplot(data = catch.summ) +
  geom_point(aes(x = common_name, y = prop_zero_catch, color = survey), size = 5)


##correlation of annual cpue fjs v. bss ----
sprmn.corr <- lapply(spp.l, function(x) {
  cor.dt <- spp.summ.all[common_name == x]
  cor.dt <- dcast(cor.dt, common_name + year ~ survey, value.var = "mn_cpue")
  
  cor.res <- cor.test(cor.dt$bss, cor.dt$fjs, method = "spearman")
  
  dt <- data.table(species = x,
                   rho = cor.res$estimate,
                   p_val = cor.res$p.value)
  
  return(dt)
  })

sprmn.corr <- rbindlist(sprmn.corr)
sprmn.corr[, p_val := round(p_val, digits = 3)]

fwrite(sprmn.corr, paste0(PATH, "/inshore_offshore/multi-species_spearman_corr_mean_cpue.csv")) 

# compare predicted v. observed density ----
get_pred_obs_dens <- function(model.save.location, model.save.name, vast.mod = NULL) {
  
  if(is.null(vast.mod)) {
    #read in model
    vast.mod <- readRDS(paste(model.save.location, model.save.name, sep = "/"))
    # vast.mod <- reload_model(vast.mod)
  }
  
  out <- data.table(
    lat = vast.mod$data_frame$Lat_i,
    long = vast.mod$data_frame$Lon_i,
    time = vast.mod$data_frame$t_i,
    gear = vast.mod$catchability_data$gear_def,
    pred_dens = vast.mod$Report$D_i,
    pred_pos_cat = vast.mod$Report$R2_i,
    pred_encounter = vast.mod$Report$R1_i,
    obs_catch = strip_units(vast.mod$data_frame$b_i),
    area_swept = strip_units(vast.mod$data_frame$a_i)
  )
  
  out[, orig_surv := fifelse(grepl("Beach Seine", gear), "bss", "fjs")]
  
  return(out)
}

po.spp.l <- lapply(spp.l, function(spp) {
  
  spp.wd <- paste0(PATH, "/inshore_offshore/", model.run.location, gsub(" ", "_", spp))
  
  mod.name <- mod.top[species == spp, model_name]
  
  po.io <- vector("list", length(mod.name))
  for(i in 1:length(mod.name)) {
    po.io[[i]] <- get_pred_obs_dens(model.save.location = paste0(spp.wd, "/vast_models"), model.save.name = mod.name[i])
    
    for(surv in c("fjs", "bss", "all")) {
      if(grepl(surv, mod.name[i])) { po.io[[i]]$survey_set <- surv }
    }
  }
  
  po.io <- rbindlist(po.io)
  
  po.io$species <- spp
  
  return(po.io)
  
})

###### get pred v. obs relationship ----
##for annual summary
lm.res.all <- vector("list", length(po.spp.l))
for(i in seq_along(po.spp.l)) {
  
  lm.res <- lapply(c("all", "fjs", "bss"), function(s) {
    tmp <- po.spp.l[[i]][, .(obs_catch, area_swept, pred_dens, time, species, survey_set)][survey_set == s]
    tmp[, `:=` (obs_dens = obs_catch / (area_swept))] 
    tmp <- tmp[, lapply(.SD, mean, na.rm = TRUE), by = .(time, species, survey_set), .SDcols = c("obs_catch", "obs_dens", "pred_dens")]
    
    lm.out <- lm(log(pred_dens) ~ log(obs_dens), tmp)
    
    summary_fit <- summary(lm.out)
    r2 <- summary_fit$r.squared
    p_value <- summary_fit$coefficients[2, 4]
    
    ci <- confint(lm.out)
    
    lm.res <- data.table(species = unique(tmp$species),
                         survey_set = s,
                         intercept = coef(lm.out)["(Intercept)"],
                         intercept_lowerci = ci["(Intercept)", "2.5 %"],
                         intercept_upperci = ci["(Intercept)", "97.5 %"],
                         slope = coef(lm.out)["log(obs_dens)"],
                         slope_lowerci = ci["log(obs_dens)", "2.5 %"],
                         slope_upperci = ci["log(obs_dens)", "97.5 %"],
                         r2 = summary(lm.out)$r.squared,
                         p_val = summary(lm.out)$coefficients[2,4]
    )
    lm.res[, names(which(sapply(lm.res, is.numeric))) := lapply(.SD, round, digits = 3), .SDcols = is.numeric]
    
    return(lm.res)
  })

  lm.res <- rbindlist(lm.res)
  
  lm.res.all[[i]] <- lm.res
  
}

lm.res.all <- rbindlist(lm.res.all)

fwrite(lm.res.all, paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_obsvpred_lm_results_", model.run.type, ".csv"))

po.spp.l <- rbindlist(po.spp.l)

fwrite(po.spp.l, paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_density_estimates_", model.run.type, ".csv"))


# compare abundance indices ----
## for each species, read in the top model for each index, pull out index values
source("scripts/inshore_offshore/00_functions.R")

ai.io.spp.l <- lapply(spp.l, function(spp) {
  
  spp.wd <- paste0(PATH, "/inshore_offshore/", model.run.location, gsub(" ", "_", spp))
  
  mod.name <- mod.top[species == spp, model_name]
  
  ai.io <- vector("list", length(mod.name))
  for(i in 1:length(mod.name)) {
    ai.io[[i]] <- get_index(model.save.location = paste0(spp.wd, "/vast_models"), model.save.name = mod.name[i])
    
    for(surv in c("fjs", "bss", "all")) {
      if(grepl(surv, mod.name[i])) { ai.io[[i]]$survey_set <- surv }
    }
  }
  
  ai.io <- rbindlist(ai.io)
  
  ai.io$species <- spp
  
  return(ai.io)
  
})

ai.io.spp.l <- rbindlist(ai.io.spp.l)

fwrite(ai.io.spp.l, paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_abundance_estimates_", model.run.type, ".csv"))

# comparing spatial density estimates -------
get_knot_density <- function(model.save.location, model.save.name, vast.mod = NULL) {
  
  if(is.null(vast.mod)) {
    #read in model
    vast.mod <- readRDS(paste(model.save.location, model.save.name, sep = "/"))
    # vast.mod <- reload_model(vast.mod)
  }
  
  knots.ll <- as.data.table(vast.mod$spatial_list$latlon_g) #or latlong_x (its the same, not sure if this is specific to our use case though)
  
  d_gct <- as.data.table(vast.mod$Report$D_gct) |> drop_units()

  d_gct$Site <- as.character(seq(1, nrow(d_gct), 1))
  d_gct$long <- knots.ll$Lon
  d_gct$lat <- knots.ll$Lat
  
  d_gct <- melt(d_gct, measure.vars = grep("\\.", names(d_gct), value = T), variable.name = "year", value.name = "pred_dens")
  d_gct[, `:=` (year = as.integer(gsub("1\\.", "", year)))]
  
  r1_gct <- as.data.table(vast.mod$Report$R1_gct)[, .(Site, year = as.integer(Time), r1_gct = value)]
  r2_gct <- as.data.table(vast.mod$Report$R2_gct)[, .(Site, year = as.integer(Time), r2_gct = value)]
  
  d_gct <- d_gct[r1_gct, on = c("Site", "year")]
  d_gct <- d_gct[r2_gct, on = c("Site", "year")]
  
  d_gct[, Site := NULL]
  
  return(d_gct)

  }


d.knot.all <- lapply(spp.l, function(spp) {
  
  spp.wd <- paste0(PATH, "/inshore_offshore/", model.run.location, gsub(" ", "_", spp))
  
  mod.name <- mod.top[species == spp, model_name]
  
  d.gct.io <- vector("list", length(mod.name))
  for(i in 1:length(mod.name)) {
    d.gct.io[[i]] <- get_knot_density(model.save.location = paste0(spp.wd, "/vast_models"), model.save.name = mod.name[i])
    
    for(surv in c("fjs", "bss", "all")) {
      if(grepl(surv, mod.name[i])) { d.gct.io[[i]]$survey_set <- surv }
    }
  }
  
  d.gct.io <- rbindlist(d.gct.io)
  
  d.gct.io$species <- spp
  
  return(d.gct.io)
  
})

d.knot.all <- rbindlist(d.knot.all)

fwrite(d.knot.all, paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_knot_spatial_density_estimates_", model.run.type, ".csv"))

###### get 'residuals' of vast types ----
d.knot.wide <- dcast(data = d.knot.all, species + year + lat + long ~ survey_set, value.var = "pred_dens")
d.knot.wide[, paste0("norm_", c("all", 
                                "bss", 
                                "fjs")) := lapply(.SD, function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))), 
            by = .(species, year), .SDcols = c("all", 
                                               "bss", 
                                               "fjs")]

d.knot.wide[, `:=` (abs_diff_all_bss = all - bss,
                    abs_diff_all_fjs = all - fjs,
                    abs_diff_bss_fjs = bss - fjs,
                    norm_diff_all_bss = norm_all - norm_bss,
                    norm_diff_all_fjs = norm_all - norm_fjs,
                    norm_diff_bss_fjs = norm_bss - norm_fjs
                    )]

fwrite(d.knot.wide, paste0(PATH,  "/inshore_offshore/", model.run.location, "multi-species_knot_spatial_density_residuals_", model.run.type, ".csv"))

# compare index uncertainty ----
######  comparing CV and CIs ----
ai.spp <- fread(paste0(PATH,  "/inshore_offshore/", model.run.location, "multi-species_abundance_estimates_", model.run.type, ".csv"))
setnames(ai.spp, c("time", "std_error_for_estimate"), c("year", "std_error"))
ai.spp[, type := "VAST"]

cpue <- spp.summ.all[, .(species = common_name, year, estimate = mn_cpue, std_error = se_cpue, 
                         type = "CPUE", survey_set = survey,
                         upper = mn_cpue + 1.96 * se_cpue, lower = mn_cpue - 1.96 * se_cpue)]
cpue[, mn_estimate := mean(estimate), by = c("species", "survey_set")]
cpue[, `:=` (norm_estimate = estimate / mn_estimate,
             norm_std_error = std_error / mn_estimate)]
cpue[, `:=` (norm_upper = norm_estimate + 1.96 * norm_std_error,
             norm_lower = norm_estimate - 1.96 * norm_std_error)]

ai.dec <- fread(paste0(PATH, "/inshore_offshore/multi-species_abundance_estimates_vDEC.csv")) ##no cluster version
ai.dec <- ai.dec[, .(species, year, estimate = gm, std_error = jack_se, upper = jack_uci, lower = jack_lci, 
                     type = "DEC geomean", survey_set = "bss")]
ai.dec <- ai.dec[year %in% range(ai.spp$year)[1]:range(ai.spp$year)[2]]
ai.dec[, mn_estimate := mean(estimate, na.rm = TRUE), by = c("species", "survey_set")]
ai.dec[, `:=` (norm_estimate = estimate / mn_estimate,
               norm_std_error = std_error / mn_estimate)]
ai.dec[, `:=` (norm_upper = norm_estimate + 1.96 * norm_std_error,
               norm_lower = norm_estimate - 1.96 * norm_std_error)]


all.estimates <- rbindlist(list(ai.spp, cpue, ai.dec), use.names = T, fill = TRUE)

all.estimates[, `:=` (cv_yr = std_error / estimate,
                      norm_cv_yr = norm_std_error / norm_estimate, #this will be basically equivalent to cv_yr ... btw
                      rel_ciw_yr = (upper - lower) / estimate)]

fwrite(all.estimates, paste0(PATH,  "/inshore_offshore/", model.run.location, "multi-species_simple_uncertainty_estimates_", model.run.type, ".csv"))

###### test for normality, variance homogeneity ----
cv.est <- fread(paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_simple_uncertainty_estimates_", mod.type, ".csv"))

norm.lt.dt <- data.table()
for(s in spp.l) {
  x <- cv.est[species == s, cv_yr]
  x.log <- log(cv.est[species == s, cv_yr] + 1)
  
  #normality
  s.x <- shapiro.test(x)
  s.x.log <- shapiro.test(x.log)
  
  #variance homogeneity
  lt <- car::leveneTest(x ~ factor(cv.est[species == s, survey_set]))
  lt.log <- car::leveneTest(x.log ~ factor(cv.est[species == s, survey_set]))
  
  tmp <- data.table(species = s, 
                    SW_stat = s.x$statistic,
                    SW_pval = s.x$p.value,
                    SW_log_stat = s.x.log$statistic,
                    SW_log_pval = s.x.log$p.value,
                    LT_fstat = lt$`F value`[1],
                    LT_pval = lt$`Pr(>F)`[1],
                    LT_log_fstat = lt$`F value`[1],
                    LT_log_pval = lt$`Pr(>F)`[1])
  
  tmp <- tmp[, .(species, round(.SD, digits = 3)), .SDcols = is.numeric]
  
  norm.lt.dt <- rbindlist(list(norm.lt.dt, tmp)) 
}

fwrite(norm.lt.dt, file = paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_CV_test_results_shapiro_levene_", mod.type, ".csv"))

###### compare group cvs ----
kt.l <- lapply(spp.l, function(s) {
  set.seed(234)
  
  ## create matrix with only selected variables (ie doesn't make sense to do raw concentrations with ratios)
  x <- cv.est[species == s & type == "VAST", .(species, cv_yr, survey_set)]
  
  kt <- kruskal.test(x$cv_yr ~ x$survey_set)
  pw.wt <- pairwise.wilcox.test(x$cv_yr, x$survey_set, p.adjust.method = "BH")
  pw.wt <- as.data.table(pw.wt$p.value, keep.rownames = T)
  
  tmp <- data.table(species = s,
                    kt.stat = kt$statistic,
                    p.val = kt$p.value,
                    pw_allvbss = pw.wt[rn == "bss", all],
                    pw_allvfjs = pw.wt[rn == "fjs", all],
                    pw_bssvfjs = pw.wt[rn == "fjs", bss])
  tmp <- tmp[, .(species, round(.SD, digits = 3)), .SDcols = is.numeric]
  
  return(tmp)
  }
)

kt.l <- rbindlist(kt.l)
fwrite(kt.l, file = paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_CV_test_results_kruskal_", mod.type, ".csv"))

#comparison of covariate coefficients ----
library(effects); library(splines)

cov.fig.save.path <- file.path(fig.save.path, "covariate_plots")

for(i in 1:nrow(mod.top)) {
  ##set vars
  spp = mod.top[i, species]
  spp.wd <- paste0(PATH, "/inshore_offshore/", model.run.location, gsub(" ", "_", spp))
  model.save.location = paste0(spp.wd, "/vast_models")
  survey = mod.top[i, survey_set]
  mod.name = mod.top[i, model_name]
  
  ##load model and assoc. data
  vast.mod <- readRDS(file.path(model.save.location, mod.name))
  
  ### must add data-frames to global environment (from wiki)
  covariate_data_full = vast.mod$effects$covariate_data_full
  catchability_data_full = vast.mod$effects$catchability_data_full
  
  ##fix benthic cat factor error
  if(grepl("factor\\(benth_cat\\)", as.character(vast.mod$X1_formula[2]))) {
    
    vast.mod$X1_formula <- update(vast.mod$X1_formula,  ~ . - factor(benth_cat) + benth_cat)
    vast.mod$X2_formula <- update(vast.mod$X2_formula,  ~ . - factor(benth_cat) + benth_cat)
    
  }
  
  ##plots
  for(pred.type in c("X1", "X2", "Q1", "Q2")) {
    
    if(grepl("1", pred.type)) { y.axis = "effect on encounter probability" } 
    if(grepl("2", pred.type)) { y.axis = "effect on positive density" } 
    
    if(grepl("X", pred.type)) {
      covs <- names(covariate_data_full)
      covs <- covs[covs %in% c("riv_dpth", "benth_cat")]
    }
    if(grepl("Q", pred.type)) {
      covs <- names(catchability_data_full)
      covs <- covs[covs %in% c("gear", "solar_noon_diff", "sam_dpth")]
    }

    if(is.null(covs)) { next }
    
    for(d in covs) {
      
      pred <- Effect.fit_model(vast.mod,
                               focal.predictors = d,
                               which_formula = pred.type,
                               xlevels = 100,
                               transformation = list(link=identity, inverse=identity))
      
      png(paste0(cov.fig.save.path, "/", gsub(" ", "_", spp), "_", survey, "_", d, "_", pred.type, ".png"), 
          width = 800, height = 800)
      p <- plot(pred, xlab = d, 
                main = paste0(spp, ", ", survey, ", ", pred.type, ", ", d),
                ylab = y.axis)
      print(p)
      dev.off()
      
    }
    
    pred <- Effect.fit_model(vast.mod,
                             focal.predictors = covs,
                             which_formula = pred.type,
                             xlevels = 100,
                             transformation = list(link=identity, inverse=identity))
    
    png(paste0(cov.fig.save.path, "/", gsub(" ", "_", spp), "_", survey, "_interaction_", pred.type, ".png"), 
        width = 800, height = 800)
    p <- plot(pred, main = paste0(spp, ", ", survey, ", ", pred.type, ", interaction"), ylab = y.axis)
    print(p)
    dev.off()
    
  }
  
  cat(round(i/nrow(mod.top)*100, 2), "%")
}

#comparison metrics from Cacciopaglia et al. 2024 ----
###### ratio ----
ratio <- dcast(all.estimates[!type %in% c("CPUE")], formula = species + year ~ type + survey_set, value.var = "estimate")

#ratio to gm
ratio[, paste0("ratio_", gsub("VAST_", "", grep("VAST", names(ratio), value = TRUE)), "_gm") := 
        lapply(.SD, function(x) x / `DEC geomean_bss`), 
      .SDcols = grep("VAST", names(ratio), value = TRUE)]

#ratio to combined ind.
ratio[, paste0("ratio_", c("bss", "fjs"), "_combined") := 
        lapply(.SD, function(x) x / VAST_all), 
      .SDcols = c("VAST_bss", "VAST_fjs")]

#summarize 
ratio.summ <- ratio[, lapply(.SD, mean, na.rm = T), .SDcols = grep("ratio", names(ratio), value = TRUE), by = c("species")]

fwrite(ratio, paste0(PATH,  "/inshore_offshore/", model.run.location, "multi-species_index_ratios_raw_", model.run.type, ".csv"))
fwrite(ratio.summ, paste0(PATH,  "/inshore_offshore/", model.run.location, "multi-species_index_ratios_summ_", model.run.type, ".csv"))

##comparing old and new results ----
if(FALSE) {
  library(ggplot2)
  s <- "all"
  spp <- "American Shad"
  ind <- get_index(vast.mod = fit)
  ind[, `:=` (cv_yr = std_error_for_estimate / estimate,
              norm_cv_yr = norm_std_error / norm_estimate, #this will be basically equivalent to cv_yr ... btw
              rel_ciw_yr = (upper - lower) / estimate)]
  
  orig_res <- fread(paste0(PATH, "/inshore_offshore/multi-species_simple_uncertainty_estimates.csv"))
  tmp <- orig_res[survey_set == s & species == spp & type == "VAST"]
  
  mean(ind$cv_yr)
  mean(tmp$cv_yr)
  
  p1 <- ggplot() +
    geom_path(data = tmp, aes(x = year, y = norm_estimate)) +
    geom_path(data = ind, aes(x = time, y = norm_estimate, color = strata)) +
    ggtitle("blue is new vast model results", subtitle = "knots = 200,\ndecreased res. of extrapolation grid")
  
  p2 <- ggplot() +
    geom_path(data = tmp, aes(x = year, y = estimate)) +
    geom_path(data = ind, aes(x = time, y = estimate, color = strata))
  
  p3 <- ggplot() +
    geom_point(data = tmp, aes(x = year, y = cv_yr)) +
    geom_point(data = ind, aes(x = time, y = cv_yr, color = strata)) +
    ggtitle(paste0("new mean cv = ", round(mean(ind$cv_yr), 3), ";\nold mean cv = ", round(mean(tmp$cv_yr), 3)))
  
  ggpubr::ggarrange(p1, p2, p3, nrow = 1, align = "v", common.legend = T)
  ggsave(paste0(PATH, "/figures/inshore_offshore_index/example_updated_vast_model_results_", gsub(" ", "_", spp), "_", s, "_strata.png"), width = 12, height = 4)
  
}

