### parameter comparisons ###
# compare estimates and cv of all models successfully run

library(data.table); library(VAST); library(ggplot2)
source("scripts/inshore_offshore/00_functions.R")
source("scripts/inshore_offshore/XX_colors.R")

PATH <- getwd()

spp.l = c("Alewife", "Striped Bass", "American Shad", "Blueback Herring")

model.run.location = "amarel_cluster/" #""

param.theme <- list(
  theme_bw(),
  theme(axis.text = element_text(size = 12, color = "black"), 
        axis.title.x = element_blank(),
        legend.text = element_text(size = 12, color = "black"),
        legend.position = "bottom")
)

#obs. model comparison ----
obs.mod <- fread(paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_omsens_results_20250918.csv"))
top.obs.mod <- fread(paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_omsens_results_filtered_20250918.csv"))

ai.all <- data.table()
for(i in seq_len(nrow(obs.mod))) {
  
  spp.wd <- paste0(PATH, "/inshore_offshore/", model.run.location, gsub(" ", "_", obs.mod[i, species]))
  
  if(file.exists(paste0(spp.wd, "/vast_models/", obs.mod[i, save_name])) & 
     obs.mod[i, convergence_check] == "There is no evidence that the model is not converged") {
    
    vast.mod <- readRDS(paste0(spp.wd, "/vast_models/", obs.mod[i, save_name]))
    
    ai <- get_index(vast.mod = vast.mod)
    
    ai[, `:=` (species = obs.mod[i, species], obs_mod = obs.mod[i, obs_mod])]
    
    if(obs.mod[i, save_name] %in% top.obs.mod$save_name) {
      ai[, selected_model := TRUE]
    } else {
      ai[, selected_model := FALSE]
    }
    
    for(surv in c("fjs", "bss", "all")) {
      if(grepl(surv, obs.mod[i, save_name])) { ai[, survey_type := surv] }
    }
    
    ai.all <- rbindlist(list(ai.all, ai))
    
  }
  
  print_progress(iter_num = i, total_iter = nrow(obs.mod))
  
}

fwrite(ai.all, paste0(PATH, "/inshore_offshore/", model.run.location, "param_comp_obs_model_abundance_indices.csv"))


## summarize ----
ai.all[, max_est := max(estimate), by = c("species", "survey_type")]
ai.summ <- ai.all[, .(sd = sd(estimate), mn = mean(estimate), cv = sd(estimate) / mean(estimate), max_est), by = c("species", "survey_type", "time")] |> unique()
ai.summ <- ai.summ[, .(mn_sd = mean(sd), mn_cv = mean(cv), max_est), by = .(species, survey_type)] |> unique()
ai.summ[, c("mn_sd", "mn_cv") := lapply(.SD, function(x) round(x, 2)), .SDcols = c("mn_sd", "mn_cv")]

## for plotting
ai.summ[, `:=` (x = 2010, mn_cv_label = paste0("mean CV = ", mn_cv))]

## plot abundance ----
ggplot(data = ai.all, aes(x = time, y = estimate)) +
  geom_path(aes(color = obs_mod, linetype = selected_model)) +
  geom_text(data = ai.summ, aes(x = x, y = max_est, label = mn_cv)) +
  facet_wrap(species ~ survey_type, scales = "free_y", nrow = 4) +
  scale_color_viridis_d(end = 0.8, name = "ObsModel setting") +
  scale_linetype_discrete(name = "selected model?") +
  scale_x_continuous(n.breaks = 3) +
  labs(x = "", y = "abundance index") +
  param.theme
ggsave(paste0(PATH, "/inshore_offshore/", model.run.location, "/figures/obs_model_abundance_index_comparisons.png"), width = 8, height = 8, bg = "white")

##sensitivity of obs model choice for estimate ----
obs.mod.corr <- data.table()
for(spp in spp.l) {
  for(surv in  c("fjs", "bss", "all")) {
    
    tmp <- ai.all[species == spp & survey_type == surv, .(time, estimate, obs_mod)]
    tmp <- dcast(tmp, time ~ obs_mod, value.var = "estimate")
    
    for(c.type in c("pearson", "spearman")) {
      
      tmp.cor <- Hmisc::rcorr(as.matrix(tmp[, -"time"]), type = c.type)
      
      tmp.cor.res <- as.data.table(tmp.cor$r, keep.rownames = TRUE)
      tmp.cor.p <- as.data.table(tmp.cor$P)
      setnames(tmp.cor.p, names(tmp.cor.p), paste0(gsub(" ", "_", names(tmp.cor.p)), "_pval"))
      
      tmp.all <- cbind(tmp.cor.res, tmp.cor.p)
      
      tmp.all[, `:=` (species = spp, survey_type = surv, corr_type = c.type)]
      
      obs.mod.corr <- rbindlist(list(obs.mod.corr, tmp.all), use.names = TRUE, fill = TRUE)
    }
  }
}

fwrite(obs.mod.corr, paste0(PATH, "/inshore_offshore/", model.run.location, "param_comp_obs_model_abundance_indices_correlation_results.csv"))

## plot cv, annual and average ----
ai.all[, `:=` (cv_yr = std_error_for_estimate / estimate,
               rel_ciw_yr = (upper - lower) / estimate)]
ai.all[, max_cv := max(cv_yr), by = .(species, survey_type)]

cv.summ <- ai.all[, .(sd = sd(cv_yr), mn = mean(cv_yr), cv = sd(cv_yr) / mean(cv_yr), max_cv), by = c("species", "survey_type", "time")] |> unique()
cv.summ <- cv.summ[, .(mn_sd = mean(sd), mn_cv = mean(cv), max_cv), by = .(species, survey_type)] |> unique()
cv.summ[, c("mn_sd", "mn_cv") := lapply(.SD, function(x) round(x, 2)), .SDcols = c("mn_sd", "mn_cv")]
cv.summ[, x := 2010]

ggplot(data = ai.all, aes(x = time, y = cv_yr)) +
  geom_point(aes(fill = obs_mod, shape = selected_model)) +
  geom_text(data = cv.summ, aes(x = x, y = max_cv, label = mn_cv)) +
  facet_wrap(species ~ survey_type, scales = "free_y", nrow = 4) +
  scale_fill_viridis_d(end = 0.8, name = "ObsModel setting") +
  scale_shape_manual(values = c("TRUE" = 24, "FALSE" = 21), name = "selected model?") +
  scale_x_continuous(n.breaks = 3) +
  labs(x = "", y = "CV (SE / abundance estimate)") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)),
         shape = guide_legend(override.aes = list(size = 3))) +
  param.theme +
  theme(legend.box = "vertical")
ggsave(paste0(PATH, "/inshore_offshore/", model.run.location, "/figures/obs_model_annual_CV_comparisons.png"), width = 8, height = 8, bg = "white")

cv.summ <- ai.all[, .(mn_cv = mean(cv_yr)), by = c("species", "survey_type", "obs_mod")]

ggplot(data = cv.summ, aes(x = survey_type, y = mn_cv, shape = obs_mod, color = species)) +
  geom_point(size = 5) +
  scale_color_manual(values = spp.pal) +
  facet_wrap(~ species, nrow = 1) +
  theme_classic()
ggsave(paste0(PATH, "/inshore_offshore/", model.run.location, "/figures/obs_model_mean_CV_comparisons.png"), width = 8, height = 4, bg = "white")


## aic, rmse
# for(i in seq_len(nrow(obs.mod))) {
#   if(obs.mod[i, save_name] %in% top.obs.mod$save_name) {
#     obs.mod[i, selected_model := TRUE]
#   } else {
#     obs.mod[i, selected_model := FALSE]
#   }
# }

# ggplot(data = obs.mod, aes(x = aic, y = rmse)) +
#   geom_point(aes(fill = obs_mod, shape = selected_model)) +
#   facet_wrap(species ~ survey, scales = "free", nrow = 4) +
#   scale_fill_viridis_d(end = 0.8, name = "ObsModel setting") +
#   scale_shape_manual(values = c("TRUE" = 24, "FALSE" = 21), name = "selected model?") +
#   param.theme

###not really sure this is a useful plot....

# cov comp comparison ----
cov.mod <- fread(paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_cov_comp_results.csv"))
cov.mod[, obsmod4 := FALSE]
cov.mod4 <- fread(paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_cov_comp_4_results.csv"))
cov.mod4[, obsmod4 := TRUE]
cov.mod <- rbindlist(list(cov.mod, cov.mod4))

top.cov.mod <- fread(paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_cov_comp_results_filtered.csv"))
top.cov.mod[, obsmod4 := FALSE]
top.cov.mod4 <- fread(paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_cov_comp_4_results_filtered.csv"))
top.cov.mod4[, obsmod4 := TRUE]
top.cov.mod <- rbindlist(list(top.cov.mod, top.cov.mod4))

ai.all <- data.table()
for(i in seq_len(nrow(cov.mod))) {
  
  spp.wd <- paste0(PATH, "/inshore_offshore/", model.run.location, gsub(" ", "_", cov.mod[i, species]))
  
  if(file.exists(paste0(spp.wd, "/vast_models/", cov.mod[i, model_name])) & 
     cov.mod[i, convergence_check] == "There is no evidence that the model is not converged") {
    
    vast.mod <- readRDS(paste0(spp.wd, "/vast_models/", cov.mod[i, model_name]))
    
    ai <- get_index(vast.mod = vast.mod)
    
    ai[, `:=` (species = cov.mod[i, species], 
               obs_mod = cov.mod[i, obs_mod],
               survey_type = cov.mod[i, survey_set],
               catch_cov = cov.mod[i, catch_cov],
               dens_cov = cov.mod[i, dens_cov],
               obsmod4_run = cov.mod[i, obsmod4])]
    
    if(cov.mod[i, model_name] %in% top.cov.mod$model_name) {
      ai[, selected_model := TRUE]
    } else {
      ai[, selected_model := FALSE]
    }
    
    ai.all <- rbindlist(list(ai.all, ai))
    
  }
  
  print_progress(iter_num = i, total_iter = nrow(cov.mod))
  
}

fwrite(ai.all, paste0(PATH, "/inshore_offshore/", model.run.location, "param_comp_cov_abundance_indices.csv"))

##plot abundance ----
ai.all[, catch_dens_covs := paste0(catch_cov, "\n", dens_cov)]
ai.all[, catch_dens_covs := gsub("bs\\(log\\(|\\+ 1\\),| df = 3\\)|ns\\(|, df = 5\\)|factor\\(|\\)", "", catch_dens_covs)]
ai.all[, catch_dens_covs := gsub("  ", " ", catch_dens_covs)]
ai.all[, catch_dens_covs := gsub("dpth $", "dpth", catch_dens_covs)]
unique(ai.all$catch_dens_covs)

ai.all[, select_obs := paste0(selected_model, ", ", obsmod4_run)]

ggplot(data = ai.all, aes(x = time, y = estimate)) +
  geom_path(aes(color = catch_dens_covs, linetype = select_obs)) +
  # geom_text(data = ai.summ, aes(x = x, y = max_est, label = mn_cv)) +
  facet_wrap(species ~ survey_type, scales = "free_y", nrow = 4) +
  # scale_color_viridis_d(end = 0.8, name = "covariates") +
  scale_linetype_manual(name = "selected model? obs 4 run?", 
                        values = c("FALSE, FALSE" = "dotted", "FALSE, TRUE" = "dashed", "TRUE, FALSE" = "solid", "TRUE, TRUE" = "longdash")) +
  scale_y_continuous(transform = "log1p") +
  scale_x_continuous(n.breaks = 3) +
  labs(x = "", y = "abundance index") +
  param.theme +
  theme(legend.box = "vertical")
ggsave(paste0(PATH, "/inshore_offshore/", model.run.location, "/figures/covariate_abundance_index_comparisons.png"), width = 16, height = 12, bg = "white")

##sensitivity of obs model choice for estimate ----
ai.all[, comp_col := paste0(catch_dens_covs, " ", select_obs)]
cov.corr <- data.table()
for(spp in spp.l) {
  for(surv in  c("fjs", "bss", "all")) {
    
    tmp <- ai.all[species == spp & survey_type == surv, .(time, estimate, comp_col)]
    tmp <- dcast(tmp, time ~ comp_col, value.var = "estimate")
    
    for(c.type in c("pearson", "spearman")) {
      
      tmp.cor <- Hmisc::rcorr(as.matrix(tmp[, -"time"]), type = c.type)
      
      tmp.cor.res <- as.data.table(tmp.cor$r, keep.rownames = TRUE)
      tmp.cor.p <- as.data.table(tmp.cor$P)
      setnames(tmp.cor.p, names(tmp.cor.p), paste0(gsub(" ", "_", names(tmp.cor.p)), "_pval"))
      
      tmp.all <- cbind(tmp.cor.res, tmp.cor.p)
      
      tmp.all[, `:=` (species = spp, survey_type = surv, corr_type = c.type)]
      
      cov.corr <- rbindlist(list(cov.corr, tmp.all), use.names = TRUE, fill = TRUE)
    }
  }
}

fwrite(cov.corr, paste0(PATH, "/inshore_offshore/", model.run.location, "param_comp_cov_abundance_indices_correlation_results.csv"))

View(cov.corr[grepl("TRUE, FALSE|TRUE, TRUE", rn), .SD, 
              .SDcols = c("rn", "species", "survey_type", "corr_type", 
                          grep("TRUE, FALSE|TRUE, TRUE", names(cov.corr), value = T))])

## plot cv, annual and average ----
ai.all[, `:=` (cv_yr = std_error_for_estimate / estimate,
               rel_ciw_yr = (upper - lower) / estimate,
               select_obs = factor(select_obs, levels = c("TRUE, TRUE", "TRUE, FALSE",
                                                          "FALSE, TRUE", "FALSE, FALSE")))]
ai.all[, max_cv := max(cv_yr), by = .(species, survey_type)]

cv.summ <- ai.all[, .(sd = sd(cv_yr), mn = mean(cv_yr), cv = sd(cv_yr) / mean(cv_yr), max_cv), by = c("species", "survey_type", "time")] |> unique()
cv.summ <- cv.summ[, .(mn_sd = mean(sd), mn_cv = mean(cv), max_cv), by = .(species, survey_type)] |> unique()
cv.summ[, c("mn_sd", "mn_cv") := lapply(.SD, function(x) round(x, 2)), .SDcols = c("mn_sd", "mn_cv")]
cv.summ[, x := 2010]

ggplot(data = ai.all, aes(x = time, y = cv_yr)) +
  geom_point(aes(fill = catch_dens_covs, shape = select_obs, color = select_obs)) +
  # geom_text(data = cv.summ, aes(x = x, y = max_cv, label = mn_cv)) +
  facet_wrap(species ~ survey_type, scales = "free_y", nrow = 4) +
  scale_fill_discrete(name = "covariates") +
  scale_shape_manual(name = "selected model? obs model 4 run?",
                     values = c("FALSE, FALSE" = 23, "FALSE, TRUE" = 24, "TRUE, FALSE" = 21, "TRUE, TRUE" = 22)) +
  scale_color_manual(name = "selected model? obs model 4 run?",
                     values = c("FALSE, FALSE" = "black", "FALSE, TRUE" = "black", "TRUE, FALSE" = "red", "TRUE, TRUE" = "red")) +
  scale_x_continuous(n.breaks = 3) +
  scale_y_continuous(transform = "log10") +
  labs(x = "", y = "CV (SE / abundance estimate)") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)),
         shape = guide_legend(override.aes = list(size = 3))) +
  param.theme +
  theme(legend.box = "vertical")
ggsave(paste0(PATH, "/inshore_offshore/", model.run.location, "/figures/cov_comp_annual_CV_comparisons.png"), width = 16, height = 12, bg = "white")

cv.summ <- ai.all[, .(mn_cv = mean(cv_yr)), by = c("species", "survey_type", "catch_dens_covs", "select_obs")]

ggplot(data = cv.summ, aes(x = survey_type, y = mn_cv, fill = catch_dens_covs, 
                           shape = select_obs, color = select_obs, alpha = select_obs)) +
  geom_point(size = 5) +
  scale_fill_discrete(name = "covariates") +
  scale_shape_manual(name = "selected model? obs model 4 run?",
                     values = c("FALSE, FALSE" = 23, "FALSE, TRUE" = 24, "TRUE, FALSE" = 21, "TRUE, TRUE" = 22)) +
  scale_color_manual(name = "selected model? obs model 4 run?",
                     values = c("FALSE, FALSE" = "black", "FALSE, TRUE" = "black", "TRUE, FALSE" = "red", "TRUE, TRUE" = "red")) +
  scale_alpha_manual(name = "selected model? obs model 4 run?",
                     values = c("FALSE, FALSE" = 0.5, "FALSE, TRUE" = 0.5, "TRUE, FALSE" = 1, "TRUE, TRUE" = 1)) +
  facet_wrap(~ species, nrow = 1) +
  scale_y_continuous(transform = "log10") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)),
         shape = guide_legend(override.aes = list(size = 3))) +
  theme_classic()
ggsave(paste0(PATH, "/inshore_offshore/", model.run.location, "/figures/cov_comp_mean_CV_comparisons.png"), width = 9, height = 8, bg = "white")


#compare knot ----
which_obs_run = "" #or "_4"
knot.mod <- fread(paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_knot_comp", which_obs_run, "_results.csv"))
top.knot.mod <- fread(paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_knot_comp", which_obs_run, "_results_filtered.csv"))

ai.all <- data.table()
for(i in seq_len(nrow(knot.mod))) {
  
  spp.wd <- paste0(PATH, "/inshore_offshore/", model.run.location, gsub(" ", "_", knot.mod[i, species]))
  
  if(file.exists(paste0(spp.wd, "/vast_models/", knot.mod[i, model_name])) & 
     knot.mod[i, convergence_check] == "There is no evidence that the model is not converged") {
    
    vast.mod <- readRDS(paste0(spp.wd, "/vast_models/", knot.mod[i, model_name]))
    
    ai <- get_index(vast.mod = vast.mod)
    
    ai[, `:=` (species = knot.mod[i, species], 
               obs_mod = knot.mod[i, obs_mod],
               survey_type = knot.mod[i, survey_set],
               catch_cov = knot.mod[i, catch_cov],
               dens_cov = knot.mod[i, dens_cov],
               knots = knot.mod[i, knots])]
    
    if(knot.mod[i, model_name] %in% top.knot.mod$model_name) {
      ai[, selected_model := TRUE]
    } else {
      ai[, selected_model := FALSE]
    }
    
    ai.all <- rbindlist(list(ai.all, ai))
    
  }
  
  print_progress(iter_num = i, total_iter = nrow(knot.mod))
  
}

fwrite(ai.all, paste0(PATH, "/inshore_offshore/", model.run.location, "param_comp_knot_abundance_indices.csv"))

##plot abundance ----
ai.all[, knots := factor(knots, levels = unique(ai.all$knots)[order(unique(ai.all$knots))])]
ggplot(data = ai.all, aes(x = time, y = estimate)) +
  geom_path(aes(color = factor(knots), linetype = selected_model)) +
  # geom_text(data = ai.summ, aes(x = x, y = max_est, label = mn_cv)) +
  facet_wrap(species ~ survey_type, scales = "free_y", nrow = 4) +
  scale_color_viridis_d(end = 0.8, name = "knots") +
  scale_linetype_manual(name = "selected model?", values = c("TRUE" = "solid", "FALSE" = "dashed")) +
  # scale_y_continuous(transform = "log1p") +
  scale_x_continuous(n.breaks = 3) +
  labs(x = "", y = "abundance index") +
  param.theme +
  theme(legend.box = "vertical")
ggsave(paste0(PATH, "/inshore_offshore/", model.run.location, "/figures/knot_abundance_index_comparisons", which_obs_run, ".png"), width = 10, height = 10, bg = "white")

##sensitivity of obs model choice for estimate ----
knot.corr <- data.table()
for(spp in spp.l) {
  for(surv in  c("fjs", "bss", "all")) {
    
    tmp <- ai.all[species == spp & survey_type == surv, .(time, estimate, knots)]
    tmp <- dcast(tmp, time ~ knots, value.var = "estimate")
    
    for(c.type in c("pearson", "spearman")) {
      
      tmp.cor <- Hmisc::rcorr(as.matrix(tmp[, -"time"]), type = c.type)
      
      tmp.cor.res <- as.data.table(tmp.cor$r, keep.rownames = TRUE)
      tmp.cor.p <- as.data.table(tmp.cor$P)
      setnames(tmp.cor.p, names(tmp.cor.p), paste0(gsub(" ", "_", names(tmp.cor.p)), "_pval"))
      
      tmp.all <- cbind(tmp.cor.res, tmp.cor.p)
      
      tmp.all[, `:=` (species = spp, survey_type = surv, corr_type = c.type)]
      
      knot.corr <- rbindlist(list(knot.corr, tmp.all), use.names = TRUE, fill = TRUE)
    }
  }
}

fwrite(knot.corr, paste0(PATH, "/inshore_offshore/", model.run.location, "param_comp_knot_abundance_indices_correlation_results.csv"))

## plot cv, annual and average ----
ai.all[, `:=` (cv_yr = std_error_for_estimate / estimate,
               rel_ciw_yr = (upper - lower) / estimate)]
ai.all[, max_cv := max(cv_yr), by = .(species, survey_type)]

cv.summ <- ai.all[, .(sd = sd(cv_yr), mn = mean(cv_yr), cv = sd(cv_yr) / mean(cv_yr), max_cv), by = c("species", "survey_type", "time")] |> unique()
cv.summ <- cv.summ[, .(mn_sd = mean(sd), mn_cv = mean(cv), max_cv), by = .(species, survey_type)] |> unique()
cv.summ[, c("mn_sd", "mn_cv") := lapply(.SD, function(x) round(x, 2)), .SDcols = c("mn_sd", "mn_cv")]
cv.summ[, x := 2010]

ggplot(data = ai.all, aes(x = time, y = cv_yr)) +
  geom_point(aes(fill = knots, shape = selected_model, color = selected_model)) +
  # geom_text(data = cv.summ, aes(x = x, y = max_cv, label = mn_cv)) +
  facet_wrap(species ~ survey_type, scales = "free_y", nrow = 4) +
  scale_fill_viridis_d(name = "knots") +
  scale_shape_manual(name = "selected model?",
                     values = c("TRUE" = 21, "FALSE" = 22)) +
  scale_color_manual(name = "selected model?",
                     values = c("FALSE" = "black","TRUE" = "red")) +
  scale_x_continuous(n.breaks = 3) +
  scale_y_continuous(transform = "log10") +
  labs(x = "", y = "CV (SE / abundance estimate)") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)),
         shape = guide_legend(override.aes = list(size = 3))) +
  param.theme +
  theme(legend.box = "vertical")
ggsave(paste0(PATH, "/inshore_offshore/", model.run.location, "/figures/knot_annual_CV_comparisons.png"), width = 10, height = 10, bg = "white")

cv.summ <- ai.all[, .(mn_cv = mean(cv_yr)), by = c("species", "survey_type", "knots", "selected_model")]

ggplot(data = cv.summ, aes(x = survey_type, y = mn_cv, fill = knots, shape = selected_model, color = selected_model)) +
  geom_point(size = 5) +
  scale_fill_viridis_d(name = "knots") +
  scale_shape_manual(name = "selected model?",
                     values = c("TRUE" = 21, "FALSE" = 22)) +
  scale_color_manual(name = "selected model?",
                     values = c("FALSE" = "black","TRUE" = "red")) +
  facet_wrap(~ species, nrow = 1) +
  scale_y_continuous(transform = "log10") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)),
         shape = guide_legend(override.aes = list(size = 3))) +
  theme_classic()
ggsave(paste0(PATH, "/inshore_offshore/", model.run.location, "/figures/knot_mean_CV_comparisons.png"), width = 8, height = 5, bg = "white")



#compare trajectory of results as parameter decisions are made ----
