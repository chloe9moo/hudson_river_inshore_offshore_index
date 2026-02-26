### DEC index for comparison to vast indices ###
library(data.table); library(bootstrap)

PATH <- getwd()

#load in dec bss data (see 00_bss_data_prep.R)
bss.alosa <- fread(paste0(PATH, "/bio_data/dec_bss_alosa_yoy_counts.csv"))
bss.sb <- fread(paste0(PATH, "/bio_data/dec_bss_striped_bass_yoy_counts.csv"))

#make list of each species w/ common count col name
spp.v <- c("Shad", "Ale", "BB")
yoy.cts <- vector("list", length = 4)
for(i in 1:length(spp.v)) {
  spp <- spp.v[[i]]
  yoy.cts[[i]] <- bss.alosa[, .SD, .SDcols = c("Batch", "YEAR", spp)]
  setnames(yoy.cts[[i]], spp, "yoy")
  yoy.cts[[i]]$species <- spp
  yoy.cts[[i]][, species := fifelse(species == "Shad", "American Shad",
                                    fifelse(species == "Ale", "Alewife",
                                            fifelse(species == "BB", "Blueback Herring",
                                                    species)))]
}

yoy.cts[[4]] <- bss.sb[, .(Batch, YEAR, yoy = SB, species = "Striped Bass")]

#set functions ----
##traditional gm method ----
traditiona_geomean_calc <- function(DT, id_col, count_col, time_col){
  x <- copy(DT)
  
  x[, `:=` (yoy = get(count_col),
            uniq_id = get(id_col))]
  x[, yoy1 := yoy + 1]
  
  gmci <- x[, .(hauls = length(yoy1), 
                catch = sum(yoy), 
                am = round(mean(yoy), 3),
                am.sd = sd(yoy), 
                ml = mean(log(yoy1)), 
                sdl = sd(log(yoy1))), by = time_col]
  
  traditional <- gmci[, .(hauls, catch, am, 
                          am.se = round(am.sd/sqrt(hauls), 3),
                          uci = round(exp((ml + (1.96*sdl)/sqrt(hauls)))-1, 3), #upper confidence int
                          lci = round(exp((ml - (1.96*sdl)/sqrt(hauls)))-1, 3), #lower conf int
                          gm = round(exp(ml) - 1, 3)), #geomean
                      by = time_col]
  
  #get zero hauls
  zeros <- subset(x, yoy == 0)
  zeros <- zeros[, .(zero.hauls = length(uniq_id)), by = time_col]
  # print(zeros)
  
  traditional <- zeros[traditional, on = time_col]
  traditional[, zero.hauls := fifelse(is.na(zero.hauls), 0, zero.hauls)]
  
  return(traditional)
}

##testing...
# traditiona_geomean_calc(bss.alosa, id_col = "Batch", count_col = "BB", time_col = "YEAR")

## jackknife method -----
## from Gary apparently..
jackknife_se_calc <- function(DT,  id_col, count_col, time_col) {
  
  back <- function(x) { exp(mean(x))-1 }
  jk_log <- function(x) { 
    out <- jackknife(log(x), back)
    out <- data.table(jack_se = round(out$jack.se, 3),
                      jack_bias = round(out$jack.bias, 3))
    return(out)
    }
  
  s <- copy(DT)
  
  s[, yoy1 := get(count_col) + 1]
  
  s.split <- split(s, s[, get(time_col)])
  
  jk.out <- lapply(s.split, function(x) {
    jk_log(x$yoy1)
  })
  
  time.v <- names(jk.out)
  jk.out <- rbindlist(jk.out)
  jk.out$time <- as.numeric(time.v)
  setnames(jk.out, "time", time_col)
  
  return(jk.out)
}

#testing..
# jackknife_se_calc(bss.alosa, id_col = "Batch", count_col = "BB", time_col = "YEAR")

## bootstrap method ----
##also from Gary..

bootstrap_se_calc <- function(DT, nboot = 500, id_col, count_col, time_col) {
  
  geomean_est <- function(x) { exp(mean(sample(log(x), length(x), replace=TRUE))) - 1 }
  
  s <- copy(DT)
  
  s[, yoy1 := get(count_col) + 1]
  
  bs <- s[, .(t = (replicate(nboot, geomean_est(yoy1)))), by = time_col]
  
  bs.se <- bs[, .(boot_se = round(sd(t), 3)), by = time_col]
  
  return(bs.se)
  
}

#testing...
# bootstrap_se_calc(bss.alosa, id_col = "Batch", count_col = "BB", time_col = "YEAR")

#full dec gm ai run through ----
## (from original code): Note that I am using the inf t value of 1.96 rather than using the qt function in 
## the CI calcs... other than that the CI calcs are basically equation 10 in the write-up

get_dec_ai <- function(DT, nboot = 500, id_col, count_col, time_col) {
  
  s <- copy(DT)
  
  s.gm <- traditiona_geomean_calc(s, id_col = id_col, count_col = count_col, time_col = time_col)
  jk.out <- jackknife_se_calc(s, id_col = id_col, count_col = count_col, time_col = time_col)
  bs.out <- bootstrap_se_calc(s, nboot = nboot, id_col = id_col, count_col = count_col, time_col = time_col)
  
  ai <- merge(s.gm, jk.out, by = time_col)
  ai <- merge(ai, bs.out, by = time_col)
  
  ai[, `:=` (jack_uci = gm + 1.96 * jack_se,
             jack_lci = gm - 1.96 * jack_se,
             boot_uci = gm + 1.96 * boot_se,
             boot_lci = gm - 1.96 * boot_se,
             species = unique(s$species))]
  
  setcolorder(ai, c("species", time_col))
  
  return(ai)
  
}

yoy.ai <- lapply(yoy.cts, get_dec_ai, id_col = "Batch", count_col = "yoy", time_col = "YEAR")
yoy.ai <- rbindlist(yoy.ai)

setnames(yoy.ai, c("YEAR"), c("year"))

#save
fwrite(yoy.ai, paste0(PATH, "/inshore_offshore/multi-species_abundance_estimates_vDEC.csv"))

# asa fjs design-based index ----
##averaged density (number of individuals divided by the volume of water sampled) over all surveyed regions, strata, and weeks



# #plot check
# library(ggplot2)
# 
# ggplot(data = yoy.ai, aes(x = year, y = gm)) +
#   geom_path() +
#   geom_point() +
#   facet_wrap(~ species, scales = "free_y")


