### STANDARDIZE AND COMBINE FJS AND BSS DATA ###
# generally meant to make adjustments to fall juvenile survey and dec beach seine data so they can be joined appropriately 
# also make some comparisons

# set up
library(data.table); library(sf)

PATH <- getwd()
if(grepl("scripts", PATH)) { PATH <- gsub("/scripts", "", PATH) }

# load in cleaned data for each + prep
##fall juvenile survey
fjs <- fread(paste0(PATH, "/bio_data/fall_juvenile_survey_combined.csv"))
#remove records with problematic use codes (sampling problems) and incomplete sample years
##note, use code 4 is not in the data dictionary, but it always occurs when there was no survey/catch data avail.
fjs <- fjs[!use_code %in% c(2, 4, 5) & 
             !year %in% c(2015, 2016) &
             !gear_def == "2m Beam Trawl"] ##only 67 unique tows using this gear and its not in any metadata..?

##beach seine
bss <- fread(paste0(PATH, "/bio_data/beach_seine_plus_combined.csv"))

#keep only striped bass or alosine beach seine data
bss <- bss[prog_code %in% c(4, 7) & 
             #remove flagged performance issues; remove aborted runs (no effort / no catch) as per comments
             flag_perf != TRUE & flag_siteabort != TRUE]

# View(data.frame(survey = c(rep("fjs", length(names(fjs))), rep("bss", length(names(bss)))),
#                 name = c(names(fjs), names(bss))))

# standardize species names ----
## make a species list and save
spp.list <- lapply(list(fjs, bss), function(x) {
  x[, unique(.SD), .SDcols = c("taxon_code", "species_name")]
})

spp.list[[1]] <- spp.list[[1]][, survey := "fjs"]
spp.list[[2]] <- spp.list[[2]][, survey := "bss"]

spp.list <- do.call(rbind, spp.list)

##do backbone search for full name to match
if(!file.exists(paste0(PATH, "/bio_data/species_list_fjs_bss.csv"))) {
  library(rgbif)
  
  taxa_nam <- NULL
  for(i in seq_len(nrow(spp.list))) {
    # for(i in 127:nrow(spp.list)) {
    
    Sys.sleep(0.5) ##looking up names too fast seems to throw errors sometimes
    
    taxa1 <- name_lookup(query = spp.list[i, species_name], status = "ACCEPTED", isExtinct = FALSE)
    tmp <- taxa1 ##for testing
    
    taxa1 <- data.table(taxa1[["data"]])
    
    if(nrow(taxa1) == 0) {   #skip if no records returned
      df <- spp.list[i][, common_name := species_name]
      taxa_nam <- rbindlist(list(taxa_nam, df[, .(common_name, survey, taxon_code)]), fill = TRUE)
      next 
    }
    
    cols <- intersect(names(taxa1), c("kingdom", "phylum", "class", "order", "family", "genus", "species")) #get columns that exist
    taxa1 <- taxa1[, ..cols][,
                             `:=` (common_name = spp.list[i, species_name],
                                   survey = spp.list[i, survey],
                                   taxon_code = spp.list[i, taxon_code])]
    taxa1 <- unique(taxa1)
    
    taxa_nam <- rbindlist(list(taxa_nam, taxa1), fill = TRUE)
  }
  
  #filter and adjust to save
  taxa_nam <- taxa_nam[!is.na(taxon_code) & !taxon_code %in% c(604, 999, 33)]
  taxa_nam[, class := as.character(class)]
  
  #save and manually check for correct names
  fwrite(taxa_nam, paste0(PATH, "/bio_data/species_list_fjs_bss.csv"))
  
}

#read in species list with final names selected (on google drive)
taxa_nam <- fread(paste0(PATH, "/bio_data/species_list_fjs_bss.csv"))

taxa_nam <- taxa_nam[action == "keep"][, action := NULL]

#check for any potentially (unexpected) overlooked species:
# spp.list[!taxon_code %in% taxa_nam$taxon_code]
setnames(taxa_nam, "common_name", "name_orig")
taxa_nam[, taxa_name := fifelse(species != "", species, 
                                fifelse(genus != "", paste0(genus, " spp."),
                                        fifelse(family != "", family, name_orig)))]
taxa_nam[taxa_name == "Unidentifiable", taxa_name := ""] ##fix this to be empty, not as a name

#fix common names so they match across surveys too
taxa_nam[, common_name := name_orig[survey == "fjs"][1], by = taxa_name]
taxa_nam[is.na(common_name), common_name := name_orig]
taxa_nam[species == "Morone saxatilis" & !taxon_code %in% c(579, 888), common_name := "Striped Bass"] #this corrected to hatchery striped bass when it is not
taxa_nam[taxa_name == "Cyprinidae", common_name := "Carp Minnow Family"]

#save
fwrite(taxa_nam, paste0(PATH, "/bio_data/species_list_fjs_bss_clean.csv"))


## add spp to data and save update ----
fjs <- taxa_nam[survey == "fjs"][fjs, on = c("taxon_code", "name_orig" = "species_name")]
fjs <- fjs[, survey := "fjs"] #add survey name to all rows (some will be missing)
fwrite(fjs, paste0(PATH, "/bio_data/fall_juvenile_survey_combined_spp_updated.csv"))

bss <- taxa_nam[survey == "bss"][bss, on = c("taxon_code", "name_orig" = "species_name")]
bss <- bss[, survey := "bss"] #add survey to all rows (some will be missing)

#update empty names
bss[, common_name := fifelse(flag_nofish, "no catch", 
                             fifelse(flag_uncountedfish, "uncounted fish", 
                                     fifelse(taxon_code == 604, "unknown fish", common_name)))]

fwrite(bss, paste0(PATH, "/bio_data/beach_seine_combined_spp_updated.csv"))

# check assigned river mile in data against shapefile ----
##obv can only do this with coordinate data, but still worth doing where possible
rmi <- read_sf(paste0(PATH, "/env_data/gis_layers/NOAA_rivermile_polygon_edited.shp")) |>
  subset(select = -c(PERIMETER, Shape_Leng, POLY_AREA)) |>
  st_transform(4326)

comb.surv <- lapply(list(fjs, bss), function(dt) {
  dt.sf <- dt[!is.na(long)][, unique(.SD), .SDcols = c("riv_mile", "long", "lat", "uniq_id", "survey")]
  #need to do this step bc rmi shapefile doesn't have strata earlier than rmi 12:
  dt.sf[riv_mile < 12, riv_mile := 0]
  dt.sf <- st_as_sf(dt.sf, coords = c("long", "lat"), crs = 4326)
})

comb.surv <- do.call(rbind, comb.surv)

comb.surv <- st_join(comb.surv, rmi)

#find distance of points outside river shapefile to the shapefile
##NOTE: THIS IS CURRENTLY REALLY INEFFICIENT, COULD STREAMLINE AT SOME POINT>>
if(file.exists(paste0(PATH, "/inshore_offshore/sample_dist2_river_mile.csv"))) {
  sample.pt.dist <- fread(paste0(PATH, "/inshore_offshore/sample_dist2_river_mile.csv"))
  # sample.pt.dist[, geometry := NULL]
  
} else {
  sample.pt.dist <- data.table()
  
  for(i in seq_len(nrow(comb.surv))) {
    
    tmp <- comb.surv[i, ]
    
    #do original river mile first
    tmp.riv_mile <- rmi[rmi$Rmile == tmp$riv_mile,]
    
    tmp$dist2riv_mile <- st_distance(tmp, tmp.riv_mile)
    
    #do shapefile assigned river mile:
    
    if(!is.na(tmp$Rmile)) {
      
      if(tmp$riv_mile == tmp$Rmile) { ## if its the same, make distance the same value:
        
        tmp$dist2Rmile <- tmp$dist2riv_mile
        
      } else { ## otherwise find distance to assigned shapefile river mile
        
        tmp.Rmile <- rmi[rmi$Rmile == tmp$Rmile, ]
        tmp$dist2Rmile <- st_distance(tmp, tmp.Rmile)
        
      }
      
    } else { ## if it's outside the shapefile (NA), determine which is closest and how far
      
      dist.rm <- st_distance(tmp, rmi)
      tmp$dist2Rmile <- dist.rm[which.min(dist.rm)]
      tmp$Rmile_new <- rmi[which.min(dist.rm), ]$Rmile #save which one is closest
      
    }
    
    if(!is.na(tmp$riv_mile) & !is.na(tmp$Rmile) & as.numeric(paste0("1", tmp$riv_mile)) == tmp$Rmile) {
      if(tmp$dist2riv_mile > tmp$dist2Rmile) {
        tmp$riv_mile <- tmp$Rmile
        tmp$dist2riv_mile <- tmp$dist2Rmile
      }
      if(tmp$dist2riv_mile < tmp$dist2Rmile) {
        tmp$Rmile <- tmp$riv_mile
        tmp$dist2Rmile <- tmp$dist2riv_mile
      }
    }
    
    sample.pt.dist <- rbindlist(list(sample.pt.dist, tmp), fill = TRUE)
    
    cat(round(i/nrow(comb.surv)*100, 2), "% ...\n")
  }
  
  sample.pt.dist[, geometry := NULL]
  fwrite(sample.pt.dist, paste0(PATH, "/inshore_offshore/sample_dist2_river_mile.csv"))
  
}

comb.surv <- subset(comb.surv, select = -riv_mile)
comb.surv <- sample.pt.dist[as.data.table(comb.surv), on = intersect(names(comb.surv), names(sample.pt.dist))]

##summarize discrepancies in case there are any patterns ----
cs.rivermile.comp <- as.data.table(comb.surv)[, .(ttl_rm2rm = .N), by = .(survey, riv_mile, Rmile, Rmile_new)]
setkeyv(cs.rivermile.comp, c("survey", "riv_mile"))

fwrite(cs.rivermile.comp, paste0(PATH, "/inshore_offshore/combined_survey_rivermile_comparison_ct.csv"))

##flag coordinates that are issues, save all for checking out later in qgis ----
comb.surv[is.na(comb.surv$Rmile) | comb.surv$riv_mile != comb.surv$Rmile, problem_coord := TRUE]
comb.surv[, problem_coord := fifelse(is.na(problem_coord), FALSE, TRUE)]
comb.surv$X <- st_coordinates(comb.surv$geometry)[, "X"]
comb.surv$Y <- st_coordinates(comb.surv$geometry)[, "Y"]

fwrite(comb.surv[, geometry := NULL], paste0(PATH, "/inshore_offshore/combined_survey_rivermile_comparison_pts.csv"))

##FOR NOW: use river mile assignment from river mile shapefile ----
##remove points that are way out of the spatial domain (likely incorrect coordinates)
# comb.surv[, geometry := NULL]
setnames(comb.surv, "riv_mile", "riv_mile_new")

tmp <- lapply(list(fjs, bss), function(x) {
  dt <- copy(x)
  
  tmp <- comb.surv[dt, on = c("uniq_id", "survey")]
  
  #correct join on river mile
  tmp[, riv_mile_new := fifelse(is.na(riv_mile_new), riv_mile, riv_mile_new)]
  tmp[, riv_mile := NULL]
  setnames(tmp, "riv_mile_new", "riv_mile")
  
  #remove all points greater than 250 meters away from any type of river mile (these are way wrong, will need additional info to fix)
  tmp <- tmp[(dist2Rmile <= 100 & dist2riv_mile <= 100) | is.na(dist2Rmile)]
  
  #replace river mile with the one it falls within (if possible)
  ##note: this basically just keeps the original river mile (with the removal of very distant points above), this is more done as a check
  tmp[, riv_mile := fifelse(is.na(lat), riv_mile, #if no coordinates, have to use original riv mile
                            fifelse(!is.na(Rmile) & !is.na(riv_mile) & abs(Rmile - riv_mile) == 1, riv_mile, #if only 1 river mile apart, keep original
                                    fifelse(!is.na(Rmile) & !is.na(riv_mile) & Rmile == riv_mile, riv_mile, #if the same, keep original
                                            fifelse(is.na(Rmile) & (Rmile_new == riv_mile | abs(Rmile_new - riv_mile) == 1), riv_mile,
                                                    9999))))]
  tmp[, c("Rmile", "dist2riv_mile", "dist2Rmile", "Rmile_new", "problem_coord", "X", "Y") := NULL]
  
  #fix na common name (for adding zero catch later)
  tmp[, common_name := fifelse(is.na(common_name), name_orig, common_name)]
  
  return(tmp)
})

# fjs <- tmp[[1]]
# bss <- tmp[[2]]

# combine surveys ----
## select or make common columns, make corrections ----
# intersect(names(fjs), names(bss))
#fjs
tmp[[1]] <- tmp[[1]][, .(uniq_id, riv_mile, survey, common_name, year, ct_yoy, date, time, gear_def, lat, long,
                         strata_new, volume, riv_dpth, sam_dpth, month, jul_day, jul_week, site)]

##fix fjs zero volumes
# mn_vol <- tmp[[1]][volume != 0, .(mn_vol = mean(volume, na.rm = TRUE)), by = c("gear_def", "strata_new", "year")]
# zero_vol <- tmp[[1]][volume == 0]
# 
# zero_vol <- mn_vol[zero_vol, on = intersect(names(mn_vol), names(zero_vol))]
# zero_vol[, volume := mn_vol][,
#                              mn_vol := NULL]
# 
# # row_check <- nrow(tmp[[1]])
# # test <- rbindlist(list(tmp[[1]][volume != 0], zero_vol), use.names = TRUE)
# # row_check == nrow(test)
# 
# tmp[[1]] <- rbindlist(list(tmp[[1]][volume != 0], zero_vol), use.names = TRUE)

##need to add cols to bss
tmp[[2]] <- tmp[[2]][, .(uniq_id, riv_mile, survey, common_name, year, ct_yoy, date, time, gear_def, lat, long,
                         month, jul_day, jul_week, shore_code)]

tmp[[2]][, `:=` (strata_new = "beach",
                 site = shore_code,
                 riv_dpth = 1,
                 sam_dpth = 1, 
                 ##volume estimated from average for each gear in 2017 (see bss_area_swept.R) asssuming height is 1 m
                 volume  = fifelse(grepl("200", gear_def), 488.6766, 470.3452))]
tmp[[2]][, shore_code := NULL]

## join ----
comb.surv <- rbindlist(tmp, use.names = TRUE)
comb.surv[, uniq_id := as.character(uniq_id)]
#need to do this step bc rmi shapefile doesn't have strata earlier than rmi 12:
comb.surv[riv_mile < 12, riv_mile := 0]
fwrite(comb.surv, paste0(PATH, "/inshore_offshore/combined_surveys_fjs_bss.csv"))

comb.surv <- fread(paste0(PATH, "/inshore_offshore/combined_surveys_fjs_bss.csv"),
                   colClasses = c("uniq_id" = "character"))

## get river mile range of each species (could be useful later) ----
tmp <- comb.surv[, unique(.SD), .SDcols = c("common_name", "riv_mile", "survey")]

tmp[, survey := fifelse(survey == "fjs", 1, 2)]
tmp <- tmp[, common_name := gsub("\\.", "" ,gsub(" ", "_", common_name))]

tmp <- dcast(tmp, riv_mile ~ common_name, value.var = "survey", fun.aggregate = sum) ##1 = FJS, 2 = BSS, 3 = both

fwrite(tmp, paste0(PATH, "/bio_data/combined_survey_species_presence_in_rivmiles.csv"))


#attach pseudocoordinates -----
if(file.exists(paste0(PATH, "/inshore_offshore/combined_surveys_fjs_bss_coord_only.csv"))) {
  
  w_psc <- fread(paste0(PATH, "/inshore_offshore/combined_surveys_fjs_bss_coord_only.csv"), colClasses = c("uniq_id" = "character"))

  comb.surv <- w_psc[comb.surv, on = intersect(names(w_psc), names(comb.surv))]
  
  fwrite(comb.surv, paste0(PATH, "/inshore_offshore/combined_surveys_fjs_bss.csv"))
  
} else {
  ##this is also really inefficiently written..
  ps.coords <- fread(paste0(PATH, "/inshore_offshore/pseudo_coordinate_options.csv"))
  ps.coords <- st_as_sf(ps.coords, coords = c("X", "Y"), crs = 26918) |> 
    st_transform(4326) #convert to lat long
  
  ps.comb.surv <- data.table()
  for(i in 1:length(unique(comb.surv$uniq_id))) {
    
    id <- unique(comb.surv$uniq_id)[i]
    
    tmp <- comb.surv[uniq_id == id]
    
    if(length(unique(tmp$riv_mile)) > 1 | length(unique(tmp$strata_new)) > 1) { break } #stop check, should only be 1 per uniq_id
    
    if(unique(tmp$strata_new) %in% c("Shoals", "Bottom") & !is.na(unique(tmp$site)) & unique(tmp$site) != 5) { hab <- "mid" }
    if(unique(tmp$strata_new) == "Channel" | isTRUE(unique(tmp$site) == 5)) { hab <- "center" }
    if(unique(tmp$strata_new) == "beach" & is.na(unique(tmp$site))) { hab <- "edge" }
    if(unique(tmp$strata_new) == "beach" & !is.na(unique(tmp$site)) & unique(tmp$site) != 2) { hab <- "edge" }
    if(unique(tmp$strata_new) == "beach" & isTRUE(unique(tmp$site) == 2)) { hab <- "center" }
    # if(unique(tmp$strata_new) == "beach" & isTRUE(unique(tmp$site) %in% c(4, 6))) { break } #just going to use 4,6 code from fjs?
    
    #pull appropriate region
    tmp.psc <- ps.coords[ps.coords$riv_mile == unique(tmp$riv_mile) & ps.coords$pt_location == hab,]
    
    #if coordinates exist
    if(!all(is.na(tmp$lat))) {
      tmp.sf <- st_as_sf(tmp[, unique(.SD), .SDcols = c("long", "lat")], coords = c("long", "lat"), crs = 4326)
      psc.dist <- st_distance(tmp.sf, tmp.psc)
      
      tmp.psc <- tmp.psc[which.min(psc.dist),]
      
      tmp$dist2psc_m <- psc.dist[which.min(psc.dist)]
    }
    
    #if missing coordinates, attach by river mile and in river location 
    if(all(is.na(tmp$lat))) {
      if(hab %in% c("mid", "edge")) { #two options here:
        if(!unique(tmp$site) %in% c(1, 3, 4, 6)) { break } #checking
        
        #west of channel
        if(unique(tmp$site) %in% c(1, 4)) { 
          tmp.psc <- tmp.psc[which.min(st_coordinates(tmp.psc)[, "X"]),]
        }
        #east of channel
        if(unique(tmp$site) %in% c(3, 6)) { 
          tmp.psc <- tmp.psc[which.max(st_coordinates(tmp.psc)[, "X"]),]
        }
      }
    }
    
    #now add info to survey data
    tmp$long_psc <- st_coordinates(tmp.psc)[, "X"]
    tmp$lat_psc <- st_coordinates(tmp.psc)[, "Y"]
    
    ps.comb.surv <- rbindlist(list(ps.comb.surv, tmp), fill = TRUE)
    
    cat(round(nrow(ps.comb.surv)/nrow(comb.surv)*100, 2), "% ...\n")
    rm(hab, tmp.psc, tmp)
    
  }
  
  #save
  fwrite(ps.comb.surv, paste0(PATH, "/inshore_offshore/combined_surveys_fjs_bss.csv"))
  
  #for checking coordinates
  coord.only <- ps.comb.surv[, .(uniq_id, lat, long, lat_psc, long_psc)] |>
    unique()
  fwrite(coord.only, paste0(PATH, "/inshore_offshore/combined_surveys_fjs_bss_coord_only.csv"))
  
  ##compare original with add coordinate dt
  # test <- ps.comb.surv[, 1:6]
  # test2 <- comb.surv[, 1:6]
  # setkeyv(test, "uniq_id")
  # setkeyv(test2, "uniq_id")
  # all.equal(test, test2)
}

#extract bottom sediment type ----
comb.surv <- fread(paste0(PATH, "/inshore_offshore/combined_surveys_fjs_bss.csv"),
                   colClasses = c("uniq_id" = "character"))
id_check <- unique(comb.surv$uniq_id)

##get unique sites (for faster extraction)
comb.surv.ll <- comb.surv[, .(uniq_id, lat, long, lat_psc, long_psc)] |>
  unique()

##read in bottom type data
benth <- read_sf(paste0(PATH, "/env_data/gis_layers/Hudson_River_Estuary_Sediment_Type.shp"))
benth <- benth[, !names(benth) %in% c("AREA", "LEN")]
names(benth) <- c("benth_cat", "geometry")

#set function for duplicate extraction
get_benth_dat <- function(coord.dat, is_psc = TRUE) {
  
  ##extract for true // pseudocoordinates
  if(!is_psc) { coord.dat <- coord.dat[!is.na(lat), .(uniq_id, long, lat)] }
  if(is_psc)  { coord.dat <- coord.dat[,  .(uniq_id, long_psc, lat_psc)] }
  
  ##turn sf
  cs.ll.sf <- st_as_sf(coord.dat, coords = names(coord.dat)[grep("long|lat", names(coord.dat))], crs = 4326) |> 
    st_transform(st_crs(benth))
  
  ##get bottom type for all coords
  with_benth <- st_intersection(cs.ll.sf, benth)
  
  ##fill in those that are just outside the shapefile
  no_benth <- cs.ll.sf[!cs.ll.sf$uniq_id %in% with_benth$uniq_id, ]
  
  benth.ind <- st_nearest_feature(no_benth, benth)
  no_benth$benth_cat <- benth[benth.ind,]$benth_cat
  
  ##join, fix names, return
  coord.dat <- do.call(rbind, list(with_benth, no_benth))
  
  if(is_psc) { names(coord.dat)[names(coord.dat) == "benth_cat"] <- "benth_cat_psc" }
  
  return(st_drop_geometry(coord.dat))
  
}

cs.ll.sf <- get_benth_dat(comb.surv.ll, is_psc = FALSE)
cs.ps.ll.sf <- get_benth_dat(comb.surv.ll, is_psc = TRUE)

##add in var to survey dat
comb.surv <- merge(comb.surv, cs.ll.sf, all.x = TRUE)
comb.surv <- merge(comb.surv, cs.ps.ll.sf, all.x = TRUE)

# comb.surv[!is.na(lat) & is.na(benth_cat), flag_benth_ll := 1]
# comb.surv[is.na(benth_cat_psc), flag_benth_psc := 1]
# check <- comb.surv[, .(uniq_id, lat, long, lat_psc, long_psc, flag_benth_ll, flag_benth_psc)] |>
#   unique()
# fwrite(check, paste0(PATH, "/inshore_offshore/check_benthic.csv"))
# comb.surv[, c("flag_benth_ll", "flag_benth_psc") := NULL]

fwrite(comb.surv, paste0(PATH, "/inshore_offshore/combined_surveys_fjs_bss.csv"))

# get diff from solar noon ----
library(suntools)
comb.surv <- fread(paste0(PATH, "/inshore_offshore/combined_surveys_fjs_bss.csv"))

##get avg coord for river
ll <- comb.surv[, lapply(.SD, mean, na.rm = TRUE), .SDcols = c("long", "lat")]
ll <- as.matrix(ll[1])

cs.time <- comb.surv[, unique(.SD), .SDcols = c("date")]
cs.time[, date := as.POSIXct(date, tz = "America/New_York")]

cs.time <- as.data.table(solarnoon(ll, cs.time$date, POSIXct.out = TRUE))

cs.time[, `:=` (date = as.Date(time), day_frac = NULL)]
setnames(cs.time, "time", "solar_noon")

comb.surv <- cs.time[comb.surv, on = "date"]

comb.surv[time != "", date_time := as.POSIXct(paste0(date, " ", time), tz = "America/New_York")]
comb.surv[time != "", solar_noon_diff := difftime(date_time, solar_noon, units = "hours")]
comb.surv[, solar_noon_diff := as.numeric(solar_noon_diff)] #for some reason it doesn't want to convert if done within the same line above
comb.surv[, `:=` (solar_noon_diff = round(solar_noon_diff, 3), date_time = NULL)]

setcolorder(comb.surv, c(setdiff(names(comb.surv), "solar_noon"), "solar_noon"))

##all time columns need to be saved as character (otherwise can do some funky things based on time zone)
comb.surv[, c("time", "solar_noon") := lapply(.SD, as.character), .SDcols = c("time", "solar_noon")]

fwrite(comb.surv, paste0(PATH, "/inshore_offshore/combined_surveys_fjs_bss.csv"))

# estimate area swept ----
comb.surv <- fread(paste0(PATH, "/inshore_offshore/combined_surveys_fjs_bss.csv"))

##separate surveys bc way it is estimated is diff
bss <- comb.surv[survey == "bss"]
fjs <- comb.surv[survey == "fjs"]
# fjs <- fjs[!gear_def == "2m Beam Trawl"]

##bss area swept = volume bc those were estimated from area to begin with, and we assumed a height of 1 (see bss_area_swept.R)
bss[, area_swept := volume]

##fjs area swept = tow speed * tow duration * gear mouth width
dist.m <- fread(paste0(PATH, "/hrbmp_database/FSS_Level3.csv"))
dist.m <- dist.m[, .(Uniq_ID, DATE, GEAR, TOW_SPD, DURATION)]
names(dist.m) <- c("uniq_id", "date", "gear_code", "tow_speed", "duration")

#add gear name
gear.code <- fread(paste0(PATH, "/hrbmp_database/FSS_Gear_Codes.csv"))
setnames(gear.code, c("Code", "Gear"), c("gear_code", "gear_def"))
dist.m <- gear.code[dist.m, on = "gear_code"]

#add gear width, widths confirmed in HRBMP final report
dist.m[, gear_width := fifelse(gear_def == "3m Beam Trawl", 3,
                               fifelse(gear_def == "1m Tucker Trawl", 1,
                                       fifelse(gear_def == "2m Beam Trawl", 2,
                                               fifelse(gear_def == "1m Epibenthic Sled", 1, 
                                                       NA))))]

## fix NAs from missing tow speed and/or duration
dist.m[, date := {
  date_sub = sub(" .*", "", date)
  as.Date(date_sub, format = "%m/%d/%y")
}]
dist.m[, `:=` (year = year(date),
               month = month(date), 
               jul_day = yday(date))][, 
                                      jul_week := ceiling(jul_day/7)]

summ_levels <- list(
  c("date", "gear_def"),
  c("jul_week", "year", "gear_def"),
  c("month", "year", "gear_def"),
  c("year", "gear_def"),
  c("gear_def")
)

for (g in summ_levels) {
  dist.m[, tow_speed := fifelse(is.na(tow_speed) | is.nan(tow_speed),
                                mean(tow_speed, na.rm = TRUE),
                                tow_speed),
         by = g]
}
##the only ones left are 2m beam trawls

fwrite(dist.m, paste0(PATH, "/inshore_offshore/fjs_area_swept_estimates.csv"))

#calc area swept
fjs <- dist.m[fjs, on = "uniq_id"]
fjs[, area_swept := {
  dist_traveled = tow_speed * duration * 60 #(m/s) * (min) * (60 s / min) = m
  dist_traveled * gear_width # m * m
}]

fjs <- fjs[, .SD, .SDcols = names(bss)]

#add back together
comb.surv <- rbindlist(list(fjs, bss))
fwrite(comb.surv, paste0(PATH, "/inshore_offshore/combined_surveys_fjs_bss.csv"))

if(FALSE) {
  ##compare survey, gear area swepts
  s.comb.surv <- comb.surv[, unique(.SD), .SDcols = c("uniq_id", "survey", "riv_mile", "year", "gear_def", "area_swept", "volume")]
  
  summ_levels = list(
    c("survey"),
    c("survey", "year"),
    c("survey", "riv_mile"),
    c("survey", "riv_mile", "year"),
    c("gear_def"),
    c("gear_def", "year"),
    c("gear_def", "riv_mile"),
    c("gear_def", "riv_mile", "year")
  )
  
  summ_swept <- vector("list", length(summ_levels))
  for(i in seq_along(summ_levels)) {
    
    summ_swept[[i]] <- s.comb.surv[, .(mn_area = mean(area_swept),
                                       sd_area = sd(area_swept),
                                       mn_vol = mean(volume),
                                       sd_vol = sd(volume)),
                                 by = eval(summ_levels[[i]])]
    
    summ_swept[[i]]$summ_level <- paste(summ_levels[[i]], collapse = ", ")
    
  }
  
  summ_swept_dt <- rbindlist(summ_swept, fill = T)
  
  fwrite(summ_swept_dt, paste0(PATH, "/inshore_offshore/area_swept_summaries.csv"))
  
  ggplot(data = s.comb.surv, aes(x = gear_def, y = area_swept)) +
    geom_boxplot(outliers = T) #+
    # geom_jitter(width = 0.1, alpha = 0.6)
  
  ggplot(data = summ_swept[[6]]) +
    geom_path(aes(x = year, y = mn_area, colour = gear_def)) +
    geom_point(aes(x = year, y = mn_area, color = gear_def))
  ggsave(paste0(PATH, "/figures/inshore_offshore_index/area_swept_annual_gear_comp.png"), width = 8, height = 5)
  
}

# get river depth from DEM ----
library(terra)
comb.surv <- fread(paste0(PATH, "/inshore_offshore/combined_surveys_fjs_bss.csv"))

bath <- rast(paste0(PATH, "/env_data/gis_layers/30m_grid/dem/estuary.dem"))

if(any(grepl("dpth_est", names(comb.surv)))) {
  comb.surv[, `:=` (riv_dpth_est = NULL, riv_dpth_est_psc = NULL, 
                    flag_coord_issue = NULL, flag_coord_issue_psc = NULL, 
                    bath_cell_dist_m = NULL, bath_cell_dist_m_psc = NULL)]
}

name.order <- names(comb.surv)

o.coords <- comb.surv[, unique(.SD), .SDcols = c("uniq_id", "survey", "lat", "long")][!is.na(long)]
p.coords <- comb.surv[, unique(.SD), .SDcols = c("uniq_id", "survey", "lat_psc", "long_psc")]

##for fixing coordinates outside the depth layer
bath.pts <- as.points(bath, na.rm = TRUE)
bath.pts <- st_as_sf(bath.pts)

dpth <- lapply(list(o.coords, p.coords), function(DT) {
  DT.sf <- st_as_sf(DT, coords = c(grep("long", names(DT), value = T), grep("lat", names(DT), value = T)), crs = 4269) |>
    st_transform(st_crs(bath))
  
  bath.out <- extract(bath, vect(DT.sf), method = "simple")
  DT$riv_dpth_est <- abs(bath.out$estuary)
  # DT[, riv_dpth_est := fifelse(is.na(riv_dpth_est), 0, riv_dpth_est)]
  DT[, flag_coord_issue := fifelse(is.na(riv_dpth_est), TRUE, FALSE)]
  
  ##deal with coordinates outside the data layer
  DT.na <- DT[flag_coord_issue == TRUE]
  DT.notna <- DT[flag_coord_issue == FALSE]
  
  DT.sf <- st_as_sf(DT.na, coords = c(grep("long", names(DT), value = T), grep("lat", names(DT), value = T)), crs = 4269) |>
    st_transform(st_crs(bath))
  
  #get nearest cells (previously converted to points)
  near.indx <- st_nearest_feature(DT.sf, bath.pts)
  
  #get matching raster points and values
  near.cells <- bath.pts[near.indx, ]
  
  #get distance
  dists <- st_distance(DT.sf, near.cells, by_element = TRUE)
  
  #get depth, distance, and add back in
  DT.na$riv_dpth_est <- abs(near.cells$estuary)
  DT.na$bath_cell_dist_m <- as.numeric(dists)
  
  if(nrow(DT.na) + nrow(DT.notna) != nrow(DT)) { stop("Something went wrong. Number of rows for fixing missing depth is not what is expected.") }
  
  DT <- rbindlist(list(DT.na, DT.notna), fill = TRUE)
  
  if(any(grepl("_psc", names(DT)))) {
    setnames(DT, c("riv_dpth_est", "flag_coord_issue", "bath_cell_dist_m"), c("riv_dpth_est_psc", "flag_coord_issue_psc", "bath_cell_dist_m_psc"))
  }
  
  return(DT)
})

#attach to full data
for(i in seq_along(dpth)) {
  comb.surv <- dpth[[i]][comb.surv, on = intersect(names(dpth[[i]]), names(comb.surv))]
}

##check differences
if(FALSE) {
  tmp <- comb.surv[, unique(.SD), .SDcols = c("uniq_id", "survey", "lat", "long", "lat_psc", "long_psc", "strata_new", 
                                              "riv_dpth", "riv_dpth_est", "riv_dpth_est_psc", 
                                              "flag_coord_issue", "flag_coord_issue_psc", "bath_cell_dist_m", "bath_cell_dist_m_psc")]
  tmp[, `:=` (diff_dpth = abs(riv_dpth - riv_dpth_est),
              diff_dpth_psc = abs(riv_dpth - riv_dpth_est_psc))]
  tmp[, .(mn_dpth = mean(riv_dpth_est, na.rm = T),
          min_dpth = min(riv_dpth_est, na.rm = T),
          max_dpth = max(riv_dpth_est, na.rm = T)), by = survey]
  fwrite(tmp, paste0(PATH, "/inshore_offshore/river_dpth_bath_coord_check.csv"))
  # tmp[survey == "fjs" & riv_dpth_est == 0] ##these are the biggest problems
}

##adjust river depth and sample depth where necessary (right now only bss)
comb.surv[, `:=` (riv_dpth = fifelse(survey == "bss", riv_dpth_est, riv_dpth),
                  sam_dpth = fifelse(survey == "bss" & riv_dpth_est <= 5, riv_dpth_est, 
                                     fifelse(survey == "bss" & riv_dpth_est > 5, 5, sam_dpth)))]
comb.surv[, `:=` (riv_dpth = fifelse(is.na(riv_dpth), riv_dpth_est_psc, riv_dpth),
                  sam_dpth = fifelse(survey == "bss" & is.na(sam_dpth) & riv_dpth_est_psc <= 5, riv_dpth_est_psc,
                                     fifelse(survey == "bss" & is.na(sam_dpth) & riv_dpth_est_psc > 5, 5, sam_dpth)))]

#reorganize
setcolorder(comb.surv, name.order)

#save
fwrite(comb.surv, paste0(PATH, "/inshore_offshore/combined_surveys_fjs_bss.csv"))

# SPECIES SPECIFIC DATAPREP ----
#load in cleaned/filtered data (all steps above applied)
io.dat <- fread(paste0(PATH, "/inshore_offshore/combined_surveys_fjs_bss.csv"))

##subset years ----
# io.dat[, .(min_yr = min(year), max_yr = max(year)), by = c("survey", "gear_def")]
io.dat <- io.dat[year <= 2014 & year >= 1980] ##most consistent overlap gears x years

##set FJS tows with no assigned river depth (coordinate errors) to be pseudocoords ----
io.dat <- io.dat[(survey == "fjs" & flag_coord_issue == TRUE), `:=` (lat = NA, long = NA)]

##subset individual survey types ----
io.dat.l <- split(io.dat, io.dat$survey)

io.dat.l[[3]] <- io.dat
names(io.dat.l)[3] <- "all"

##set target spp ----
target.species <- c("Alewife", "American Shad", "Blueback Herring", "Striped Bass")
target.surveys <- names(io.dat.l)

#apply data filters to each species for each survey set
for(s in target.species) {
  
  io.dat.spp <- lapply(io.dat.l, function(DT, spp = s) {
    
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
    
    ## use pseudo coords for missing lat long // extracted vars ----
    DT.spp[, `:=` (lat = fifelse(is.na(lat), lat_psc, lat),
                       long = fifelse(is.na(long), long_psc, long),
                       benth_cat = fifelse(benth_cat == "", benth_cat_psc, benth_cat))]
    DT.spp[, c("lat_psc", "long_psc", "benth_cat_psc") := NULL]
    
    #return finalized species dataset
    return(DT.spp)
    
  })
  
  spp.dir <- paste0(PATH, "/inshore_offshore/", gsub(" ", "_", s))
  
  if(!dir.exists(spp.dir)) { dir.create(spp.dir) }
  
  for(surv in target.surveys) {
    
    ##add some checks JIC
    if(all(io.dat.spp[[surv]]$common_name != s)) { stop("Something went wrong subsetting the species data for ", surv, "...") }
    
    if(surv != "all") {
      if(all(io.dat.spp[[surv]]$survey != surv)) { stop("Something went wrong subsetting the surveys for fjs or bss...") }
    } else {
      if(!all(io.dat.spp[[surv]]$survey %in% c("fjs", "bss"))) { stop("Either FJS or BSS is missing from combined dataset...") }
    }
    
    ##fix tows with more than 1 row per species ----
    if(length(unique(io.dat.spp[[surv]]$uniq_id)) != nrow(io.dat.spp[[surv]])) {

      keep_cols <- setdiff(names(io.dat.spp[[surv]]), c("uniq_id", "ct_yoy"))  # everything else
      
      dup_rows <- io.dat.spp[[surv]][io.dat.spp[[surv]][, .I[.N > 1], by = uniq_id]$V1]
      
      dup_rows <- dup_rows[, c(
        list(ct_yoy = sum(ct_yoy)), 
        lapply(.SD, first)            
      ), by = uniq_id, .SDcols = keep_cols]
      
      io.dat.spp[[surv]] <- rbindlist(list(io.dat.spp[[surv]][!uniq_id %in% dup_rows$uniq_id], dup_rows), use.names = TRUE)
      
    }
    
    fwrite(io.dat.spp[[surv]], paste0(spp.dir, "/", gsub(" ", "_", s), "_catch_dat_", surv, ".csv"))
    
    message(surv, " data for ", s, " saved...")
  }
  
}

rm(list = ls())

# final data set - up (mostly remove missing covariate rows) ----
library(data.table)

PATH <- getwd()

target.species <- c("Alewife", "American Shad", "Blueback Herring", "Striped Bass")
target.surveys <- c("all", "bss", "fjs")

##set covariates ----
cov.col.names <- c("gear_def", "sam_dpth", "solar_noon_diff", "riv_dpth", "benth_cat")
 
summ_datasets <- data.table()

for(spp in target.species) {
  
  spp.wd <-  paste0(PATH, "/inshore_offshore/", gsub(" ", "_", spp))
  
  for(surv in target.surveys) {
    io.dat.spp <- fread(paste0(spp.wd, "/", gsub(" ", "_", spp), "_catch_dat_", surv, ".csv"))
    
    # if(any(io.dat.spp[, is.na(.SD), .SDcols = intersect(names(io.dat.spp), cov.col.names)])) { cat("CHECK SPECIES DATA, WHY ARE THERE NAs IN THE COVs?\n")}
    missing_dat_rows <- io.dat.spp[!complete.cases(io.dat.spp[, ..cov.col.names])] # returns rows with at least one NA in cov columns
    
    ## remove obs if missing data ----
    io.dat.spp <- na.omit(io.dat.spp, cols = cov.col.names)
    
    ## update benthic categories ----
    ## benthic category causes errors in the model, need to combine rare levels
    #just doing three cats
    if(any(grepl("benth_cat", cov.col.names))) {
      io.dat.spp[, benth_cat := fifelse(benth_cat %in% c("sandy mud", "gravelly mud", "mud"), "mud", 
                                        fifelse(benth_cat %in% c("gravelly sand", "muddy sand", "sand"), "sand",
                                                fifelse(benth_cat %in% c("sandy gravel", "muddy gravel", "gravel"), "gravel",
                                                        "other")))] ##none should be in this category otherwise something is missing
      
      
      #check for too few per group again
      min_n <- 10
      n_group <- io.dat.spp[, .N, by = benth_cat][order(-N)]
      
      #get small groups
      small_groups <- n_group[N < min_n]$benth_cat
      
      benth_check <- length(small_groups) > 0
      
      while(length(small_groups) != 0) {
        if(length(small_groups) == 1) {
          
          next_grp <- n_group$benth_cat[nrow(n_group)-1]
          
          warning("Making new benthic category: ", paste(c(next_grp, small_groups), collapse = "-"))
          
          io.dat.spp[benth_cat %in% c(small_groups, next_grp), benth_cat := paste(c(next_grp, small_groups), collapse = "-")]
          
        } else {
          
          warning("Making new benthic category: ", paste(small_groups, collapse = "-"))
          io.dat.spp[benth_cat %in% small_groups, benth_cat := paste(small_groups, collapse = "-")]
          
        } 
        
        #check again
        n_group <- io.dat.spp[, .N, by = benth_cat][order(-N)]
        
        #get small groups
        small_groups <- n_group[N < min_n]$benth_cat
      }
    }
    
    #summarize final counts
    tmp <- data.table(species = spp,
                      survey = surv,
                      total_tows = nrow(io.dat.spp),
                      mn_pos_catch = mean(io.dat.spp[ct_yoy > 0]$ct_yoy),
                      total_zero_catch = nrow(io.dat.spp[ct_yoy == 0]),
                      tows_per_fjs = nrow(io.dat.spp[survey == "fjs"]),
                      tows_per_bss = nrow(io.dat.spp[survey == "bss"]),
                      rows_removed_for_na = nrow(missing_dat_rows),
                      collapsed_benth_cat = benth_check)
    
    summ_datasets <- rbindlist(list(summ_datasets, tmp))
    
    fwrite(io.dat.spp, paste0(spp.wd, "/", gsub(" ", "_", spp), "_vast_prep_dat_", surv, ".csv"))
  }
}

fwrite(summ_datasets, paste0(PATH, "/inshore_offshore/species_datasets_final_counts_for_model_input.csv"))







