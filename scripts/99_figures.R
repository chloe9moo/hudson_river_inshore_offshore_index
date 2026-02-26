## main figures ##

library(data.table); library(sf); library(ggplot2)

PATH <- getwd()
fig.save.path <- paste0(PATH, "/inshore_offshore/amarel_cluster/figures/")
#which run location
model.run.location = "amarel_cluster/" #""
#which type
model.run.type = "knot comp"
#target spp
spp.l = c("Alewife", "Striped Bass", "American Shad", "Blueback Herring")


#load in full dataset
io.dat <- fread(paste0(PATH, "/inshore_offshore/combined_surveys_fjs_bss.csv"))
io.dat[is.na(long), `:=` (long = long_psc, lat = lat_psc)]

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

#which sections to run?
fig_n1 <- FALSE
fig_0 <- F
fig_1 <- T
fig_2 <- T
fig_3 <- T
fig_4 <- T
fig_s1 <- T
pred_splines_check <- F
dharma <- T

#gg themes ----
io_theme <- theme(panel.border = element_rect(color = "black", fill = NA),
                  panel.background = element_rect(fill = "white"),
                  panel.grid.major = element_line(color = "lightgray"),
                  strip.background = element_rect(color = "black", fill = "gray"),
                  axis.text.x = element_text(angle = 47, vjust = 1, hjust = 1),
                  axis.text = element_text(size = 12, color = "black"),
                  strip.text = element_text(size = 12, face = "bold"),
                  axis.title = element_text(size = 12),
                  legend.text = element_text(size = 12, color = "black"))
source(paste0(PATH, "/scripts/inshore_offshore/XX_colors.R"))

#misc
if(grepl("knot", model.run.type)) {
  mod.type = "knot_comp"
}
if(grepl("cov", model.run.type)) {
  mod.type = "cov_comp"
}

# figure -1: survey distribution map ----
if(fig_n1) {
#river polygon 
# hr.poly <- read_sf(paste0(PATH, "/env_data/gis_layers/NOAA_rivermile_polygon_edited.shp"))
hr.poly <- read_sf(paste0(PATH, "/env_data/gis_layers/HudsonRiverKms_Poly.shp"))
  
#river line (for directionality)
nhd <- read_sf(paste0(PATH, "/env_data/gis_layers/hudson_watershed_stream_clip.shp"))
nhd <- nhd[nhd$GNIS_NAME == "Hudson River" & !is.na(nhd$GNIS_NAME), ] |>
  st_transform(st_crs(hr.poly)) |>
  st_zm()

io.id <- io.dat[, unique(.SD), .SDcols = c("uniq_id", "long", "lat", "riv_mile", "survey")]
io.id[, coord_id := .GRP, by = c("long", "lat")]

io.sf <- io.id[, unique(.SD), .SDcols = c("coord_id", "long", "lat")]
io.sf <- st_as_sf(io.sf, coords = c("long", "lat"), crs = 4326) |>
  st_transform(st_crs(hr.poly))

get_rel_widths <- function(pt, crs.in = st_crs(hr.poly)) {
  #get relative location along width of river
  if(FALSE) {
    i <- 29605
    pt <- st_geometry(io.sf[i,])
    crs.in = st_crs(hr.poly)
    # targ.rm <- io.sf[i,]$riv_mile
  }

  pt <- st_sfc(pt, crs = crs.in)
  
  # get nhd stream in river mile/river km
  sub.hr <- st_filter(hr.poly, pt)
  if(nrow(sub.hr) == 0) {
    n.i <- st_nearest_feature(pt, hr.poly)
    sub.hr <- hr.poly[n.i,]
  }
  
  if(nrow(sub.hr) > 1 | nrow(sub.hr) == 0) { stop("Too many/too few river polygons pulled.") }
  if(any(names(sub.hr) == "Rmile")) {
    targ.rm <- sub.hr$Rmile
  }
  if(any(names(sub.hr) == "RKm")) {
    targ.rm <- sub.hr$RKm
  }
  if(is.null(targ.rm)) {
    stop("Something went wrong pulling target river mile/km")
  }
  
  if(any(names(sub.hr) == "Rmile")) {
    
    #cut nhd stream to river mile polygon
    nhd.rmi <- suppressWarnings(st_intersection(nhd, sub.hr)) |>
      st_union() |>
      st_cast("MULTILINESTRING") |>
      st_line_merge()
    
    # get start and end of the nhd line
    start <- st_coordinates(lwgeom::st_startpoint(nhd.rmi))
    end <- st_coordinates(lwgeom::st_endpoint(nhd.rmi))
    
    # turn those points into a single straight line
    nhd.c <- start[which(!start %in% end)]
    end <- end[which(!end %in% start)]
    
    nhd.c <- rbind(nhd.c, end)
    line.sf <- st_sfc(st_linestring(nhd.c), crs = crs.in) ##make spatial
    
    # get perpendicular line (90 d rotation) to nhd
    delta.x <- nhd.c[2, 1] - nhd.c[1, 1]
    delta.y <- nhd.c[2, 2] - nhd.c[1, 2]
    
    # find nearest point along line from target point
    pt.line <- st_nearest_points(pt, line.sf)
    pt.line <- st_cast(pt.line, "POINT")[2] |>
      st_coordinates()
    
    # use this to make a perpendicular transect
    p.90 <- rbind(
      pt.line[, c("X", "Y")] + 150000 * c(-delta.y, delta.x) / sqrt(delta.x^2 + delta.y^2),
      pt.line[, c("X", "Y")] - 150000 * c(-delta.y, delta.x) / sqrt(delta.x^2 + delta.y^2)
    )
    
    transect <- st_sfc(st_linestring(p.90), crs = crs.in)
  }
  
  if(any(names(sub.hr) == "RKm")) {
    
    #make transect along center of river km polygon (these are divided differently than the river miles)
    start <- c(st_coordinates(pt)[, "X"] - 10000, st_coordinates(pt)[, "Y"])
    end <- c(st_coordinates(pt)[, "X"] + 10000, st_coordinates(pt)[, "Y"])
    
    transect <- st_sfc(st_linestring(rbind(start, end)), crs = crs.in) ##make spatial
    
  }

  #intersect transect with polygon to get the two shores
  cross <- suppressWarnings(st_intersection(sub.hr, transect))

  #find the endpoints (shore points) on the transect
  shore_pts <- suppressWarnings(st_cast(cross, "POINT"))

  if(nrow(shore_pts) != 2 & nrow(shore_pts) != 0) { 
    shore_pts <- shore_pts[c(which.min(st_coordinates(shore_pts)[,"X"]),
                             which.max(st_coordinates(shore_pts)[,"X"])),]
  }
  if(nrow(shore_pts) == 0) {
    out <- data.table(assigned_riv_km = targ.rm, river_width = NA, rel_position = NA)
    return(out)
  }
  
  #order from west to east
  coords <- st_coordinates(shore_pts)
  order_x <- order(coords[, "X"])
  
  shore_pts <- shore_pts[order_x,]

  #get width and rel. distance along width for coordinate
  dist_total <- as.numeric(st_distance(shore_pts[1,], shore_pts[2,], by_element = TRUE))
  dist_pt <- as.numeric(st_distance(shore_pts[1,], pt))
  rel_position <- dist_pt / dist_total
  
  ##testing
  if(FALSE) {
    # tmp <- st_as_sf(io.pos[assigned_riv_km == 175], coords = c("long", "lat"), crs = 4326)
    ggplot() +
      # geom_sf(data = hr.poly[hr.poly$RKm %in% c(56, 57, 58),]) +
      geom_sf(data = sub.hr) +
      # geom_sf(data = nhd.rmi, color = "blue") +
      # geom_sf(data = line.sf) +
      # geom_sf(data = pt.line) +
      geom_sf(data = pt, shape = 4) +
      # geom_sf(data = tmp, shape = 4) +
      # geom_sf(data = transect) +
      geom_sf(data = shore_pts)
  }
  
  if(any(names(sub.hr) == "Rmile")) {
    out <- data.table(assigned_riv_mile = targ.rm, river_width = dist_total, rel_position = rel_position)
  }
  if(any(names(sub.hr) == "RKm")) {
    out <- data.table(assigned_riv_km = targ.rm, river_width = dist_total, rel_position = rel_position)
  }

  return(out)
}

# io.pos <- lapply(st_geometry(io.sf), get_rel_widths)

##looping to track progress instead..
if(FALSE) {
  io.pos <- vector("list", length = nrow(io.sf))
  for(i in 1:nrow(io.sf)) {
    io.pos[[i]] <- get_rel_widths(st_geometry(io.sf[i,]))
    
    cat(sprintf("\r%-20s", sprintf("%.0f%% complete...", 100 * i / nrow(io.sf))))
    
    if(i == nrow(io.sf)) { cat("\nDone!\n")}
  }
  
  io.pos <- rbindlist(io.pos)
  
  io.pos[which(round(io.pos$rel_position, 2) > 1)]
  io.pos[is.na(io.pos$rel_position)]
  
  fwrite(io.pos, paste0(PATH, "/inshore_offshore/io_survey_relative_widths_rkm.csv"))
  io.pos <- fread(paste0(PATH, "/inshore_offshore/io_survey_relative_widths_rkm.csv"))
  
  io.sf$river_width <- io.pos$river_width
  io.sf$rel_position <- io.pos$rel_position
  
  if(any(names(hr.poly) == "RKm")) {
    io.sf$assigned_riv_km <- io.pos$assigned_riv_km
    fwrite(st_drop_geometry(io.sf), paste0(PATH, "/inshore_offshore/io_survey_relative_widths_rkm.csv"))
  } else {
    io.sf$assigned_riv_mile <- io.pos$assigned_riv_mile
    fwrite(st_drop_geometry(io.sf), paste0(PATH, "/inshore_offshore/io_survey_relative_widths.csv"))
  }
}


io.pos <- fread(paste0(PATH, "/inshore_offshore/io_survey_relative_widths_rkm.csv"))

##fix the 0 - 12 river miles to be individual sections ----
if(any(names(hr.poly) == "Rmile")) {
  hr.0 <- hr.poly[hr.poly$Rmile == 0,]
  io.0 <- io.sf[io.sf$assigned_riv_mile == 0,]
  
  #get bounding box of zero section
  bbox <- st_bbox(hr.0)
  
  #make 12 equal size breaks up river
  lat.breaks <- seq(bbox$ymin, bbox$ymax, length.out = 13)  # 12 bands = 13 edges
  
  #get river width (ish)
  x.min <- bbox$xmin
  x.max <- bbox$xmax
  
  #make polygons to cut up the river with
  slices <- lapply(seq_len(length(lat.breaks) - 1), function(i) {
    st_polygon(list(rbind(
      c(x.min, lat.breaks[i]),
      c(x.max, lat.breaks[i]),
      c(x.max, lat.breaks[i+1]),
      c(x.min, lat.breaks[i+1]),
      c(x.min, lat.breaks[i])
    )))
  })
  
  slices.sf <- st_sf(geometry = st_sfc(slices, crs = st_crs(hr.0)))
  
  #cut river
  hr.0.new <- st_intersection(hr.0, slices.sf)
  hr.0.new$Rmile <- seq_len(nrow(hr.0.new)) - 1  #assign river mile (ish!)
  
  #assign new river mile to io.sf obj
  for(i in seq_len(nrow(hr.0.new))) {
    tmp <- hr.0.new[i,]
    
    tmp.io <- st_filter(io.0, tmp)
    
    io.pos[, assigned_riv_mile := fifelse(coord_id %in% tmp.io$coord_id, unique(tmp$Rmile), assigned_riv_mile)]
  }
  
  #add new split polygon back to full river
  hr.poly <- do.call(rbind, list(hr.poly[hr.poly$Rmile != 0,], hr.0.new))
  
  # tmp <- io.pos[rel_position > 100, uniq_id]
  # tmp <- io.sf[io.sf$uniq_id %in% tmp,]
  
  ggplot() +
    # geom_sf(data = slices.sf) +
    # geom_sf(data = hr.0, fill = NA) +
    # geom_sf(data = hr.0.new, aes(fill = as.factor(Rmile))) +
    # geom_sf(data = hr.poly[hr.poly$Rmile == 87,], aes(fill = as.factor(Rmile))) +
    # geom_sf(data = tmp[tmp$riv_mile == 87,]) +
    geom_sf(data = io.0) +
    geom_sf(data = tmp, color = "red", fill = NA) +
    geom_sf(data = tmp.io, color = "red") +
    # geom_sf(data = io.sf[io.sf$uniq_id %in% io.pos[riv_mile == 11, uniq_id],], color = "blue")
    theme(legend.position = "none")
}

##summarize relative positions for plotting ----
io.pos <- io.pos[io.id, on = "coord_id"]

#subset years
# io.dat[, .(min_yr = min(year), max_yr = max(year)), by = c("survey", "gear_def")]
io.dat.sub <- io.dat[year <= 2014 & year >= 1980, uniq_id]
io.pos <- io.pos[uniq_id %in% io.dat.sub]

#'bin' the data
io.pos[, rel_position := fifelse(rel_position > 1, 1, 
                                 fifelse(rel_position < 0, 0, rel_position))]
io.pos[, width_bin := round(rel_position, digits = 1)]

#summarize
io.summ <- io.pos[, .(num_hauls = .N), by = c("survey", "assigned_riv_km", "width_bin")]
io.summ[, survey := fifelse(survey == "fjs", "offshore", "inshore")]

##plot ----
p2 <- ggplot(data = io.summ, aes(x = width_bin, y = assigned_riv_km, fill = num_hauls)) +
  geom_tile() +
  scale_y_continuous(n.breaks = 10, expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0),    
                     breaks = seq(0, 1, by = 0.25),
                     labels = c("0", "", "0.5", "", "1")) +
  scale_fill_viridis_c(trans = "log", breaks = c(1, 10, 100, 1000), limits = c(1, 1000),
                       labels = scales::comma_format(),  # or use label = log10 if you prefer
                       name = "Total\nhauls\n(log scale)") +
  facet_wrap(~ survey) +
  labs(x = "relative position along\nwidth of the river", y = "river kilometer") +
  theme(panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "lightgray"),
        strip.background = element_rect(color = "black", fill = "gray"),
        # axis.text.x = element_text(angle = 47, vjust = 1, hjust = 1),
        axis.text = element_text(size = 12, color = "black"),
        strip.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12, color = "black")) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 0),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 16, vjust = 1, hjust = 0),
        legend.title = element_text(size = 16))
p2

bbox <- st_bbox(st_transform(hr.poly, crs = 4326))
xmin <- bbox[1] 
xmax <- bbox[3]
xbreaks <- seq(xmin, xmax, length.out = 3)
xbreaks <- round(xbreaks, 2)

#make labels for river polygon
hr.poly <- st_transform(hr.poly, 4326)
hr.poly$r_km_label <- fifelse(hr.poly$RKm %in% seq(0, 225, 25), hr.poly$RKm, NA)
hr.center <- st_centroid(hr.poly)
hr.center$X <- st_coordinates(hr.center)[, 'X']
hr.center$Y <- st_coordinates(hr.center)[, 'Y']
hr.center <- hr.center[which(!is.na(hr.center$r_km_label)),]

p1 <- ggplot() +
  geom_sf(data = hr.poly) +
  ggrepel::geom_text_repel(data = hr.center, aes(x = X, y = Y, label = r_km_label), 
                           xlim = c(bbox[3]-0.05, bbox[3]+0.1), min.segment.length = 0.1, size = 5.3) +
  scale_x_continuous(breaks = xbreaks, limits = c(bbox[1], bbox[3]+0.1)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "", title = " ") +
  # coord_sf(xlim = c(xmin, xmax), expand = FALSE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text = element_text(size = 16, color = "black"),
        panel.border = element_rect(color = "#d62828", fill = NA, linewidth = 1.5),
        plot.margin = margin(0, 0, 0, 0, "cm"))

plower <- ggpubr::ggarrange(p1, p2, nrow = 1, labels = c("B", "C"), font.label = list(size = 21), widths = c(0.7, 1), heights = c(1.2, 1))
plower
ggsave(paste0(fig.save.path, "/study_extent_with_summ_tows.png"), plot = plower, width = 10, height = 10, bg = "white")

### inset ----
north_america <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
north_america <- subset(north_america, continent == "North America") |> 
  st_transform(., crs = 4269)
north_america <- north_america[north_america$subunit %in% c("United States", "Canada"),] |>
  st_union()

ny.sf <- st_as_sf(maps::map("state", c("new york"), fill = T, plot = F)) |> 
  st_transform(., crs = 4269)
ny.bb <- st_bbox(ny.sf)
usa.sf <- st_as_sf(maps::map("state", fill = T, plot = F)) |> 
  st_transform(., crs = 4269)

# state.sf <- st_as_sf(maps::map("state", c("new york", "new jersey", "connecticut", "massachusetts", "vermont", "pennsylvania"), fill = T, plot = F)) |> 
#   st_transform(., crs = 4269)
bbox.sf <- st_as_sfc(bbox)

labels <- data.table(name = c("Canada", "New York, USA"),#, "Atlantic Ocean"),
                     X = c(-79.5, -75.8),#, -72),
                     Y = c(44.5, 42.8))#, 40.5))

p3 <- ggplot() +
  geom_sf(data = north_america, fill = "lightgray") +
  geom_sf(data = usa.sf, fill = "lightgrey") +
  geom_sf(data = ny.sf[ny.sf$ID == "new york",], fill = "darkgrey") +
  geom_sf(data = bbox.sf, fill = NA, color = "#d62828", linewidth = 1.5) +
  geom_sf(data = hr.poly, fill = "blue", color = "blue") +
  # geom_text(data = labels, aes(x = X, y = Y, label = name), fontface = "italic") +
  coord_sf(xlim = c(ny.bb[1]-1, ny.bb[3]+1), ylim = c(ny.bb[2], ny.bb[4])) +
  theme_bw() + labs(x = "longitude", y = "latitude") +
  theme(axis.text = element_text(color = "black", size = 16),
        axis.title = element_text(size = 16))
ggsave(paste0(fig.save.path, "/study_extent_state_inset.png"), plot = p3, width = 5, height = 5)

##all combined ----
p.se <- ggpubr::ggarrange(p3, p1, p2, nrow = 1, labels = "AUTO", font.label = list(size = 21))
ggsave(paste0(fig.save.path, "/study_extent_full_plot.png"), plot = p.se, width = 14.5, height = 9, bg = "white")

 } #end section

## looking at whether upper river bss is lower effort than fjs ----
io.all <- lapply(c("Alewife", "American Shad", "Blueback Herring", "Striped Bass"), load_spp_dat)
io.all <- rbindlist(io.all)

io.tmp <- io.all[, unique(.SD), .SDcols = c("uniq_id", "survey", "riv_mile")]
io.tmp <- io.tmp[, .(ttl_hauls = .N), by = c("survey", "riv_mile")]
io.tmp[, in_or_out := fifelse(riv_mile >= 77 & riv_mile <= 108, "in", "out")]
io.tmp[, .(mn_hauls_per_survey_riverwide = mean(ttl_hauls)), by = "survey"]
io.tmp[, .(mn_hauls_per_survey_section = mean(ttl_hauls)), by = c("survey", "in_or_out")]
# haul.sum <- io.pos[, .(ttl_hauls = .N), by = c("survey", "assigned_riv_km")]
# 
# mn.haul <- haul.sum[, .(mn_hauls = sum(ttl_hauls)/245), by = "survey"]
# sect.mn.haul <- haul.sum[assigned_riv_km >= 125 & assigned_riv_km <= 175, .(section_mn_hauls = sum(ttl_hauls)/70), by = "survey"]
# 
# sect.mn.haul[survey == "fjs", section_mn_hauls] / mn.haul[survey == "fjs", mn_hauls]
# sect.mn.haul[survey == "bss", section_mn_hauls] / mn.haul[survey == "bss", mn_hauls]
# 
# ggplot(data = haul.sum, aes(x = assigned_riv_km, y = ttl_hauls)) +
#   geom_col() +
#   scale_x_continuous(n.breaks = 10) +
#   facet_wrap(~ survey)


# figure 0: cpue fjs v. bss ----
if(fig_0) {
  
spp.summ.all <- fread(paste0(PATH, "/inshore_offshore/multi-species_summarized_cpue.csv"))
spp.summ.all[, mn_cpue := mn_cpue * 1000]
spp.summ.all[survey == "fjs", survey := "offshore"]
spp.summ.all[survey == "bss", survey := "inshore"]
spp.summ.wide <- dcast(data = spp.summ.all, common_name + year ~ survey, value.var = "mn_cpue")

options(scipen = 999)
p1 <- ggplot(data = spp.summ.all, aes(x = year, y = mn_cpue, fill = survey, color = survey)) +
  geom_line(show.legend = FALSE) +
  geom_point(size = 3, shape = 21, color = "black") +
  scale_y_continuous(transform = "log10", n.breaks = 6, name = bquote(mean~annual~CPUE~"*"~10^~3~(catch~"/"~m^2))) +
  facet_wrap(~ common_name, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = survey.pal) +
  scale_color_manual(values = survey.pal) +
  scale_x_continuous(n.breaks = 5) +
  io_theme +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        axis.title = element_text(size = 14),
        axis.title.y = element_text(margin = margin(r = 5)),
        legend.text = element_text(size = 14, hjust = 0, margin = margin(r = 7)),
        legend.title = element_text(size = 14, hjust = 1),
        legend.key.size = unit(1, 'cm'),
        legend.position = "bottom")

p2 <- ggplot(data = spp.summ.wide, aes(x = inshore, y = offshore)) +
  geom_smooth(method = "lm", formula = y ~ x, color = "black") +
  geom_point(fill = "#00798c", alpha = 0.8, shape = 21, size = 3, color = "black") +
  # scale_fill_viridis_c(breaks = seq(1980, 2015, 10), labels = c("'80", "'90", "'00", "'10")) +
  facet_wrap(~ common_name, scales = "free", ncol = 1) +
  scale_x_continuous(transform = "log10", n.breaks = 6, name = bquote(inshore~CPUE~"*"~10^~3)) +
  scale_y_continuous(transform = "log10", n.breaks = 6, name = bquote(offshore~CPUE~"*"~10^~3)) +
  io_theme +
  theme(axis.title.x = element_text(margin = margin(t = 5)),
        axis.title.y = element_text(margin = margin(r = 5)),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14, vjust = 0.9, hjust = 1),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "none")

p3 <- ggpubr::ggarrange(p1, p2, nrow = 1, common.legend = T, widths = c(1, 0.6), align = "v", labels = c("A", "B"), legend = "bottom")

p3

ggsave(paste0(fig.save.path, "/raw_cpue_survey_comparison_combined.png"), plot = p3, width = 9, height = 9, bg = "white")
# ggsave(paste0(fig.save.path, "/raw_cpue_survey_comparison_timeseries.png"), plot = p1, width = 10, height = 6)
# ggsave(paste0(fig.save.path, "/raw_cpue_survey_comparison_1v1.png"), plot = p2, width = 10, height = 6)

}

#figure 1: pred v. obs ----
if(fig_1) {
  po.spp.l <- fread(paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_density_estimates_", mod.type, ".csv"))
  po.spp.l[survey_set == "all", survey_set := "inshore+offshore"]
  po.spp.l[survey_set == "bss", survey_set := "inshore"]
  po.spp.l[survey_set == "fjs", survey_set := "offshore"]
  po.spp.l[, survey_set := factor(survey_set, levels = c("inshore", "offshore", "inshore+offshore"))]
  p <- ggplot(data = po.spp.l[obs_catch > 0], aes(x = obs_catch/(area_swept), y = pred_dens, color = survey_set)) +
    geom_point(alpha = 0.3) +
    scale_color_manual(values = survey.pal) +
    geom_smooth(method = "glm", formula = y ~ x, color = "black") +
    scale_x_continuous(transform = "log10") +
    scale_y_continuous(transform = "log10") +
    facet_wrap(~ species + survey_set, ncol = 3, scales = "free") +
    labs(x = "observed density", y = "predicted density") +
    # facet_grid(rows = vars(survey_set), cols = vars(species), scales = "free") +
    io_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          plot.margin = unit(c(1,1,1,1), "cm"),
          legend.position = "none")
  p
  ggsave(paste0(fig.save.path, "/obs_v_pred_density_nosumm_", mod.type, ".png"), plot = p, width = 12, height = 10)
  
  p <- ggplot(data = po.spp.l[obs_catch == 0]) +
    geom_histogram(aes(x = pred_dens+1, fill = survey_set), color = "black", linewidth = 0.1) +
    scale_fill_manual(values = survey.pal) +
    facet_wrap(~ species + survey_set, ncol = 3, scales = "free") +
    scale_x_continuous(transform = "log1p", n.breaks = 4) +
    # scale_y_continuous(transform = "log1p", n.breaks = 3) +
    labs(x = "predicted density at observed catch = 0", y = "frequency") +
    io_theme +
    theme(legend.position = "none")
  p
  ggsave(paste0(fig.save.path, "/obs_v_pred_density_zero_catch_hist_", mod.type, ".png"), plot = p, width = 10, height = 12)
  
  tmp <- po.spp.l[, .(obs_catch, area_swept, pred_dens, time, species, survey_set)]
  tmp[, `:=` (obs_dens = obs_catch / (area_swept))] 
  tmp <- tmp[, lapply(.SD, mean, na.rm = TRUE), by = .(time, species, survey_set), .SDcols = c("obs_catch", "obs_dens", "pred_dens")]
  
  options(scipen = 99)
  p <- ggplot(data = tmp, aes(x = obs_dens, y = pred_dens, fill = survey_set, color = survey_set)) +
    geom_point(alpha = 0.5, size = 3, color = "black", shape = 21) +
    geom_smooth(method = "lm", formula = y~x, show.legend = F) +
    scale_fill_manual(values = survey.pal, name = "survey used") +
    scale_color_manual(values = survey.pal, name = "survey used") +
    scale_x_continuous(transform = "log10", name = expression(paste(bold("observed"), " mean annual density"))) +
    scale_y_continuous(transform = "log10", name = expression(paste(bold("predicted"), " mean annual density"))) +
    facet_wrap(~ species, scales = "free") +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 5))) +
    io_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "bottom",
          plot.margin = margin(5, 20, 5, 5))
  p
  ggsave(paste0(fig.save.path, "/obs_v_pred_density_annual_summ_", mod.type, ".png"), plot = p, width = 8, height = 7)
}

#figure 2: index comparison ----
if(fig_2) {
  
  ai.io <- fread(paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_abundance_estimates_", mod.type, ".csv"))
  ai.io[survey_set == "all", survey_set := "inshore+offshore"]
  ai.io[survey_set == "bss", survey_set := "inshore"]
  ai.io[survey_set == "fjs", survey_set := "offshore"]
  ai.dec <- fread(paste0(PATH, "/inshore_offshore/multi-species_abundance_estimates_vDEC.csv"))
  
  #modify dec to fit with io
  tmp <- ai.dec[year <= max(ai.io$time) & year >= min(ai.io$time), .(estimate = gm, survey_set = "DEC BSS geomean", time = year, species)]
  
  ai.io <- rbindlist(list(ai.io, tmp), fill = TRUE)
  ai.io[survey_set == "DEC BSS geomean", survey_set := "DEC inshore geomean"]
  
  ai.io[, survey_set := factor(survey_set, levels = c("DEC inshore geomean", "inshore", "offshore", "inshore+offshore"))]
  
  # options(scipen = 99)
  # p1 <- ggplot(data = ai.io[survey_set != "DEC BSS geomean"], aes(x = time, y = estimate, color = survey_set)) +
  #   # geom_ribbon(aes(ymin = lower, ymax = upper, colour = survey_set)) +
  #   geom_path() +
  #   geom_point() +
  #   facet_wrap(~ species, scales = "free_y") +
  #   scale_color_manual(values = survey.pal, name = "index type") +
  #   scale_y_continuous(transform = "log10", name = bquote(abundance~index~(count~"/"~m^2)~on~log[10]~scale)) +
  #   scale_x_continuous(n.breaks = 10, name = "year") +
  #   io_theme +
  #   theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  # p1
  # ggsave(paste0(fig.save.path, "/abundance_index_survey_comparison_nogeomean_", mod.type, ".png"), plot = p1, width = 12, height = 6)
  
  p2 <- ggplot(data = ai.io, aes(x = time, y = estimate, color = survey_set, fill = survey_set)) +
    # geom_ribbon(aes(ymin = lower, ymax = upper, colour = survey_set)) +
    geom_path(show.legend = F) +
    geom_point(shape = 21, alpha = 0.7, color = "black", size = 2) +
    facet_wrap(~ species, scales = "free_y", ncol = 2) +
    scale_color_manual(values = survey.pal, name = "index type") +
    scale_fill_manual(values = survey.pal, name = "index type") +
    scale_y_continuous(transform = "log10", name = "abundance index (count or geomean)") +
    scale_x_continuous(n.breaks = 5, name = "year") +
    io_theme +
    guides(color = guide_legend(override.aes = list(size = 6))) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5), 
          legend.title = element_text(hjust = 0, size = 14), legend.position = "bottom")
  p2
  ggsave(paste0(fig.save.path, "/abundance_index_survey_comparison_wgeomean_", mod.type, ".png"), plot = p2, width = 9, height = 6)
  
  ai.io[, mn_estimate := mean(estimate), by = c("species", "survey_set")]
  ai.io[, norm_estimate := estimate / mn_estimate]

  ai.wide <- dcast(ai.io, species + time ~ survey_set, value.var = "norm_estimate")
  
  #get r2 and slope values
  lm.res.dt <- data.table()
  for(i in seq_along(spp.l)) {
    
    spp <- spp.l[[i]]
    
    for(survey.y in c("inshore+offshore", 
                      "DEC inshore geomean", 
                      "inshore"
                      )) {
      for(survey.x in c("inshore", 
                        "offshore", 
                        "DEC inshore geomean")) {
        
        if(survey.x == survey.y) { next }
        dat <- ai.wide[species == spp]
        
        f <- paste0("log1p(`", survey.y, "`) ~ log1p(`", survey.x, "`)")
        # f <- paste0("`", survey.y, "` ~ `", survey.x, "`")
        
        l.res <- lm(as.formula(f), data = dat)
        ci <- confint(l.res)
        
        s.cor <- cor.test(log1p(dat[[survey.x]]), log1p(dat[[survey.y]]), method = "spearman", exact = F)
        p.cor <- cor.test(log1p(dat[[survey.x]]), log1p(dat[[survey.y]]), method = "pearson")
        
        tmp <- data.table(species = spp,
                          comparison = f,
                          intercept = coef(l.res)[1],
                          intercept_lowerci = ci["(Intercept)", "2.5 %"],
                          intercept_upperci = ci["(Intercept)", "97.5 %"],
                          slope = coef(l.res)[2],
                          slope_lowerci = ci[2, "2.5 %"],
                          slope_upperci = ci[2, "97.5 %"],
                          adj_r_sq = summary(l.res)$adj.r.squared,
                          lm_p_val = summary(l.res)$coef[2, 4],
                          sp_corr = s.cor$estimate,
                          sp_corr_p_val = s.cor$p.value,
                          p_corr = p.cor$estimate,
                          p_corr_p_val = p.cor$p.value)
        
        lm.res.dt <- rbindlist(list(lm.res.dt, tmp))
        
      }
    }
  }
  
  lm.res.dt[, comparison := gsub("`", "", comparison)]
  lm.res.dt[, (names(lm.res.dt)[sapply(lm.res.dt, is.numeric)]) := lapply(.SD, round, 3), .SDcols = is.numeric]
  
  #remove repeat
  lm.res.dt <- lm.res.dt[comparison != "log1p(inshore) ~ log1p(DEC inshore geomean)"]
  lm.res.dt[, `:=` (intercept_formatted = paste0(intercept, " (", intercept_lowerci, ", ", intercept_upperci, ")"),
                    slope_formatted = paste0(slope, " (", slope_lowerci, ", ", slope_upperci, ")"))]
  lm.res.dt[, `:=` (mn_corr = mean(sp_corr),
                    min_corr = min(sp_corr),
                    max_corr = max(sp_corr)), by = c("species")]
  
  fwrite(lm.res.dt, paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_index_comparison_lm_results_", mod.type, ".csv"))
  
  ##I can't think of a better way to do this right now :O
  lm.res.dt[grep("DEC", comparison), .(mn_rho = mean(sp_corr))]
  lm.res.dt[grep("inshore)", comparison), .(mn_rho = mean(sp_corr))]
  lm.res.dt[grep("offshore)", comparison), .(mn_rho = mean(sp_corr))]
  lm.res.dt[grep("inshore\\+offshore", comparison), .(mn_rho = mean(sp_corr))]
  
  ai_1v1_theme <- list(
    # scale_fill_viridis_c(breaks = seq(1980, 2015, 10)),
    scale_fill_manual(values = survey.pal, name = "index type"),
    scale_color_manual(values = survey.pal, name = "index type"),
    geom_point(shape = 21, size = 3, color = "black", alpha = 0.4),
    geom_smooth(aes(color = variable), method = "lm", se = F, alpha = 0.2, formula = y~x, show.legend = FALSE, linewidth = 1),
    # scale_y_continuous(transform = "log10"),
    # scale_x_continuous(transform = "log10"),
    facet_wrap(~ species, scales = "free", ncol = 1),
    io_theme,
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  )
  
  # #make dec v. offshore and inshore
  # ai.d.f.b <- melt(ai.wide, id.vars = c("species", "time", "DEC inshore geomean"), measure.vars = c("offshore", 
  #                                                                                                   "inshore", 
  #                                                                                                   "inshore+offshore"))
  # 
  # ggplot(data = ai.d.f.b, aes(x = `DEC inshore geomean`, y = value, fill = variable)) +
  #   ai_1v1_theme +
  #   scale_x_continuous(transform = "log1p") +
  #   scale_y_continuous(transform = "log1p") +
  #   labs(x = "normalized DEC inshore geomean", y = "normalized VAST abundance estimate")
  # ggsave(paste0(fig.save.path, "/abundance_index_survey_comp_dec_v_inshore-offshore_", mod.type, ".png"), width = 8, height = 6)
  
  #make offshore v. inshore
  # ggplot(data = ai.wide, aes(x = inshore, y = offshore, color = species, shape = species, fill = species)) +
  #   labs(x = "VAST inshore", y = "VAST offshore") +
  #   geom_smooth(method = "lm", alpha = 0.2, formula = y~x, linetype = "dashed", show.legend = FALSE) +
  #   geom_point(size = 3, color = "black") +
  #   scale_fill_manual(values = spp.pal) +
  #   scale_color_manual(values = spp.pal) +
  #   scale_shape_manual(values = spp.shape.pal) +
  #   scale_x_continuous(transform = "log1p") +
  #   scale_y_continuous(transform = "log1p") +
  #   facet_wrap(~ species, scales = "free") +
  #   io_theme +
  #   theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none") 
  # ggsave(paste0(fig.save.path, "/abundance_index_survey_comp_offshore_v_inshore_", mod.type, ".png"), width = 6, height = 6)
  
  #make combined v. offshore and inshore
  ai.c.f.b <- melt(ai.wide, id.vars = c("species", "time", "inshore+offshore"), measure.vars = c("offshore", 
                                                                                                 "inshore", 
                                                                                                 "DEC inshore geomean"))
  
  p4 <- ggplot(data = ai.c.f.b, aes(x = `inshore+offshore`, y = value, fill = variable)) +
    ai_1v1_theme +
    scale_x_continuous(transform = "log1p") +
    scale_y_continuous(transform = "log1p") +
    labs(x = "inshore + offshore\nnormalized index", y = "single survey normalized index")
  # p4
  # ggsave(paste0(fig.save.path, "/abundance_index_survey_comp_combined_v_inshore-offshore_", mod.type, ".png"), width = 8, height = 6)

  p3 <- ggplot(data = ai.io, aes(x = time, y = norm_estimate, color = survey_set, group = survey_set)) +
    # geom_ribbon(aes(ymin = lower, ymax = upper, colour = survey_set)) +
    # geom_errorbar(aes(ymin = norm_estimate - norm_std_error, ymax = norm_estimate + norm_std_error)) +
    geom_path(show.legend = F, linewidth = 1) +
    # geom_smooth(method = "loess") + ##doesn't really show the trends well enough
    geom_point(size = 2) +
    facet_wrap(~ species, scales = "free_y", ncol = 1) +
    scale_color_manual(values = survey.pal, name = "index type") +
    scale_y_continuous(name = "normalized abundance index") +
    scale_x_continuous(n.breaks = 5, name = "year") +
    guides(color = guide_legend(override.aes = list(size = 6))) +
    io_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.key.size = unit(1, 'cm'),
          legend.position = "none",
          legend.title = element_text(hjust = 1, size = 12),
          legend.text = element_text(hjust = 0, margin = margin(r = 7)))
  # p3
  # ggsave(paste0(fig.save.path, "/norm_abundance_index_survey_comparison_wgeomean_", mod.type, ".png"), plot = p3, width = 6, height = 8)
  
  ggpubr::ggarrange(p3, p4, common.legend = TRUE, legend = "bottom", widths = c(1, 0.5), labels = c("A", "B"), align = "h", nrow = 1)
  ggsave(paste0(fig.save.path, "/norm_abundance_index_comparison_combined_time_1v1_", mod.type, ".png"), width = 9, height = 9, bg = "white")
  # 
  # p4 <- ggplot(data = ai.io[survey_set != "DEC inshore geomean"], aes(x = time, y = norm_estimate, color = survey_set)) +
  #   # geom_ribbon(aes(ymin = lower, ymax = upper, colour = survey_set)) +
  #   geom_path() +
  #   geom_point() +
  #   facet_wrap(~ species, scales = "free_y") +
  #   scale_color_manual(values = survey.pal, name = "index type") +
  #   scale_y_continuous(name = "normalized abundance index") +
  #   scale_x_continuous(n.breaks = 10, name = "year") +
  #   io_theme +
  #   theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  # ggsave(paste0(fig.save.path, "/norm_abundance_index_survey_comparison_nogeomean.png"), plot = p4, width = 12, height = 6)
}


# figure 3: spatial density comparisons ----
##this became complicated, so moving to its own script (99_figures_spatial_density.R)

# figure 4: cv comparison ----
if(fig_4) {
  
  cv.est <- fread(paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_simple_uncertainty_estimates_", mod.type, ".csv"))
  cv.est[survey_set == "all", survey_set := "inshore+offshore"]
  cv.est[survey_set == "bss", survey_set := "inshore"]
  cv.est[survey_set == "fjs", survey_set := "offshore"]
  cv.est[, survey_set := factor(survey_set, levels = c("inshore", "offshore", "inshore+offshore"))]
  cv.est[, label := paste0(survey_set, "\n", type)]
  
  cv.summ <- cv.est[type != "CPUE", .(mn_cv = mean(cv_yr)), by = c("species", "survey_set", "type")]
  cv.summ <- cv.summ[, `:=` (type = fifelse(type == "DEC geomean", "DEC\ngeomean", type),
                             cv_outlier = fifelse(mn_cv > 1, TRUE, FALSE))]
  cv.summ <- cv.summ[order(survey_set)]
  # cv.summ[, label := factor(label, levels = label[match(c("DEC geomean", "CPUE", "VAST"), type)])]
  
  cv_theme <- list(
    facet_grid(~ type, scales = "free_x", space = "free_x", switch = "both"),
    scale_fill_manual(values = spp.pal, name = ""),
    scale_color_manual(values = spp.pal, name = ""),
    scale_shape_manual(values = spp.shape.pal, name = ""),
    labs(x = "abundance index estimate type & data source"),
    # facet_wrap(~ type, scales = "free_x"),
    theme_classic(),
    theme(axis.text = io_theme$axis.text,
          strip.text = io_theme$strip.text,
          axis.title.y = io_theme$axis.title,
          legend.text = io_theme$legend.text,
          strip.placement = "outside",
          # strip.background = element_blank(),
          # legend.title = io_theme$legend.text,
          panel.grid.major.y = element_line(color = "grey"))
  )
  
  ggplot(data = cv.summ) +
    geom_path(aes(x = survey_set, y = mn_cv, color = species, group = species), show.legend = F) +
    geom_point(aes(x = survey_set, y = mn_cv, fill = species, shape = species), size = 4, alpha = 0.9) +
    cv_theme +
    scale_y_continuous(expand = c(0,0), n.breaks = 10, name = "mean annual CV (SE / abundance estimate)") +
    coord_cartesian(ylim = c(0, 0.75))
  
  ggsave(paste0(fig.save.path, "/cv_survey_comparison_lessthan1_path_", mod.type, ".png"), width = 8, height = 6)
  
  ggplot(data = cv.summ) +
    geom_path(aes(x = survey_set, y = mn_cv, color = species, group = species), show.legend = F) +
    geom_point(aes(x = survey_set, y = mn_cv, fill = species, shape = species), size = 4, alpha = 0.9) +
    # scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5), guide = "none") +
    cv_theme +
    scale_y_continuous(transform = "log10", n.breaks = 10, 
                       name = bquote(mean~annual~CV~(SE~"/"~abundance~estimate)))
  
  ggsave(paste0(fig.save.path, "/cv_survey_comparison_all_path_", mod.type, ".png"), width = 8, height = 6)
  
  ##annual cv est
  ggplot(data = cv.est[!type %in% c("CPUE", "DEC geomean")], aes(x = year, y = cv_yr)) +
    geom_path(aes(color = survey_set), show.legend = FALSE) +
    geom_point(aes(fill = survey_set), shape = 21, size = 3, alpha = 0.8) +
    scale_fill_manual(values = survey.pal, name = "index type") +
    scale_color_manual(values = survey.pal) +
    scale_y_continuous(transform = "log10", name = bquote(CV~(SE~"/"~abundance))) +
    facet_wrap(~ species, scales = "free_y") +
    io_theme
  ggsave(paste0(fig.save.path, "/cv_annual_survey_comparison_all_path_", mod.type, ".png"), width = 8, height = 6)
  
  ggplot(data = cv.est[!type %in% c("CPUE")], aes(x = survey_set, y = cv_yr)) +
    # geom_violin(aes(fill = species)) +
    geom_boxplot(aes(fill = species), color = "black", outliers = FALSE) +
    geom_point(aes(fill = species), shape = 21, position = position_jitterdodge(jitter.width = 0.06, dodge.width = 0.75),
               show.legend = FALSE, alpha = 0.6) +
    scale_fill_manual(values = spp.pal) +
    scale_color_manual(values = spp.pal) +
    scale_y_continuous(transform = "log10", n.breaks = 10, 
                       name = bquote(annual~CV~(SE~"/"~abundance))) +
    cv_theme
  ggsave(paste0(fig.save.path, "/cv_survey_comparison_all_boxplot_", mod.type, ".png"), width = 8, height = 6)
  
  ggplot(data = cv.est[!type %in% c("CPUE", "DEC geomean")], aes(x = species, y = cv_yr)) +
    # geom_violin(aes(fill = species)) +
    geom_boxplot(aes(fill = survey_set), color = "black", outliers = FALSE) +
    geom_point(aes(fill = survey_set), shape = 21, position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
               show.legend = FALSE, alpha = 0.6, size = 2) +
    scale_fill_manual(values = survey.pal, name = "index type") +
    # scale_color_manual(values = survey.pal) +
    scale_y_continuous(transform = "log10", n.breaks = 10, 
                       name = bquote(annual~CV~(SE~"/"~abundance~estimate))) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "bottom",
          panel.grid.major.y = element_line(color = "lightgray"),
          panel.grid.major.x = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          panel.background = element_rect(fill = "white"),
          strip.background = element_rect(fill = "white"),
          axis.text = element_text(size = 12, color = "black"),
          strip.text = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12, color = "black")) +
    facet_wrap(~ species, scales = "free_x", nrow = 1, strip.position = "bottom") +
    labs(x = "")
  ggsave(paste0(fig.save.path, "/cv_survey_comparison_all_boxplot_byspecies_", mod.type, ".png"), width = 10, height = 5)
  
}

#figure s1: rmse over knots ----
if(fig_s1) {
  
  mod.res <- fread(paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_knot_comp_results.csv"))
  mod.res.top <- fread(paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_knot_comp_results_filtered.csv"))
  
  cv.est <- fread(paste0(PATH, "/inshore_offshore/", model.run.location, "param_comp_knot_abundance_indices.csv"))
  cv.est[, cv_yr := std_error_for_estimate / estimate]
  cv.summ <- cv.est[, .(mn_cv = mean(cv_yr)), by = c("species", "survey_type", "knots")]
  
  mod.res <- cv.summ[mod.res, on = c("species", "survey_type" = "survey_set", "knots")]
  mod.res[survey_type == "all", survey_type := "inshore+offshore"]
  mod.res[survey_type == "bss", survey_type := "inshore"]
  mod.res[survey_type == "fjs", survey_type := "offshore"]
  mod.res[, selected_mod := fifelse(model_name %in% mod.res.top$model_name, TRUE, FALSE)]
  mod.res <- mod.res[order(knots)]

  ggplot(data = mod.res[!is.na(aic)], aes(x = knots, y = rmse)) +
    # geom_smooth(method = "loess", se = F,) +
    geom_path(color = "blue", linewidth = 1.1) +
    geom_point(aes(fill = selected_mod, size = mn_cv), shape = 21) +
    scale_fill_manual(name = "selected model?", values = c("TRUE" = "red", "FALSE" = "gray")) +
    scale_size_continuous(name = "mean index CV", range = c(1,5), breaks = c(seq(0, 1, 0.2), 2, 4)) +
    facet_wrap(species ~ survey_type, scales = "free_y", nrow = 4) +
    theme_classic() +
    theme(axis.text = element_text(size = 12, color = "black"))
  ggsave(paste0(fig.save.path, "/knot_rmse_comparison.png"), width = 10, height = 6)
  
}

#pseudocoordinate example ----
if(psc_ex) {
  
  psc <- fread(file.path(PATH, "inshore_offshore", "pseudo_coordinate_options.csv"))
  hr.poly <- read_sf(paste0(PATH, "/env_data/gis_layers/NOAA_rivermile_polygon_edited.shp"))
  
  rmile = 100
  
  psc.sf <- psc[riv_mile == rmile]
  psc.sf <- st_as_sf(psc.sf, coords = c("X", "Y"), crs = 26918)
  
  ggplot() +
    geom_sf(data = hr.poly[hr.poly$Rmile == rmile,]) +
    geom_sf(data = psc.sf, aes(color = pt_location), size = 3) +
    scale_color_discrete(name = "pseudo-coord location") +
    scale_y_continuous(n.breaks = 4) +
    io_theme
  ggsave(file.path(fig.save.path, "pseudo_coord_example_rmile100.png"), width = 7, height = 5)
}

#dharma diagnostics ----
if(dharma) {
  library(VAST)
  
  if(!dir.exists(paste0(PATH, "/inshore_offshore/", model.run.location, "dharma_plots"))) {
    dir.create(paste0(PATH, "/inshore_offshore/", model.run.location, "dharma_plots"))
  }
  
  mod.top <- fread(paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_", gsub(" ", "_", model.run.type), "_results_filtered.csv"))
  
  for(s in spp.l) {
    
    m <- mod.top[species == s]
    
    spp.wd <- paste0(PATH, "/inshore_offshore/", model.run.location, gsub(" ", "_", s))
    
    for(i in 1:nrow(m)) {
      
      if(file.exists(paste0(PATH, "/inshore_offshore/", model.run.location, "dharma_plots/", 
                            gsub(".rds", ".png", m[i, model_name])))) {
        next
      }
      
      vast.mod <- readRDS(paste0(spp.wd, "/vast_models/", m[i, model_name]))
      vast.mod <- reload_model(vast.mod)
      
      vm.res <- summary(vast.mod, what = "residuals")
      
      png(file=paste0(PATH, "/inshore_offshore/", model.run.location, "dharma_plots/", 
                      gsub(".rds", ".png", m[i, model_name])), 
          width=8, height=4, res=200, units='in')
      plot_dharma(vm.res)
      dev.off()
    }
    
  }
  
}


# looking at relative catch x time of day ----
spp.l <- c("Alewife", "American Shad", "Blueback Herring", "Striped Bass")

p_list <- vector("list", 4)
names(p_list) <- spp.l
for(spp in spp.l) {
  spp.wd = paste0(PATH, "/inshore_offshore/", gsub(" ", "_", spp))
  io.dat.spp <- fread(paste0(spp.wd, "/", gsub(" ", "_", spp), "_vast_prep_dat_", survey2run, ".csv"))
  
  io.dat.spp[, noon_diff_cat := as.factor(round(solar_noon_diff))]
  io.summ <- io.dat.spp[, .(ttl_catch = sum(ct_yoy), ttl_haul = .N), by = c("survey", "noon_diff_cat")]
  io.summ[, ttl_cpue := ttl_catch / ttl_haul]
  io.summ[, `:=` (min_ttl_cpue = min(ttl_cpue), max_ttl_cpue = max(ttl_cpue)), by = survey]
  io.summ[, norm_cpue := (ttl_cpue - min_ttl_cpue) / (max_ttl_cpue - min_ttl_cpue)]
  
  p_list[[spp]] <- ggplot(data = io.summ, aes(x = noon_diff_cat, y = norm_cpue)) +
    geom_col() +
    facet_wrap(~ survey, ncol = 1)
}

# looking at SB north v south catch ----
io.dat.spp <- fread(paste0(PATH, "/inshore_offshore/Striped_Bass/Striped_Bass_vast_prep_dat_all.csv"))
io.dat.spp[, riv_sect := fifelse(riv_mile > 39, "north", "south")]

io.summ <- io.dat.spp[ct_yoy != 0, .(ttl_catch = sum(ct_yoy), ttl_haul = .N), by = c("survey", "riv_sect")]
io.summ[, ttl_cpue := ttl_catch / ttl_haul]

# looking at alewife v. bbh spatial catch distribution ----
aw <- fread(paste0(PATH, "/inshore_offshore/Alewife/Alewife_vast_prep_dat_all.csv"))
bbh <- fread(paste0(PATH, "/inshore_offshore/Blueback_Herring/Blueback_Herring_vast_prep_dat_all.csv"))

ssp.dat <- rbindlist(list(aw, bbh))
ssp.dat[, `:=` (long_grp = round(long, 3), lat_grp = round(lat, 2))]
ssp.dat <- ssp.dat[, .(sum_ct = sum(ct_yoy)), by = c("common_name", "survey", "year", "long_grp", "lat_grp")]

ggplot(data = ssp.dat[common_name == "Alewife" & sum_ct > 0], aes(x = long_grp, y = lat_grp, color = log1p(sum_ct))) +
  geom_point() +
  scale_color_viridis_c() +
  facet_wrap(~ survey)
ggplot(data = ssp.dat[common_name == "Blueback Herring" & sum_ct > 0], aes(x = long_grp, y = lat_grp, color = log1p(sum_ct))) +
  geom_point() +
  scale_color_viridis_c() +
  facet_wrap(~ survey)

gg.col.p <- list(geom_col(aes(x = lat_grp, y = sum_ct)),
                 scale_y_continuous(trans = "log1p"),
                 theme(axis.text.y = element_blank()),
                 facet_wrap(~ survey, ncol = 1, scales = "free_y"))

ggplot(data = ssp.dat[common_name == "Alewife" & sum_ct > 0], ) +
  gg.col.p
ggplot(data = ssp.dat[common_name == "Blueback Herring" & sum_ct > 0]) +
  gg.col.p


wide <- dcast(ssp.dat, common_name + year + lat_grp ~ survey, fun.aggregate = sum, value.var = "sum_ct")
wide <- wide[!is.na(bss) & !is.na(fjs)]
ggplot(data = wide, aes(x = bss, y = fjs, color = factor(lat_grp))) +
  geom_point() +
  geom_smooth(method = "lm", se  = F) +
  scale_y_continuous(trans = "log1p") +
  scale_x_continuous(trans = "log1p") +
  facet_wrap(~ common_name, scales = "free") +
  theme(legend.position = "none")

#looking at comparison of sturgeon catch  ----
fjs <- fread(paste0(PATH, "/bio_data/fall_juvenile_survey_combined.csv"))
fjs <- fjs[year <= 2014 & year >= 1980] 
bss <- fread(paste0(PATH, "/bio_data/beach_seine_plus_combined.csv"))
bss <- bss[year <= 2014 & year >= 1980] 

fjs.st <- fjs[grep("sturgeon", species_name, ignore.case = T)][ct_total != 0]
bss.st <- bss[grep("sturgeon", species_name, ignore.case = T)][ct_total != 0]

fjs.st[, .(sum = sum(ct_total)), by = species_name]
bss.st[, .(sum = sum(ct_total)), by = species_name]

#some summary tables ----

mod.res <- fread(paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_", mod.type, "_results_filtered.csv"))

comb.only <- mod.res[survey_set == "all"]
comb.only <- comb.only[, .(species, catch_cov, dens_cov, max_grad, rmse)]
fwrite(comb.only, paste0(PATH, "/inshore_offshore/multi-species_combined_model_results_table.csv"))

mod.save <- mod.res[, .(species, survey_set, catch_cov, dens_cov, max_grad, rmse)]
fwrite(mod.save, paste0(PATH, "/inshore_offshore/multi-species_all_final_model_results_table.csv"))


# predictor shape plots ----
if(pred_splines_check) {
  
  preds <- c("sam_dpth", "riv_dpth", "solar_noon_diff")
  
  for(spp in spp.l) {
    
    sp.dat <- load_spp_dat(spp = spp)
    
    for(p in preds) {
      
      p1 <- ggplot(data = sp.dat[!is.na(get(p))], aes(x = log(.data[[p]]+1), y = ct_yoy)) +
        geom_point() +
        facet_wrap(~ surv_label, scales = "free") +
        theme_classic()
      ggsave(paste0(spp.wd, "/figures/pred_relationship_", p, ".png"), width = 7, height = 4)
      
    }
    
  }
  
  #looking at splines ----
  # library(splines)
  # 
  # set.seed(123)
  # x <- io.dat.spp.all$riv_dpth
  # 
  # # Simulate y: bell/hill shape with curved tails
  # y <- io.dat.spp.all$ct_yoy
  # 
  # data <- data.frame(x = x, y = y)
  # x_seq <- seq(min(x), max(x), length.out = 200)
  # 
  # # Fit spline models
  # model_ns <- lm(y ~ ns(log(x+1), df = 5), data = data)
  # model_bs <- lm(y ~ bs(log(x+1), df = 5), data = data)
  # 
  # # Plot
  # plot(x, y, main = "Comparison of bs() vs ns()", xlab = "x", ylab = "y", pch = 16, col = "gray")
  # lines(x_seq, predict(model_ns, newdata = data.frame(x = x_seq)), col = "blue", lwd = 2)
  # lines(x_seq, predict(model_bs, newdata = data.frame(x = x_seq)), col = "red", lwd = 2)
  # legend("topright", legend = c("ns() (Natural Spline)", "bs() (B-spline)"),
  #        col = c("blue", "red"), lty = 1, lwd = 2)
  # 
  # plot(resid(model_bs), main = "Residuals")
  # 
  # plot(fitted(model_ns), resid(model_ns), xlab = "Fitted", ylab = "Residuals")
  # abline(h = 0, col = "red")
  
  plot(density(resid(model_ns)))
  hist(resid(model_ns), breaks = 20)
  
  # Fit a spline model
  fit <- lm(ct_yoy ~ bs(solar_noon_diff, df = 4), data = io.dat.spp[ct_yoy > 0])
  summary(fit)
  # Make a sequence of x values over the observed range
  newdat <- data.frame(solar_noon_diff = seq(
    min(io.dat.spp$solar_noon_diff, na.rm = TRUE),
    max(io.dat.spp$solar_noon_diff, na.rm = TRUE),
    length.out = 200
  ))
  
  # Get predictions
  pred <- predict(fit, newdata = newdat, interval = "confidence")
  newdat$pred <- pred[, "fit"]
  newdat$lwr <- pred[, "lwr"]
  newdat$upr <- pred[, "upr"]
  
  # Plot data + fitted smooth
  ggplot() +
    geom_point(data = io.dat.spp[ct_yoy > 0], aes(x = solar_noon_diff, y = ct_yoy), pch = 16, color = "grey") +
    geom_line(data = newdat, aes(x = solar_noon_diff, y = pred), color = "blue") +
    geom_ribbon(data = newdat, aes(ymin = lwr, ymax = upr, x = solar_noon_diff), alpha = 0.2, fill = "blue") +
    scale_y_continuous(transform = "log10")
  
}


