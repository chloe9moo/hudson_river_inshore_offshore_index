#make maps of density across the river
library(terra); library(tidyterra); library(RColorBrewer)

#set up ----
PATH <- getwd()

fig.save.path <- paste0(PATH, "/inshore_offshore/amarel_cluster/figures/")
#which run location
model.run.location = "amarel_cluster/" #""
#which type
mod.type = "knot_comp"
#target spp
spp.l = c("Alewife", "Striped Bass", "American Shad", "Blueback Herring")
spp = "Alewife"

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

#prep polys ----
##get full river extent as raster
hr.poly <- read_sf(paste0(PATH, "/env_data/gis_layers/HudsonRiverKms_Poly.shp"))

#river line (for directionality)
nhd <- read_sf(paste0(PATH, "/env_data/gis_layers/hudson_watershed_stream_clip.shp"))
nhd <- nhd[nhd$GNIS_NAME == "Hudson River" & !is.na(nhd$GNIS_NAME), ] |>
  st_transform(st_crs(hr.poly)) |>
  st_zm()

#prep density estimates ----
knot.dens <- fread(paste0(PATH, "/inshore_offshore/", model.run.location, "multi-species_knot_spatial_density_estimates_", mod.type, ".csv"))
knot.dens[, coord_id := .GRP, by = c("long", "lat")]

##get unique coordinates and pair with rkms
u.knot.dens <- knot.dens[, unique(.SD), .SDcols = c("coord_id", "lat", "long")]

kd.sf <- st_as_sf(u.knot.dens, coords = c("long", "lat"), crs = 4326) |>
  st_transform(st_crs(hr.poly))

#transform river to a rectangle and plot ----
get_rel_widths <- function(pt, crs.in = st_crs(hr.poly)) {
  #get relative location along width of river
  if(FALSE) {
    i <- 200
    pt <- st_geometry(kd.sf[i,])
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
      geom_sf(data = hr.poly[hr.poly$RKm %in% c(54, 55, 56, 57),]) +
      geom_sf(data = kd.sf[kd.sf$coord_id == 938,])
      # geom_sf(data = sub.hr) +
      # geom_sf(data = nhd.rmi, color = "blue") +
      # geom_sf(data = line.sf) +
      # geom_sf(data = pt.line) +
      # geom_sf(data = pt, shape = 4) +
      # geom_sf(data = tmp, shape = 4) +
      # geom_sf(data = transect) +
      # geom_sf(data = shore_pts)
  }
  
  if(any(names(sub.hr) == "Rmile")) {
    out <- data.table(assigned_riv_mile = targ.rm, river_width = dist_total, rel_position = rel_position)
  }
  if(any(names(sub.hr) == "RKm")) {
    out <- data.table(assigned_riv_km = targ.rm, river_width = dist_total, rel_position = rel_position)
  }
  
  return(out)
}

##looping to track progress..
kd.pos <- vector("list", length = nrow(kd.sf))
for(i in 1:nrow(kd.sf)) {
  kd.pos[[i]] <- get_rel_widths(st_geometry(kd.sf[i,]))
  
  cat(sprintf("\r%-20s", sprintf("%.0f%% complete...", 100 * i / nrow(kd.sf))))
  
  if(i == nrow(kd.sf)) { cat("\nDone!\n")}
}

kd.pos <- rbindlist(kd.pos)

kd.pos[which(round(kd.pos$rel_position, 2) > 1)]
kd.pos[is.na(kd.pos$rel_position)] ##these are occurring in the edge in a section that is not in the poly

kd.pos[is.na(rel_position) & assigned_riv_km == 56, rel_position := 1.0]
kd.pos[is.na(rel_position) & assigned_riv_km == 59, rel_position := 1.0]
kd.pos[is.na(rel_position) & assigned_riv_km == 75, rel_position := 0.0]

kd.sf$assigned_riv_km <- kd.pos$assigned_riv_km
kd.sf$river_width <- kd.pos$river_width
kd.sf$rel_position <- kd.pos$rel_position

##add to the density table
knot.dens <- knot.dens[st_drop_geometry(kd.sf), on = "coord_id"]

##'bin' the data
knot.dens[, rel_position := fifelse(rel_position > 1, 1, 
                                    fifelse(rel_position < 0, 0, rel_position))]
knot.dens[, width_bin := round(rel_position, digits = 1)]

##summarize
kd.summ.yr <- knot.dens[, lapply(.SD, mean, na.rm = T), 
                        by = c("survey_set", "species", "width_bin", "assigned_riv_km", "year"),
                        .SDcols = c("pred_dens", "r1_gct", "r2_gct")]

kd.summ <- knot.dens[, lapply(.SD, mean, na.rm = T), 
                     by = c("survey_set", "species", "width_bin", "assigned_riv_km"),
                     .SDcols = c("pred_dens", "r1_gct", "r2_gct")]

##plot ----
##this won't really work bc the fewer knots only have one point per rkm sometimes
ggplot(data = kd.summ[species == spp], aes(x = width_bin, y = assigned_riv_km, fill = pred_dens)) +
  geom_tile() +
  scale_y_continuous(n.breaks = 10, expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0),    
                     breaks = seq(0, 1, by = 0.25),
                     labels = c("0.0", "", "0.5", "", "1.0")) +
  scale_fill_viridis_c(trans = "log10",
                       # labels = scales::comma_format(),  # or use label = log10 if you prefer
                       name = "predicted density") +
  facet_wrap(~ survey_set) +
  labs(x = "relative position along width of the river\n", y = "river kilometer") +
  io_theme

#plot as raster ----
template <- rast(vect(hr.poly), res=2000) #2 km resolution?
raster.hr <- rasterize(vect(hr.poly), template, touches = T)

get_raster_knot_dens <- function(DT, 
                                 spp = NA, yrs2plot = NA, surv2plot = c("bss", "fjs", "all"), 
                                 col2pull = "pred_dens", norm = T, resid = F) {
  ##testing
  if(FALSE) {
    DT <- copy(knot.dens)
    spp = "Alewife"; yrs2plot = c(1990, 2000)
    surv2plot = c("bss", "fjs")
    norm = T
    resid = T
    col2pull = "pred_dens"
  }
    
  if(is.na(spp) | length(spp) > 1) { stop("set species") }
  if(all(is.na(yrs2plot))) { stop("set at least 1 year") }
  if(all(is.na(surv2plot)) | (length(surv2plot) > 1 & resid == FALSE)) { stop("set 1 survey") }
  if(resid == TRUE & length(surv2plot) != 2) { stop("pick 2 surveys to pull 'residuals' from") }
  
  DT.sub <- DT[species == spp & year %in% yrs2plot & survey_set %in% surv2plot]
  
  if(norm | resid) {
    DT.sub[, paste0("norm_", col2pull) := lapply(.SD, function(x) (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))),
           by = .(species, year, survey_set), .SDcols = col2pull]
    col2pull = paste0("norm_", col2pull)
  }
  
  #split into ind. dts to plot
  if(length(yrs2plot) > 1 | length(surv2plot) > 1) {
    DT.s.l <- split(DT.sub, list(DT.sub$year, DT.sub$survey_set))
  } else {
    DT.s.l <- copy(DT.sub)
  }
  
  #for each year, make spatial >> rasterize >> fill in gaps
  dgct.r.l <- lapply(DT.s.l, function(x, w.size = 7) {
    
    d.gct.sf <- st_as_sf(x, coords = c("long", "lat"), crs = 4326) |>
      st_transform(st_crs(hr.poly))
    
    ##points 2 rast
    dgct.rast <- rasterize(vect(d.gct.sf), template, field = col2pull, fun = mean, touches = T)
    
    ##moving window average to fill in NA cells, mask to rasterized river
    dgct.rast.m <- focal(dgct.rast, w = w.size, fun = mean, na.policy = "only")
    dgct.rast.m <- mask(dgct.rast.m, raster.hr, maskvalues = NA)
    
    return(dgct.rast.m)
    
  })
    
  dgct.r.l <- rast(dgct.r.l)
  
  if(resid) {
    
    res.all <- vector("list", length = length(yrs2plot))
    for(i in seq_along(yrs2plot)) {
      
      yr <- yrs2plot[i]
      
      ind1 <- grep(paste0("(?=.*", yr, ")(?=.*", surv2plot[1], ")"), names(dgct.r.l), perl = TRUE)
      ind2 <- grep(paste0("(?=.*", yr, ")(?=.*", surv2plot[2], ")"), names(dgct.r.l), perl = TRUE)
      
      res <- dgct.r.l[[ind1]] - dgct.r.l[[ind2]]
      names(res) <- paste0(yr, ".", surv2plot[1], "-", surv2plot[2])
      
      res.all[[i]] <- res
      
    }
    
    dgct.r.l <- rast(res.all)
    
  }
  
  return(dgct.r.l)
}


##make residual, and each survey set raster stack of years ----
#loop through species
for(spp in spp.l) {
  res.all <- get_raster_knot_dens(knot.dens, spp = spp, yrs2plot = unique(knot.dens$year), surv2plot = c("bss", "fjs"), resid = T)
  fjs.res <- get_raster_knot_dens(knot.dens, spp = spp, yrs2plot = unique(knot.dens$year), surv2plot = c("fjs"), norm = T)
  bss.res <- get_raster_knot_dens(knot.dens, spp = spp, yrs2plot = unique(knot.dens$year), surv2plot = c("bss"), norm = T)
  comb.res <- get_raster_knot_dens(knot.dens, spp = spp, yrs2plot = unique(knot.dens$year), surv2plot = c("all"), norm = T)
  
  list.results <- list(res.all, 
                       fjs.res, 
                       bss.res, 
                       comb.res)
  
  list.mean <- lapply(list.results, mean, na.rm = T)
  
  
  ##plot the four summarized rasters for each species + save ----
  options(scipen = 99)
  ext <- list.mean[[1]] |>
    project("EPSG:4326", mask = TRUE) |>
    ext() %>%
    as.vector()
  
  map.theme <- list(
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3)(ext[c("xmin", "xmax")])),
    theme(panel.background = element_rect(fill = "white"),
          panel.grid.major = element_line(color = "lightgray"),
          panel.border = element_rect(color = "black", fill = NA),
          axis.text.x = element_text(angle = 47, vjust = 1, hjust = 1),
          axis.text = element_text(size = 12, color = "black"),
          strip.text = element_text(size = 12, face = "bold"),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12, color = "black"),
          legend.title = element_text(size = 12, color = "black"))
  )
  
  p4 <- ggplot() +
    geom_spatraster(data = list.mean[[1]]) +
    scale_fill_gradientn(colors = brewer.pal(10, "RdBu"), na.value = NA, 
                         breaks = c(-1, -0.5, 0, 0.5, 1), limits = c(-1.01, 1.01),
                         labels = c("offshore greater", "", "equal density", "", "inshore greater"),
                         name = "predicted density\ncomparison\n") +
    map.theme +
    theme(panel.background = element_rect(fill = "#3b3b3b"),
          panel.grid.major = element_blank(),
          plot.margin = unit(c(1, 0.01, 1, 0.5), "cm"))
  p4
  
  tmp <- c()
  for(i in 2:4) {
    tmp <- c(tmp, global(list.mean[[i]], "max", na.rm = TRUE)[[1]])
  }
  tmp <- max(tmp)
  
  dens.map.theme <- list(
    map.theme,
    scale_fill_viridis_c(transform = "log10", na.value = "transparent", name = "mean\nnormalized\npred. dens.\n",
                         # breaks = c(0, 0.001, 0.003, 0.010, 0.030, 0.1, 0.3),
                         labels = scales::label_number(accuracy = 0.0001),
                         limits = c(1e-6, 1)),
    theme(plot.margin = unit(c(0.5, 0.001, 0.5, 0.001), "cm")),
    labs(x = "", y = "")
  )
  
  #make labels for river polygon
  hr.poly <- st_transform(hr.poly, st_crs(list.mean[[2]]))
  bbox <- st_bbox(hr.poly)
  xmin <- bbox[1] 
  xmax <- bbox[3]
  xbreaks <- seq(xmin, xmax, length.out = 3)
  xbreaks <- round(xbreaks, 2)
  hr.poly$r_km_label <- fifelse(hr.poly$RKm %in% seq(0, 225, 25), hr.poly$RKm, NA)
  hr.center <- st_centroid(hr.poly)
  hr.center$X <- st_coordinates(hr.center)[, 'X']
  hr.center$Y <- st_coordinates(hr.center)[, 'Y']
  hr.center <- hr.center[which(!is.na(hr.center$r_km_label)),]
  
  #fjs
  p3 <- ggplot() +
    geom_spatraster(data = list.mean[[2]]) +
    dens.map.theme +
    ggtitle("offshore")
  
  #bss
  p2 <- ggplot() +
    geom_spatraster(data = list.mean[[3]]) +
    dens.map.theme +
    ggtitle("inshore")
  
  #combined
  p1 <- ggplot() +
    ggrepel::geom_text_repel(data = hr.center, aes(x = X, y = Y, label = r_km_label), 
                             xlim = c(bbox[3]-0.05, bbox[3]+0.1), min.segment.length = 0.01, size = 4.2) +
    geom_spatraster(data = list.mean[[4]]) +
    dens.map.theme  +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3)(ext[c("xmin", "xmax")]), 
                       limits = c(bbox[1], bbox[3]+15000)) +
    ggtitle("inshore+offshore")
  
  p.surv <- ggpubr::ggarrange(p1, p2, p3, nrow = 1, common.legend = T, widths = c(1.2, 1, 1), legend = "right")
  # p.surv
  p.all <- ggpubr::ggarrange(p.surv, p4, nrow = 1, widths = c(0.66, 0.33))
  p.all
  ggsave(paste0(fig.save.path, "predicted_density_all_", spp, ".png"), plot = p.all, width = 11, height = 7, bg = "white")
}


#multi year plots ----
spp = "Striped Bass"
surv = "all"
bss.res <- get_raster_knot_dens(knot.dens, spp = spp, yrs2plot = unique(knot.dens$year), surv2plot = surv, norm = T)
names(bss.res) <- gsub(paste0(".", surv), "", names(bss.res))
ggplot() +
  geom_spatraster(data = bss.res) +
  facet_wrap(~ lyr, nrow = 2) +
  scale_fill_viridis_c(transform = "log10", na.value = "transparent", name = "normalized\npredicted density") +
  map.theme +
  theme(axis.text.x = element_blank())
