### INSHORE OFFSHORE - MAKE EXTRAPOLATION GRID + PSEUDO COORDINATES###

# set up
library(data.table); library(sf); library(terra)#; library(ggplot2)

PATH <- getwd()
if(grepl("scripts", PATH)) { PATH <- gsub("/scripts", "", PATH) }

#load in river shapefile
##this is river km shapefile from R. Pendleton
# rkm <- read_sf(paste0(PATH, "/env_data/gis_layers/HudsonRiverKms_Poly.shp"))
# rkm <- subset(rkm, select = -c(OBJECTID_1, Join_Count, TARGET_FID, OBJECTID))

##river mile shapefile from R. Pendleton
rmi <- read_sf(paste0(PATH, "/env_data/gis_layers/NOAA_rivermile_polygon_edited.shp"))
# rmi <- subset(rmi, select = -c(OBJECTID, ID_, PERIMETER, ID, Shape_Leng, Shape_Area, POLY_AREA))
rmi <- subset(rmi, select = (Rmile))
# rmi <- rmi[rmi$Rmile != 0,]

##check it worked:
# ggplot() + geom_sf(data = rmi[rmi$Rmile != 0,], aes(fill = Rmile))


#EXTRAP GRID ----
#save shapefile crs
# st_crs(rmi) #already UTM! can get epsg code with this
prj.crs <- 26918 #project UTM (zone 18 for ny)
ll.crs <- 4326 #wgs 84 // for converting to lat long

#make the grid (smaller cell size = higher resolution = longer time to completion)
cell_size <- 250 #this is used later, so set variable here **METERS BETWEEN POINTS*
region_grid <- st_as_sf(st_make_grid(rmi, cellsize = cell_size, what = "centers"))

##check it worked:
# ggplot() + geom_sf(data = rmi, aes(fill = Rmile)) + geom_sf(data = region_grid)
# ggplot() +
#   geom_point(data = hudson_river_VAST_extrap_grid_250cs, aes(x = Lon, y = Lat), size = 0.001, shape = 3) +
#   # geom_point(data = hudson_river_VAST_extrap_grid_500cs, aes(x = Lon, y = Lat), size = 0.001, shape = 4) +
#   coord_equal()

#determine which grid points fall within the river shapefile and where
region_grid.rmi <- st_join(region_grid, rmi)

##check it worked:
# ggplot() + geom_sf(data = rmi) + geom_sf(data = tmp, aes(color = rmi))

#add depth to each grid point
bath <- rast(paste0(PATH, "/env_data/gis_layers/30m_grid/dem/estuary.dem"))

bath.out <- extract(bath, vect(region_grid.rmi), method = "simple")
region_grid.rmi$depth <- abs(bath.out$estuary)

#convert back to lon/lat coordinates and make the grid a data.frame for use in VAST
region_grid.rmi.ll <- st_transform(region_grid.rmi, crs = ll.crs)

region_grid.df <- cbind(st_drop_geometry(region_grid.rmi.ll), st_coordinates(region_grid.rmi.ll)) |>
  as.data.frame()

#prep grid dataframe
region_grid.df$Area_km2 <- (cell_size/1000)^2
names(region_grid.df)[grep("X|Y", names(region_grid.df))] <- c("Lon", "Lat")
region_grid.df$row <- 1:nrow(region_grid.df)

#remove grid that does not overlap the river polygon
##note here that one example removed grid points that were above or below a certain depth etc. so can set some other spatial restrictions it seems like
region_grid.df <- subset(region_grid.df, !is.na(Rmile))

##deal with coordinates outside the data layer (NA depth)
DT <- copy(region_grid.df)
DT.na <- DT[is.na(DT$depth),]
DT.notna <- DT[!is.na(DT$depth),]

DT.sf <- st_as_sf(DT.na, coords = c("Lon", "Lat"), crs = 4269) |>
  st_transform(st_crs(bath))

#get nearest cells (previously converted to points)
##for fixing coordinates outside the depth layer
bath.pts <- as.points(bath, na.rm = TRUE)
bath.pts <- st_as_sf(bath.pts)
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
DT[is.na(depth), depth := riv_dpth_est]
DT[, riv_dpth_est := NULL]

region_grid.df <- copy(DT)

#make strata column
region_grid.df$STRATA <- ifelse(region_grid.df$depth < 3.048, "inshore", "offshore") # shore as defined by HRBMP final report

##check it worked:
tmp <- st_as_sf(region_grid.df, coords = c("Lon", "Lat"), crs = ll.crs)
p <- ggplot() +
  geom_sf(data = rmi) +
  # geom_sf(data = region_grid.rmi.ll, aes(color = rmi))
  geom_sf(data = tmp, aes(color = STRATA)) +
  scale_x_continuous(n.breaks = 2) +
  theme_bw()
p1 <- ggplot() +
  geom_sf(data = rmi[rmi$Rmile %in% c(24, 25, 26),]) +
  geom_sf(data = tmp[tmp$Rmile %in% c(24, 25, 26),], aes(color = STRATA)) +
  theme_bw()

ggpubr::ggarrange(p, p1)
ggsave(paste0(PATH, "/inshore_offshore/extrap_grid_example.png"), bg = "white", width = 10, height = 10)

#save for use in VAST
saveRDS(region_grid.df, file = paste0(PATH, "/inshore_offshore/hudson_river_VAST_extrap_grid_", cell_size, "cs.rds"))
fwrite(region_grid.df, paste0(PATH, "/inshore_offshore/hudson_river_VAST_extrap_grid_", cell_size, "cs.csv"))

#check if VAST recognizes the grid:
out <- make_extrapolation_info(Region = "user", strata.limits = data.frame(STRATA = "offshore"), input_grid = region_grid.df)
# out <- make_extrapolation_info(Region = "eastern_bering_sea")

# PSEUDOCOORDINATES ----
##read in NHD for getting river direction:
nhd <- read_sf(paste0(PATH, "/env_data/gis_layers/hudson_watershed_stream_clip.shp"))
nhd <- nhd[nhd$GNIS_NAME == "Hudson River" & !is.na(nhd$GNIS_NAME), ] |>
  st_transform(st_crs(rmi)) |>
  st_zm()

##first make potential coordinates along the length of the river in the center of each strata (rmi or riv mile)
get_mid_coords <- function(sf.obj) {
  
  # get nhd stream in river mile
  nhd.rmi <- suppressWarnings(st_intersection(nhd, sf.obj)) |>
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
  
  # get center of line
  # strata.center <- st_centroid(st_geometry(sf.obj))
  nhd.mid <- st_line_sample(nhd.rmi, sample = 0.5) |>
    st_coordinates()
  nhd.mid <- nhd.mid[, c("X", "Y")]
  
  # get perpendicular line (90 d rotation) to nhd
  delta.x <- nhd.c[2, 1] - nhd.c[1, 1]
  delta.y <- nhd.c[2, 2] - nhd.c[1, 2]
  
  p.90 <- rbind(
    nhd.mid + 150000 * c(-delta.y, delta.x) / sqrt(delta.x^2 + delta.y^2),
    nhd.mid - 150000 * c(-delta.y, delta.x) / sqrt(delta.x^2 + delta.y^2)
    )
  
  ##get line from points, clip to river
  l.90 <- st_sfc(st_linestring(p.90), crs = st_crs(sf.obj))
  l.90 <- st_intersection(l.90, sf.obj)
  
  if(length(l.90) == 0) { stop("No intersection found for lat. midpoint line for river mile", sf.obj$Rmile) }  # if no intersection, send error
  
  # extract coordinates and make them sf points
  coords <- st_coordinates(l.90)
  
  left <- coords[which.min(coords[,1]), 1:2]
  right <- coords[which.max(coords[,1]), 1:2]
  
  # now find midway between center point and edge points
  mid.left <- c((nhd.mid[1] + left[1]) / 2, (nhd.mid[2] + left[2]) / 2)
  mid.right <- c((nhd.mid[1] + right[1]) / 2, (nhd.mid[2] + right[2]) / 2)
  
  # convert to sf obj
  strata.coords <- as.data.frame(do.call(rbind, list(left, mid.left, nhd.mid, mid.right, right)))
  strata.coords$pt_location <- c("edge", "mid", "center", "mid", "edge")
  strata.coords$riv_mile <- sf.obj$Rmile
  
  strata.coords <- st_as_sf(strata.coords, coords = c("X", "Y"), crs = st_crs(sf.obj))
  
  # check for points not in the river
  out.test <- st_join(strata.coords, sf.obj)
  idx <- which(is.na(out.test$Rmile) & strata.coords$pt_location != "edge")
  
  if(length(idx) > 0) {

    for(i in idx) {
      
      tmp <- st_nearest_points(strata.coords[i,], l.90) |>
        st_cast("POINT")
      tmp <- st_as_sf(tmp)
      
      tmp <- st_join(tmp, sf.obj)
      
      if(all(is.na(tmp$Rmile))) {
        
        remove <- which(apply(st_coordinates(tmp), 1, function(row) all(row == st_coordinates(strata.coords[i,]))))
        tmp <- tmp[-remove,]

        tmp$Rmile <- sf.obj$Rmile
        
      } else {
        
        tmp <- tmp[!is.na(tmp$Rmile), ]
        
        
      }
      
      if(nrow(tmp) > 1 | nrow(tmp) == 0) { stop("Something went wrong finding a new point in river for river mile ", sf.obj$Rmile) }
      
      tmp$pt_location <- strata.coords[i, ]$pt_location
      names(tmp)[names(tmp) == "Rmile"] <- "riv_mile"
      st_geometry(tmp) <- "geometry"
  
      strata.coords <- do.call(rbind, list(strata.coords[-i, ], tmp))
      strata.coords <- strata.coords[order(rownames(strata.coords)), ] ##have to rearrange it will mess up index pull
      
    }
  }
  
  #check w/ plot
  # ggplot() +
  #   geom_sf(data = sf.obj) +
  #   geom_sf(data = nhd.rmi) +
  #   # geom_sf(data = st_as_sf(as.data.frame(st_coordinates(nhd.rmi)), coords = c("X", "Y"), crs = st_crs(sf.obj)), color = "red") +
  #   geom_sf(data = l.90) +
  #   geom_sf(data = strata.coords, aes(color = pt_location)) +
  #   geom_sf(data = tmp, shape = 4, color = "red")

  return(strata.coords)
  
}

mid.coords <- lapply(1:nrow(rmi), function(i) get_mid_coords(rmi[i, ]))
mid.coords <- do.call(rbind, mid.coords)

##save
mid.coords$X <- st_coordinates(mid.coords)[,"X"]
mid.coords$Y <- st_coordinates(mid.coords)[,"Y"]
write.csv(st_drop_geometry(mid.coords), paste0(PATH, "/inshore_offshore/pseudo_coordinate_options.csv"))

# Plot to check
# ggplot() + geom_sf(data = rmi[rmi$Rmile == 20,]) + geom_sf(data = mid.coords[mid.coords$riv_mile == 20, ], aes(color = pt_location))





