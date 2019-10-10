setwd("H:/Desktop/Spatial/Homework 1")

### Spatial Statistics: Homework 1 ### 


# Libraries needed 

library(maps)
library(maptools)
library(sf)
library(tmap)
library(tmaptools)
library(mapproj)
library(gstat)
library(tools)
library(dplyr)
library(stringr)
library(tidyr)
library(rgeos)
library(spdep)
library(ggplot2)

################################################################################ Question 1 ##########################
# Load and read in data 

colorado_dat_0 <- read.delim("colorados-t.dat")                              # Raw .dat file

latlong2county2 <- function(pointsDF) {                                       # Function to find state value 
  
  county <- map('county', fill=TRUE, col='transparent', plot=FALSE)
  
  IDs <- sapply(strsplit(county$names, ':'), function(x) x[1])
  
  county_sp <- map2SpatialPolygons(
    county, 
    IDs=IDs,
    proj4string=CRS('+proj=longlat +datum=WGS84')
  )
  
  pointsSP <- SpatialPoints(
    pointsDF,
    proj4string=CRS('+proj=longlat +datum=WGS84'))
  
  indices <- over(pointsSP, county_sp)
  
  countyNames <- sapply(county_sp@polygons, function(x) x@ID)
  
  countyNames[indices]
  
}

names(colorado_dat_0)

colorado_dat <- st_as_sf(
  x = colorado_dat_0,
  coords = c(
    "Lon",
    "Lat"
    ),
  crs = CRS(
    "+proj=utm +zone=13"
    )
  )

colorado_dat$county <- latlong2county2(data.frame(colorado_dat_0$Lon, colorado_dat_0$Lat))

colorado_dat[, 16:17] <- str_split_fixed(colorado_dat$county ,pattern = ",", n = 2)

names(colorado_dat)[16:17]<-c("state", "county")

colorado_dat[, 15] <- NULL

state_point <- data.frame(
  state = c("Wyoming", "Nebraska", "Kansas", "Oklahoma", "New Mexico", "Colorado", "Utah"),
  Lon = c(-105.861625, -102.086648, -101.310184, -102.072771, -108.434457,-102.624180, -109.562668),
  Lat = c(41.3, 41.256259, 39.0, 36.745762, 36.180285, 40.853806, 39.187828)
)

state_point <- st_as_sf(
  x = state_point,
  coords = c(
    "Lon",
    "Lat"
  ),
  crs = CRS(
    "+proj=utm +zone=13"
  )
)


colorado_dat <- colorado_dat %>%                                             # Separates county and state     
  mutate(state = toTitleCase(state), 
         county = toTitleCase(county)
         ) %>% 
  mutate(name = paste0(county, sep = ", ", state)
         )

for(i in colorado_dat$county){                                               # Replaces duplicate county names with NA
  colorado_dat$county[duplicated(colorado_dat$county)] <- NA 
}

table(colorado_dat$county)

# Create State Map 

states_map <- map(
  database = 'county',
  regions = c(
    "Colorado",
    "Wyoming",
    "New Mexico",
    "Utah",
    "Kansas",
    "Nebraska",
    "Oklahoma",
    "Texas",
    "Arizona"
  ), 
  fill = T,
  plot = F
)

state_fill <- map(
  database = 'state',
  regions = c(
    "Colorado",
    "Wyoming",
    "New Mexico",
    "Utah",
    "Kansas",
    "Nebraska",
    "Oklahoma",
    "Texas",
    "Arizona"
  ), 
  fill = T,
  plot = F
)

states_poly <- map2SpatialPolygons(
  states_map, 
  IDs = states_map$names, 
  proj4string = CRS(
    "+proj=utm +zone=13"
    )
)

states_fill_poly <- map2SpatialPolygons(
  state_fill, 
  IDs = state_fill$names, 
  proj4string = CRS(
    "+proj=utm +zone=13"
  )
)

max_lon <- max(colorado_dat_0$Lon)
min_lon <- min(colorado_dat_0$Lon)

max_lat <- max(colorado_dat_0$Lat)
min_lat <- min(colorado_dat_0$Lat)

names(states_fill_poly)

tm_shape(
  states_fill_poly,
  ylim=c(
    min_lat - .5, 
    max_lat + .5
  ),
  xlim=c(
    min_lon - .5, 
    max_lon + .5
  )
) + 
  tm_fill(
    col = "gray99"
  ) + 
  tm_borders(
    lwd = 3,
    col = "gray"
  )  + 
tm_shape(
  states_poly,
  ylim=c(
    min_lat - .5, 
    max_lat + .5
    ),
  xlim=c(
    min_lon - .5, 
    max_lon + .5
    )
) +
  tm_borders(
    lwd = 1.5,
    col = 'gray'
  ) + 
  tm_shape(
    colorado_dat
    ) + 
  tm_symbols(
    col = "Jan",
    palette = "-RdYlBu",
    n = 7,
    style = "jenks",
    border.lwd = 0.2,
    border.col = 'gray',
    alpha = 0.9,
    scale = 1.15,
    title.col = "10's of °C"
  ) + 
  tm_text(
    "county",
    textNA = "",
    size = .65,
    just = "bottom",
    shadow = T
  ) + 
  tm_shape(
    state_point
  ) +
  tm_text(
    "state",
    textNA = "",
    remove.overlap = F,
    shadow = T,
    fontface = "bold"
    
  ) + 
  tm_legend(
    position = c(
      "left", 
      "bottom"
      ),
    legend.outside = TRUE,
    frame = F,
    main.title = 'Temperatures in Colorado and Surrounding Counties, January'
  ) 

cutoff <- .5*max(
  dist(
    cbind(
      colorado_dat_0$Lon,
      colorado_dat_0$Lat
      )
    )
  )

bins <- 30

variogram <- variogram(
  Jan ~ 1, 
  locations = ~Lat + Lon,
  data = colorado_dat_0,
  cutoff = cutoff,
  width = cutoff/bins
)

variogram

plot(
  variogram,
  ylab = expression(
    paste(
      "Average", (0.5*(Y(s[i]) - Y(s[j])^2))
      )
    ),
  xlab = "Euclidean Distance (°)",
  cex = 2,
  pch = 20
  )

covariogram <- variogram(
  Jan ~ 1, 
  locations = ~Lat + Lon,
  data = colorado_dat_0,
  cutoff = cutoff,
  width = cutoff/bins,
  covariogram = TRUE
)

covariogram

plot(
  gamma ~ dist,
  data = covariogram,
  ylab = expression(
    paste(
      0.65*(Y(s[i]) - Y(s[j])^2)
    )
  ),
  xlab = "Euclidean Distance (°)",
  cex = 2,
  pch = 20
) 
abline(
  h = 0,
  lty = 2,
  col = 'black',
  lwd = 2
)

cutoff <- .65*max(
  dist(
    cbind(
      colorado_dat_0$Lon,
      colorado_dat_0$Lat
    )
  )
)

bins_co <- 10

covar_map <- variogram(
  Jan ~ 1, 
  locations = ~Lat + Lon,
  data = colorado_dat_0,
  cutoff = cutoff,
  width = cutoff/bins_co,
  covariogram = TRUE,
  map = TRUE
)

plot(covar_map)

dir_covar <- variogram(
  Jan ~ 1, 
  locations = ~Lat + Lon,
  data = colorado_dat_0,
  cutoff = cutoff,
  width = cutoff/bins,
  covariogram = TRUE,
  alpha = c(
    0,
    45,
    90,
    135
  )
)

plot(
  dir_covar,
  main = "Directional Covariograms",
  as.table = TRUE,
  cex = 1.5,
  pch = 16
)

######################################################################################### Question 2 #################

US_map <- st_read("acs_county_us.shp")

class(US_map)

US_map <- st_transform(
  x = US_map,
  crs = "+proj=laea +lat_0=45 +lon_0=-100"
)

class(US_map)

US_map$state<-sub(".*, ", "", US_map$NAME)

table(US_map$state)

texas <- US_map %>% 
  filter(
    state == "Texas"
    ) 

texas_gather <- gather(
  data = texas, 
  key = "year",
  value = "perc_un",
  un_2012:un_2016
)

texas_gather <- texas_gather %>% 
  mutate(
    year = ifelse(
      test = year == "un_2012", 
      yes = "2012",
      no = "2016"
      )
    )

uninsured_tx <- tm_shape(
  shp = texas_gather
  ) + 
  tm_polygons(
    col = "perc_un", 
    palette = 'PuRd',
    n = 6,
    style = "fisher",
    title = "% Uninsured",
    border.lwd = 0,
    border.alpha = 0.2
  ) + 
  tm_legend(
    position = c("right", "bottom"),
    frame = T,
    outside = T,
    main.title = "% Uninsured in Texas, 2012-2016"
  ) + 
  tm_facets(
    by = 'year'
  ) + 
  tm_layout(
    bg.color = "white",
    frame.double.line = T,
    outer.bg.color = "white"
  )

uninsured_tx

# Queen 

tx_map <- map(
  'county',
  'texas',
  fill = T,
  plot = F
)

tx_poly <- st_geometry(texas)

tx_cent <- map2SpatialPolygons(
  tx_map,
  IDs = tx_map$names,
  proj4string = CRS(
    "+proj=laea +lat_0=45 +lon_0=-100"
    )
  ) 
  
tx_queen_1 <- poly2nb(
  tx_cent,
  queen = TRUE
)

tx_queen_2 <- poly2nb(
  tx_poly,
  queen = TRUE
)

tx_queen

centroids <- gCentroid(
  tx_cent, 
  byid = TRUE
)

plot(
  tx_cent,
  col = 'palegreen',
  main = "Queen Based Adjacency of Counties in Texas"
)
plot(
  tx_queen_1,
  coordinates(
    centroids
    ),
  col = "blue",
  add = TRUE
  )

moran_2012 <- moran.test(
  x = texas$un_2012, 
  listw = nb2listw(
    neighbours = tx_queen_2, 
    style = "B"
    )
  )

moran_2012

moran_2016 <- moran.test(
  x = texas$un_2016, 
  listw = nb2listw(
    neighbours = tx_queen_2, 
    style = "B"
  )
)

moran_2016

comp_df <- t(
  data.frame(
    "Uninsured 2012" = moran_2012$estimate,
    "Uninsured 2016" = moran_2016$estimate
  )
)

tx_for_lag <- nb2mat(tx_queen_2, style = "B")

lagged_2016 <-  tx_for_lag %*% texas$un_2016

ggplot(
  data = texas, 
  aes(
    x = un_2016,
    y = lagged_2016
  )
) + 
  geom_point(
    
  ) + 
  geom_smooth(
    method = "lm"
  ) + 
  theme_classic(
    base_size = 13
  ) + 
  labs(
    x = "Uninsurance Rate",
    y = "Spatially-lagged Uninsurance Rate",
    title = "Spatially Lagged Insurance Rates; Texas, 2016"
  ) 

plot(
  sp.correlogram(
    tx_queen_2,
    texas$un_2016,
    order = 10,
    style = "B",
    method = "I"
    ),
  main = "Moran's I Correlogram for % Uninsured (2016)"
  )



# The end 



