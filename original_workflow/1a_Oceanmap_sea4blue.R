#### Load dependencies ####

library(tidyverse)
library(raster)

source("Pacific_Bacterioplankton/R/Datalist_Wrangling_Functions.R")
source("Pacific_Bacterioplankton/R/Import_Data.R")

#scarico da https://oceandata.sci.gsfc.nasa.gov/file_search/ 
wget --auth-no-challenge --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies "https://oceandata.sci.gsfc.nasa.gov/ob/getfile/AQUA_MODIS.20220101_20221231.L3m.YR.SST.sst.4km.nc"

#### Read station and satellite data and formatting ####

Stations_Combined_Atlantic <- ps %>%
  subset_samples(!(sample.description %in% c("Laposnegative", "Lapospositive"))) %>%
  meta()


 






crop.vals_Atl <- c(
  lat = c(25, 45),     # 29.3007 - 36.99627, esteso a [29, 38]
  lon = c(-90, 10)     # -77.936917 - -7.874533, esteso a [-79, -7]
)



crop.vals_Atl <- c(lat = c(-90,10), 
                   lon = c(25, 45))

SST_Atl <- ncdf4::nc_open("AQUA_MODIS.20220101_20221231.L3m.YR.SST.sst.4km.nc") %>%
  oceanmap::nc2raster(., "sst") %>%
  raster::flip(., "y") %>%
  raster::crop(., raster::extent(crop.vals_Atl))

 

##### Plotting Oceanmap with Chl a backgroung ####

library(oceanmap)

 

oceanmap::v(SST_Atl, cbpos = "b", pal = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(300)),
            zlim = c(0,35), 
            cb.xlab = expression("Annual SST (°C)"),
            bwd = 0.01, grid = F, replace.na = F, border = "#504f4f",
            cex.ticks = 1, axeslabels = F, figdim = c(4,5), show.colorbar = T)

 







library(ggplot2)
library(raster)
library(sf)
library(dplyr)
library(viridis)

# 1. Converti la SST croppata in data frame per ggplot
SST_df <- as.data.frame(SST_Atl, xy = TRUE)
names(SST_df)[3] <- "sst"  # rinomina la colonna del valore

# 2. Prepara coord.df con geometria sf
points_sf <- st_as_sf(coord.df,
                      coords = c("Longitude ", "Latitude"),
                      crs = 4326)

# Aggiungi la colonna longhurst_provinces (se non c'è già)
points_sf$longhurst_provinces <- as.data.frame(result)$ProvCode



library("RColorBrewer")

ggplot() +
  # Sfondo SST
  geom_raster(data = SST_df, aes(x = x, y = y, fill = sst)) +
  scale_fill_gradientn(colors = rev(brewer.pal(11, "RdBu")),
                       name = "Annual SST (°C)",
                       limits = c(0, 35),
                       na.value = "transparent") +
  
  # Mappa dei continenti
  geom_sf(data = world, fill = "#D9D9D9", color = "gray70", size = 0.3) +
  
  # Sovrapposizione poligoni Longhurst senza fill
  geom_sf(data = longhurst_subset, color = "gray60", fill = NA, size = 0.4) +
  
  # Punti colorati per provincia
  geom_sf(data = points_sf, aes(color = longhurst_provinces), size = 3) +
  scale_color_viridis_d(name = "Longhurst Province") +
  
  # Etichette
  geom_sf_text(data = points_sf, aes(label = id), size = 3, nudge_y = 1, check_overlap = TRUE) +
  
  # Coordinate e stile
  coord_sf(crs = 4326, xlim = c(-90, 10), ylim = c(25, 45), expand = FALSE) +
  theme_minimal() +
  labs(title = "Sea Surface Temperature and Sample Points",
       subtitle = "SST from MODIS + Longhurst Provinces",
       x = "Longitude", y = "Latitude")





 
