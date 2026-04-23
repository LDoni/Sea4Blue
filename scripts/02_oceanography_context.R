source("/home/userbio/Sea4Blue/scripts/utils/helpers.R")

build_oceanography_context <- function(cfg = sea4blue_config()) {
  load_packages(c("phyloseq", "dplyr", "sf", "rnaturalearth", "raster", "ggplot2", "viridis", "RColorBrewer", "ncdf4", "oceanmap"))
  objs <- load_phyloseq_atlantic(cfg)
  ps <- objs$phyloseq_obj_css

  coord_df <- ps %>%
    subset_samples(!(sample.description %in% c("Laposnegative", "Lapospositive"))) %>%
    microbiome::meta() %>%
    dplyr::select(sample.description, Latitude, Longitude)
  colnames(coord_df) <- c("id", "lat", "lon")

  points_sf <- sf::st_as_sf(coord_df, coords = c("lon", "lat"), crs = 4326)
  longhurst <- sf::st_read(cfg$paths$longhurst_shp, quiet = TRUE)
  longhurst <- sf::st_make_valid(longhurst)
  longhurst <- sf::st_transform(longhurst, 4326)
  result <- sf::st_join(points_sf, longhurst["ProvCode"], join = sf::st_intersects)

  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  province_usate <- unique(result$ProvCode)
  longhurst_subset <- dplyr::filter(longhurst, ProvCode %in% province_usate)

  crop_vals_atl <- c(lat = c(-90, 10), lon = c(25, 45))
  sst_atl <- ncdf4::nc_open(cfg$paths$sst_nc) %>%
    oceanmap::nc2raster("sst") %>%
    raster::flip("y") %>%
    raster::crop(raster::extent(crop_vals_atl))
  sst_df <- as.data.frame(sst_atl, xy = TRUE)
  names(sst_df)[3] <- "sst"
  result$longhurst_provinces <- result$ProvCode

  sst_map <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = sst_df, ggplot2::aes(x = x, y = y, fill = sst)) +
    ggplot2::scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdBu")), name = "Annual SST (°C)", limits = c(0, 35), na.value = "transparent") +
    ggplot2::geom_sf(data = world, fill = "#D9D9D9", color = "gray70", linewidth = 0.3) +
    ggplot2::geom_sf(data = longhurst_subset, color = "gray60", fill = NA, linewidth = 0.4) +
    ggplot2::geom_sf(data = result, ggplot2::aes(color = longhurst_provinces), size = 3) +
    ggplot2::scale_color_viridis_d(name = "Longhurst Province") +
    ggplot2::geom_sf_text(data = result, ggplot2::aes(label = id), size = 3, nudge_y = 1, check_overlap = TRUE) +
    ggplot2::coord_sf(crs = 4326, xlim = c(-90, 10), ylim = c(25, 45), expand = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Sea Surface Temperature and Sample Points", subtitle = "SST from MODIS + Longhurst Provinces", x = "Longitude", y = "Latitude")

  list(points = result, longhurst = longhurst_subset, sst = sst_df, plot = sst_map)
}
