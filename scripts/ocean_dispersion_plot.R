args <- commandArgs(trailingOnly = TRUE)
particle_csv <- args[[1]]
out_pdf <- args[[2]]

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

p <- readr::read_csv(particle_csv, show_col_types = FALSE)
ord <- p |>
  dplyr::group_by(source_station) |>
  dplyr::summarise(median_disp = median(final_distance_from_station_km, na.rm = TRUE), .groups = 'drop') |>
  dplyr::arrange(median_disp) |>
  dplyr::pull(source_station)

p <- p |>
  dplyr::mutate(source_station = factor(source_station, levels = ord))

p1 <- ggplot(p, aes(final_distance_from_station_km, source_station)) +
  geom_violin(fill = '#9ecae1', color = '#225ea8', alpha = 0.8, scale = 'width') +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = 'white', color = 'black') +
  labs(x = 'Final distance from source station (km)', y = 'Source station') +
  theme_bw()

p2 <- ggplot(p, aes(source_station, final_distance_from_station_km)) +
  geom_boxplot(fill = '#fdd0a2', color = '#d94801', outlier.alpha = 0.25) +
  coord_flip() +
  labs(x = 'Source station', y = 'Final distance from source station (km)') +
  theme_bw()

pdf(out_pdf, width = 11.69, height = 8.27, useDingbats = FALSE)
print(p1 + p2 + plot_layout(widths = c(1.2, 1)))
dev.off()
