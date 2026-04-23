source('/home/userbio/Sea4Blue/scripts/utils/helpers.R')

run_oceanography_plots <- function(
  scenario_dir = '/home/userbio/Sea4Blue/results/oceanography_sensitivity/sens_60d_75km_5',
  stations_csv = '/home/userbio/Sea4Blue/parcels/stations.csv',
  prefix = 'sens60d75km5'
) {
  load_packages(c('dplyr', 'tidyr', 'ggplot2', 'readr', 'tibble', 'scales'))
  zone_colors <- sea4blue_zone_colors()

  meta <- load_phyloseq_atlantic(sea4blue_config())$sample_data |>
    tibble::rownames_to_column('station') |>
    dplyr::select(station, Zona)

  st <- readr::read_csv(stations_csv, show_col_types = FALSE) |>
    dplyr::transmute(
      station = .data[['Nominativo.campione.']],
      sample_date = .data[['data.x']],
      lat = .data[['Latitude.x']],
      lon = .data[['Longitude.x']]
    ) |>
    dplyr::left_join(meta, by = 'station') |>
    dplyr::mutate(Zona = factor(as.character(Zona), levels = names(zone_colors)))

  conn <- as.matrix(read.csv(file.path(scenario_dir, 'backward_connectivity_matrix.csv'), row.names = 1, check.names = FALSE))
  mins <- as.matrix(read.csv(file.path(scenario_dir, 'backward_min_travel_days.csv'), row.names = 1, check.names = FALSE))
  ret <- readr::read_csv(file.path(scenario_dir, 'retention_by_station.csv'), show_col_types = FALSE)
  acc <- readr::read_csv(file.path(scenario_dir, 'station_accumulation_scores.csv'), show_col_types = FALSE)

  heat_df <- as.data.frame(as.table(conn)) |>
    dplyr::rename(source = Var1, sink = Var2, connectivity = Freq) |>
    dplyr::mutate(
      source = factor(source, levels = rownames(conn)),
      sink = factor(sink, levels = colnames(conn))
    )

  idx <- which(row(conn) != col(conn) & conn > 0, arr.ind = TRUE)
  edge_df <- tibble::tibble(
    source = rownames(conn)[idx[, 1]],
    sink = colnames(conn)[idx[, 2]],
    connectivity = conn[idx],
    min_days = mins[idx]
  ) |>
    dplyr::left_join(st |> dplyr::select(station, source_lon = lon, source_lat = lat, source_zone = Zona), by = c('source' = 'station')) |>
    dplyr::left_join(st |> dplyr::select(station, sink_lon = lon, sink_lat = lat, sink_zone = Zona), by = c('sink' = 'station')) |>
    dplyr::arrange(dplyr::desc(connectivity), min_days)

  station_metrics <- st |>
    dplyr::left_join(ret, by = 'station') |>
    dplyr::left_join(acc, by = 'station') |>
    dplyr::mutate(station = factor(station, levels = station[order(lon)])) |>
    tidyr::pivot_longer(cols = c(retention, accumulation, incoming_strength), names_to = 'metric', values_to = 'value')

  heat_plot <- ggplot(heat_df, aes(sink, source, fill = connectivity)) +
    geom_tile(color = 'white', linewidth = 0.2) +
    scale_fill_gradientn(colors = c('#f7fbff', '#9ecae1', '#3182bd', '#08519c'), limits = c(0, 1)) +
    labs(x = 'Sink station', y = 'Source station', fill = 'Connectivity') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  map_plot <- ggplot() +
    geom_segment(
      data = edge_df,
      aes(x = source_lon, y = source_lat, xend = sink_lon, yend = sink_lat, linewidth = connectivity, alpha = connectivity),
      arrow = grid::arrow(length = grid::unit(0.12, 'inches')),
      color = '#1f4e79'
    ) +
    geom_point(data = st, aes(lon, lat, fill = Zona), shape = 21, color = 'black', size = 3.2, stroke = 0.25) +
    geom_text(data = st, aes(lon, lat, label = station), nudge_y = 0.35, size = 3) +
    scale_fill_manual(values = zone_colors, drop = FALSE) +
    scale_linewidth(range = c(0.3, 1.6), guide = guide_legend(order = 2)) +
    scale_alpha(range = c(0.35, 0.95), guide = 'none') +
    labs(x = 'Longitude', y = 'Latitude', fill = 'Zone', linewidth = 'Connectivity') +
    theme_bw()

  metrics_plot <- ggplot(station_metrics, aes(station, value, fill = Zona)) +
    geom_col(color = 'black', linewidth = 0.2) +
    facet_wrap(~ metric, scales = 'free_y', ncol = 1) +
    scale_fill_manual(values = zone_colors, drop = FALSE) +
    labs(x = 'Station', y = 'Value', fill = 'Zone') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  out_fig <- '/home/userbio/Sea4Blue/results/figures'
  out_tab <- '/home/userbio/Sea4Blue/results/tables'
  ggsave(file.path(out_fig, paste0('ocean_connectivity_heatmap_', prefix, '.pdf')), heat_plot, width = 11.69, height = 8.27, units = 'in')
  ggsave(file.path(out_fig, paste0('ocean_connectivity_map_', prefix, '.pdf')), map_plot, width = 11.69, height = 8.27, units = 'in')
  ggsave(file.path(out_fig, paste0('ocean_station_metrics_', prefix, '.pdf')), metrics_plot, width = 8.27, height = 11.69, units = 'in')
  write.csv(edge_df, file.path(out_tab, paste0('ocean_connectivity_edges_', prefix, '.csv')), row.names = FALSE)

  list(edges = edge_df, heat_plot = heat_plot, map_plot = map_plot, metrics_plot = metrics_plot)
}
