source('/home/userbio/Sea4Blue/scripts/utils/helpers.R')

classify_mechanism_boot <- function(bnti, rc_bc) {
  dplyr::case_when(
    bnti <= -2 ~ 'Homogeneous Selection',
    bnti >= 2 ~ 'Heterogeneous Selection',
    abs(bnti) < 2 & rc_bc >= 0.95 ~ 'Dispersal Limitation',
    abs(bnti) < 2 & rc_bc <= -0.95 ~ 'Homogenising Dispersal',
    TRUE ~ 'Drift'
  )
}

bootstrap_props <- function(df, group_cols = NULL, n_boot = 5000, seed = 42) {
  load_packages(c('dplyr', 'tidyr', 'tibble'))
  set.seed(seed)
  group_cols <- group_cols %||% character(0)
  mechs <- sort(unique(df$Mechanism))

  summarise_once <- function(dat) {
    if (length(group_cols) == 0) {
      dat |>
        dplyr::count(Mechanism, name = 'n') |>
        dplyr::mutate(prop = n / sum(n)) |>
        dplyr::select(Mechanism, prop)
    } else {
      dat |>
        dplyr::count(dplyr::across(dplyr::all_of(group_cols)), Mechanism, name = 'n') |>
        dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
        dplyr::mutate(prop = n / sum(n)) |>
        dplyr::ungroup() |>
        dplyr::select(dplyr::all_of(group_cols), Mechanism, prop)
    }
  }

  observed <- summarise_once(df) |>
    dplyr::rename(observed_prop = prop)

  boot_list <- vector('list', n_boot)
  n <- nrow(df)
  for (i in seq_len(n_boot)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    boot_list[[i]] <- summarise_once(df[idx, , drop = FALSE]) |>
      dplyr::mutate(iter = i)
  }
  boot_df <- dplyr::bind_rows(boot_list)

  if (length(group_cols) == 0) {
    ci <- boot_df |>
      dplyr::group_by(Mechanism) |>
      dplyr::summarise(
        ci_low = stats::quantile(prop, 0.025, na.rm = TRUE),
        ci_mid = stats::quantile(prop, 0.5, na.rm = TRUE),
        ci_high = stats::quantile(prop, 0.975, na.rm = TRUE),
        .groups = 'drop'
      )
  } else {
    ci <- boot_df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols)), Mechanism) |>
      dplyr::summarise(
        ci_low = stats::quantile(prop, 0.025, na.rm = TRUE),
        ci_mid = stats::quantile(prop, 0.5, na.rm = TRUE),
        ci_high = stats::quantile(prop, 0.975, na.rm = TRUE),
        .groups = 'drop'
      )
  }

  observed |>
    dplyr::full_join(ci, by = c(group_cols, 'Mechanism'))
}

run_mechanisms_bootstrap <- function(cfg = sea4blue_config(), n_boot = 5000, seed = 42) {
  load_packages(c('dplyr', 'ggplot2', 'reshape2', 'magrittr', 'geosphere', 'phyloseq', 'tibble'))
  zone_colors <- sea4blue_zone_colors()
  objs <- load_phyloseq_atlantic(cfg)
  meta <- objs$sample_data
  meta$Sample_ID <- rownames(meta)
  meta$thetao <- as.numeric(meta$thetao)

  bNTI <- as.matrix(read.csv(cfg$paths$btnti_csv, row.names = 1, check.names = FALSE))
  RC_BC <- as.matrix(read.csv(cfg$paths$rc_bray_csv, row.names = 1, check.names = FALSE))

  bNTI_mod <- magrittr::set_colnames(bNTI, rownames(bNTI))
  bNTI_mod[lower.tri(bNTI_mod)] <- NA
  diag(bNTI_mod) <- NA
  RC_BC_mod <- RC_BC
  rownames(RC_BC_mod) <- rownames(bNTI_mod)
  colnames(RC_BC_mod) <- colnames(bNTI_mod)
  RC_BC_mod[lower.tri(RC_BC_mod)] <- NA
  diag(RC_BC_mod) <- NA

  pair_df <- reshape2::melt(bNTI_mod) |>
    cbind(reshape2::melt(RC_BC_mod)) |>
    magrittr::set_colnames(c('From_Sample', 'To_Sample', 'bNTI', 'x1', 'x2', 'RC_BC')) |>
    dplyr::select(From_Sample, To_Sample, bNTI, RC_BC) |>
    dplyr::filter(!is.na(bNTI)) |>
    dplyr::mutate(Mechanism = classify_mechanism_boot(bNTI, RC_BC)) |>
    dplyr::left_join(meta |>
      dplyr::select(Sample_ID, Longitude.y, Latitude.y, Zona, thetao) |>
      dplyr::rename(From_Lon = Longitude.y, From_Lat = Latitude.y, Zone_1 = Zona, Temp_1 = thetao),
      by = c('From_Sample' = 'Sample_ID')) |>
    dplyr::left_join(meta |>
      dplyr::select(Sample_ID, Longitude.y, Latitude.y, Zona, thetao) |>
      dplyr::rename(To_Lon = Longitude.y, To_Lat = Latitude.y, Zone_2 = Zona, Temp_2 = thetao),
      by = c('To_Sample' = 'Sample_ID')) |>
    dplyr::mutate(
      geographic_km = geosphere::distHaversine(cbind(From_Lon, From_Lat), cbind(To_Lon, To_Lat)) / 1000,
      distance_bin = cut(geographic_km, breaks = c(0, 250, 500, 1000, 2000, 3000, 5000, Inf), include.lowest = TRUE, right = FALSE),
      temp_diff = abs(Temp_1 - Temp_2),
      temp_bin = cut(temp_diff, breaks = c(0, 1, 2, 4, 6, 8, Inf), include.lowest = TRUE, right = FALSE)
    )

  global_boot <- bootstrap_props(pair_df, group_cols = NULL, n_boot = n_boot, seed = seed)
  distance_boot <- bootstrap_props(pair_df |> dplyr::filter(!is.na(distance_bin)), group_cols = 'distance_bin', n_boot = n_boot, seed = seed)
  temperature_boot <- bootstrap_props(pair_df |> dplyr::filter(!is.na(temp_bin)), group_cols = 'temp_bin', n_boot = n_boot, seed = seed)

  global_plot <- ggplot(global_boot, aes(x = Mechanism, y = observed_prop * 100, fill = Mechanism)) +
    geom_col(color = 'black', linewidth = 0.2) +
    geom_errorbar(aes(ymin = ci_low * 100, ymax = ci_high * 100), width = 0.2) +
    coord_flip() +
    labs(x = NULL, y = 'Mechanism proportion (%)') +
    theme_bw() +
    theme(legend.position = 'none')

  distance_plot <- ggplot(distance_boot, aes(x = distance_bin, y = observed_prop * 100, fill = Mechanism)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.8, color = 'black', linewidth = 0.2) +
    geom_errorbar(aes(ymin = ci_low * 100, ymax = ci_high * 100), position = position_dodge(width = 0.85), width = 0.2) +
    labs(x = 'Geographic distance bin (km)', y = 'Mechanism proportion (%)') +
    theme_bw()

  temperature_plot <- ggplot(temperature_boot, aes(x = temp_bin, y = observed_prop * 100, fill = Mechanism)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.8, color = 'black', linewidth = 0.2) +
    geom_errorbar(aes(ymin = ci_low * 100, ymax = ci_high * 100), position = position_dodge(width = 0.85), width = 0.2) +
    labs(x = 'Temperature difference bin (deg C)', y = 'Mechanism proportion (%)') +
    theme_bw()

  list(
    pairwise = pair_df,
    global_boot = global_boot,
    distance_boot = distance_boot,
    temperature_boot = temperature_boot,
    plots = list(global = global_plot, distance = distance_plot, temperature = temperature_plot)
  )
}
