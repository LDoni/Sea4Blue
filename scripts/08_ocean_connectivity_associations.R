source("/home/userbio/Sea4Blue/scripts/utils/helpers.R")

get_lower_tri_df <- function(mat, value_name) {
  idx <- which(lower.tri(mat), arr.ind = TRUE)
  out <- tibble::tibble(
    sample_1 = rownames(mat)[idx[, 1]],
    sample_2 = colnames(mat)[idx[, 2]],
    value = mat[idx]
  )
  names(out)[3] <- value_name
  out
}

run_ocean_connectivity_associations <- function(cfg = sea4blue_config()) {
  load_packages(c("phyloseq", "dplyr", "vegan", "geosphere", "tibble", "ggplot2"))
  objs <- load_phyloseq_atlantic(cfg)
  ps <- objs$phyloseq_obj_css
  meta <- data.frame(sample_data(ps))
  sample_ids <- rownames(meta)

  conn_path <- file.path(cfg$ocean_dir, "backward_connectivity_matrix.csv")
  if (!file.exists(conn_path)) stop("Missing connectivity matrix: ", conn_path)
  conn <- as.matrix(read.csv(conn_path, row.names = 1, check.names = FALSE))
  conn <- conn[sample_ids, sample_ids, drop = FALSE]
  sym_conn <- (conn + t(conn)) / 2
  diag(sym_conn) <- NA_real_

  bray <- as.matrix(phyloseq::distance(ps, method = "bray"))
  unifrac <- as.matrix(phyloseq::distance(ps, method = "unifrac", weighted = FALSE))
  coords <- as.matrix(meta[, c("Longitude.y", "Latitude.y")])
  rownames(coords) <- sample_ids
  geo <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
  rownames(geo) <- colnames(geo) <- sample_ids
  temp <- as.matrix(dist(as.numeric(meta$thetao)))
  rownames(temp) <- colnames(temp) <- sample_ids
  ocean_dist <- 1 - sym_conn

  pair_df <- get_lower_tri_df(bray, "bray") |>
    dplyr::left_join(get_lower_tri_df(unifrac, "unifrac"), by = c("sample_1", "sample_2")) |>
    dplyr::left_join(get_lower_tri_df(geo, "geographic_km"), by = c("sample_1", "sample_2")) |>
    dplyr::left_join(get_lower_tri_df(temp, "delta_temp"), by = c("sample_1", "sample_2")) |>
    dplyr::left_join(get_lower_tri_df(sym_conn, "ocean_connectivity"), by = c("sample_1", "sample_2")) |>
    dplyr::left_join(get_lower_tri_df(ocean_dist, "ocean_distance"), by = c("sample_1", "sample_2")) |>
    dplyr::filter(!is.na(ocean_connectivity))

  stats <- list(
    mantel_bray_ocean = vegan::mantel(as.dist(bray), as.dist(ocean_dist), permutations = 9999),
    mantel_unifrac_ocean = vegan::mantel(as.dist(unifrac), as.dist(ocean_dist), permutations = 9999),
    partial_bray_ocean_geo = vegan::mantel.partial(as.dist(bray), as.dist(ocean_dist), as.dist(geo), permutations = 9999),
    partial_unifrac_ocean_geo = vegan::mantel.partial(as.dist(unifrac), as.dist(ocean_dist), as.dist(geo), permutations = 9999),
    partial_bray_ocean_temp = vegan::mantel.partial(as.dist(bray), as.dist(ocean_dist), as.dist(temp), permutations = 9999),
    partial_unifrac_ocean_temp = vegan::mantel.partial(as.dist(unifrac), as.dist(ocean_dist), as.dist(temp), permutations = 9999)
  )

  plot <- ggplot(pair_df, aes(ocean_connectivity, 1 - bray)) +
    geom_point(size = 2, alpha = 0.75) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    labs(x = "Symmetric oceanographic connectivity", y = "Bray-Curtis similarity") +
    theme_bw()

  list(pairwise = pair_df, stats = stats, plot = plot, connectivity = sym_conn)
}
