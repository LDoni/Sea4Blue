source('/home/userbio/Sea4Blue/scripts/utils/helpers.R')

compute_phylo_similarity <- function(phyloseq_obj, use.cores = 2, cor.use = 'na.or.complete') {
  ot <- phyloseq::otu_table(phyloseq_obj)
  if (!phyloseq::taxa_are_rows(phyloseq_obj)) ot <- t(ot)
  tree <- phyloseq::phy_tree(phyloseq_obj)
  cophen <- ape::cophenetic.phylo(tree)
  cophen <- cophen[rownames(ot), rownames(ot), drop = FALSE]
  tmp_s <- cor.par(cophen, method = 'pearson', use = cor.use, use.cores = use.cores)
  0.5 * (tmp_s + 1)
}

compute_tina_pina <- function(ps, method = c('TINA_weighted', 'PINA_weighted'), size.thresh = 1, pseudocount = 1e-6, nblocks = 25, use.cores = 2, cor.use = 'na.or.complete') {
  method <- match.arg(method)
  source('/media/shared1/atlantic_sea4blue/16S_Illumina/Lapo_18May/Pacific_Bacterioplankton/R/Functions_Similarity_Indices.R')
  if (method == 'TINA_weighted') {
    distance_wrapper2(ps, method = method, size.thresh = size.thresh, pseudocount = pseudocount, nblocks = nblocks, use.cores = use.cores, cor.use = cor.use)
  } else {
    s_phylo <- compute_phylo_similarity(ps, use.cores = use.cores, cor.use = cor.use)
    community.similarity.corr.par(phyloseq::otu_table(ps), S = s_phylo, distance = method, blocksize = nblocks, use.cores = use.cores)
  }
}

run_tina_pina <- function(cfg = sea4blue_config(), use.cores = 2) {
  load_packages(c('dplyr', 'ggplot2', 'phyloseq', 'vegan', 'geosphere', 'ape', 'Matrix', 'foreach', 'doMC', 'bigmemory', 'magrittr', 'plyr'))
  zone_colors <- sea4blue_zone_colors()
  objs <- load_phyloseq_atlantic(cfg)
  ps <- objs$physeq_normalized
  ps <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > 0, ps)
  meta <- objs$sample_data
  sample_ids <- rownames(meta)
  meta <- meta |>
    dplyr::mutate(
      Zona = factor(as.character(Zona), levels = names(zone_colors)),
      Longitude.x = as.numeric(Longitude.x),
      Latitude.x = as.numeric(Latitude.x)
    )

  tina <- compute_tina_pina(ps, method = 'TINA_weighted', use.cores = use.cores)
  pina <- compute_tina_pina(ps, method = 'PINA_weighted', use.cores = use.cores)
  rownames(tina) <- colnames(tina) <- sample_ids
  rownames(pina) <- colnames(pina) <- sample_ids

  tina_dist <- as.dist(tina)
  pina_dist <- as.dist(pina)
  bray <- as.matrix(phyloseq::distance(ps, method = 'bray'))
  coords <- as.matrix(meta[, c('Longitude.y', 'Latitude.y')])
  rownames(coords) <- sample_ids
  geo <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
  rownames(geo) <- colnames(geo) <- sample_ids

  ord_tina <- ape::pcoa(tina_dist)
  ord_pina <- ape::pcoa(pina_dist)
  tina_df <- data.frame(ord_tina$vectors[, 1:2], Zona = meta$Zona, Sample_ID = sample_ids)
  pina_df <- data.frame(ord_pina$vectors[, 1:2], Zona = meta$Zona, Sample_ID = sample_ids)
  colnames(tina_df)[1:2] <- c('Axis1', 'Axis2')
  colnames(pina_df)[1:2] <- c('Axis1', 'Axis2')

  tina_plot <- ggplot(tina_df, aes(Axis1, Axis2, color = Zona)) +
    geom_point(size = 3) +
    stat_ellipse(aes(group = Zona), linewidth = 0.6, show.legend = FALSE) +
    scale_color_manual(values = zone_colors, drop = FALSE) +
    labs(x = 'PCoA1', y = 'PCoA2', color = 'Zone') +
    theme_bw()

  pina_plot <- ggplot(pina_df, aes(Axis1, Axis2, color = Zona)) +
    geom_point(size = 3) +
    stat_ellipse(aes(group = Zona), linewidth = 0.6, show.legend = FALSE) +
    scale_color_manual(values = zone_colors, drop = FALSE) +
    labs(x = 'PCoA1', y = 'PCoA2', color = 'Zone') +
    theme_bw()

  permanova_tina <- vegan::adonis2(as.matrix(tina_dist) ~ Zona + Longitude.x + thetao, data = meta, permutations = 999)
  permanova_pina <- vegan::adonis2(as.matrix(pina_dist) ~ Zona + Longitude.x + thetao, data = meta, permutations = 999)

  lower_df <- function(mat, value_name) {
    idx <- which(lower.tri(mat), arr.ind = TRUE)
    out <- tibble::tibble(sample_1 = rownames(mat)[idx[,1]], sample_2 = colnames(mat)[idx[,2]], value = mat[idx])
    names(out)[3] <- value_name
    out
  }

  tina_pair <- lower_df(as.matrix(tina_dist), 'tina_dist') |>
    dplyr::left_join(lower_df(as.matrix(pina_dist), 'pina_dist'), by = c('sample_1','sample_2')) |>
    dplyr::left_join(lower_df(geo, 'geo_km'), by = c('sample_1','sample_2')) |>
    dplyr::left_join(lower_df(bray, 'bray'), by = c('sample_1','sample_2'))

  tina_mantel <- vegan::mantel(tina_dist, as.dist(geo), permutations = 9999)
  pina_mantel <- vegan::mantel(pina_dist, as.dist(geo), permutations = 9999)

  tina_decay_plot <- ggplot(tina_pair, aes(geo_km, tina_dist)) +
    geom_point(alpha = 0.75, size = 2) +
    geom_smooth(method = 'lm', se = TRUE, color = 'black') +
    labs(x = 'Geographic distance (km)', y = 'TINA weighted distance') +
    theme_bw()

  pina_decay_plot <- ggplot(tina_pair, aes(geo_km, pina_dist)) +
    geom_point(alpha = 0.75, size = 2) +
    geom_smooth(method = 'lm', se = TRUE, color = 'black') +
    labs(x = 'Geographic distance (km)', y = 'PINA weighted distance') +
    theme_bw()

  compare_plot <- ggplot(tina_pair, aes(bray, tina_dist)) +
    geom_point(alpha = 0.75, size = 2, color = '#2b8cbe') +
    geom_smooth(method = 'lm', se = TRUE, color = 'black') +
    labs(x = 'Bray-Curtis distance', y = 'TINA weighted distance') +
    theme_bw()

  list(
    tina = tina,
    pina = pina,
    tina_plot = tina_plot,
    pina_plot = pina_plot,
    tina_decay_plot = tina_decay_plot,
    pina_decay_plot = pina_decay_plot,
    compare_plot = compare_plot,
    pairwise = tina_pair,
    tina_permanova = permanova_tina,
    pina_permanova = permanova_pina,
    tina_mantel = tina_mantel,
    pina_mantel = pina_mantel
  )
}
