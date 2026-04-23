source(file.path('/home/userbio/Sea4Blue', 'scripts', '00_config.R'))

sea4blue_zone_colors <- function() {
  c(ANW = '#4280fc', ANC = '#ffb452', ANE = '#f7170a')
}

module_sort_key <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  out[is.na(out)] <- Inf
  out
}

load_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop('Missing packages: ', paste(missing, collapse = ', '))
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

load_phyloseq_atlantic <- function(cfg = sea4blue_config()) {
  load_packages(c('phyloseq', 'dplyr', 'microbiome', 'metagenomeSeq', 'tibble', 'readr'))
  ps <- readRDS(cfg$paths$phyloseq_rds)
  env_data <- read.csv(cfg$paths$env_csv, stringsAsFactors = FALSE)

  otu_mat <- as(phyloseq::otu_table(ps), 'matrix')
  if (!phyloseq::taxa_are_rows(ps)) {
    otu_mat <- t(otu_mat)
  }
  phyloseq::otu_table(ps) <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)

  phyloseq_obj <- ps
  tax_df <- as.data.frame(phyloseq::tax_table(phyloseq_obj), stringsAsFactors = FALSE)
  keep_taxa <-
    (is.na(tax_df$Family) | tax_df$Family != 'Mitochondria') &
    (is.na(tax_df$Class) | tax_df$Class != 'Chloroplast') &
    (is.na(tax_df$Order) | tax_df$Order != 'Chloroplast')
  phyloseq_obj <- phyloseq::prune_taxa(keep_taxa, phyloseq_obj)

  meta_df <- data.frame(phyloseq::sample_data(phyloseq_obj), stringsAsFactors = FALSE)
  keep_samples <- !(meta_df$Nominativo.campione. %in% c('Positive', 'Negative'))
  phyloseq_obj <- phyloseq::prune_samples(keep_samples, phyloseq_obj)
  phyloseq::sample_names(phyloseq_obj) <- phyloseq::sample_data(phyloseq_obj)$Nominativo.campione.

  merged_metadata <- merge(
    microbiome::meta(phyloseq::sample_data(phyloseq_obj)),
    env_data,
    by = 'Nominativo.campione.',
    all.x = TRUE,
    sort = FALSE
  )
  rownames(merged_metadata) <- merged_metadata$Nominativo.campione.
  merged_metadata <- merged_metadata[phyloseq::sample_names(phyloseq_obj), , drop = FALSE]
  phyloseq::sample_data(phyloseq_obj) <- merged_metadata

  doubleton <- phyloseq::genefilter_sample(
    phyloseq_obj,
    phyloseq::filterfun_sample(function(x) x > 1),
    A = 1
  )
  doubleton <- phyloseq::prune_taxa(doubleton, phyloseq_obj)

  metaseq <- phyloseq::phyloseq_to_metagenomeSeq(doubleton)
  p <- metagenomeSeq::cumNormStat(metaseq)
  data_cumnorm <- metagenomeSeq::cumNorm(metaseq, p = p)
  data_css <- metagenomeSeq::MRcounts(data_cumnorm, norm = TRUE, log = TRUE)

  phyloseq_obj_css <- phyloseq_obj
  phyloseq::otu_table(phyloseq_obj_css) <- phyloseq::otu_table(data_css, taxa_are_rows = TRUE)
  phyloseq::sample_names(phyloseq_obj_css) <- phyloseq::sample_data(phyloseq_obj_css)$Nominativo.campione.

  physeq_normalized <- phyloseq_obj_css
  phyloseq::otu_table(physeq_normalized) <- phyloseq::otu_table(
    as.matrix(round(phyloseq::otu_table(phyloseq_obj_css))),
    taxa_are_rows = TRUE
  )

  list(
    phyloseq_obj = phyloseq_obj,
    phyloseq_obj_css = phyloseq_obj_css,
    physeq_normalized = physeq_normalized,
    sample_data = data.frame(phyloseq::sample_data(phyloseq_obj_css))
  )
}

classify_mechanism <- function(bNTI, RC_BC) {
  dplyr::case_when(
    bNTI > 2 ~ 'Heterogeneous Selection',
    bNTI < -2 ~ 'Homogeneous Selection',
    RC_BC > 0.95 ~ 'Dispersal Limitation',
    RC_BC < -0.95 ~ 'Homogenising Dispersal',
    TRUE ~ 'Drift'
  )
}
