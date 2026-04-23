source("/home/userbio/Sea4Blue/scripts/utils/helpers.R")

run_community_composition <- function(cfg = sea4blue_config(), top_n = 8, rank = "Phylum") {
  load_packages(c("phyloseq", "dplyr", "tidyr", "ggplot2", "forcats", "scales", "tibble"))
  ps <- readRDS(cfg$paths$phyloseq_all)
  meta <- read.csv(cfg$paths$metadata, check.names = FALSE)
  keep <- meta$sample.id[grepl("^Lapo[0-9]+$", meta$sample.id)]
  ps <- phyloseq::subset_samples(ps, sample_names(ps) %in% keep)
  meta <- meta[match(sample_names(ps), meta$sample.id), ]
  otu <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) otu <- t(otu)
  rel <- sweep(otu, 2, colSums(otu), "/")
  tax <- as.data.frame(phyloseq::tax_table(ps), stringsAsFactors = FALSE)
  rank_vec <- tax[[rank]]
  rank_vec[is.na(rank_vec) | rank_vec == ""] <- "Unclassified"
  agg <- rowsum(rel, rank_vec)
  avg <- sort(rowMeans(agg), decreasing = TRUE)
  keep_taxa <- names(avg)[seq_len(min(top_n, length(avg)))]
  other <- setdiff(rownames(agg), keep_taxa)
  if (length(other)) {
    agg2 <- rbind(agg[keep_taxa, , drop = FALSE], Others = colSums(agg[other, , drop = FALSE]))
  } else {
    agg2 <- agg[keep_taxa, , drop = FALSE]
  }
  long <- as.data.frame(agg2) |>
    tibble::rownames_to_column("Taxon") |>
    tidyr::pivot_longer(-Taxon, names_to = "sample.id", values_to = "RelativeAbundance") |>
    dplyr::left_join(meta, by = "sample.id") |>
    dplyr::mutate(
      Taxon = factor(Taxon, levels = c(keep_taxa, if ("Others" %in% Taxon) "Others")),
      sample.id = factor(sample.id, levels = meta$sample.id),
      Zona = factor(Zona, levels = c("ANW", "ANC", "ANE"))
    )
  p <- ggplot(long, aes(sample.id, RelativeAbundance, fill = Taxon)) +
    geom_col(width = 0.94) +
    facet_grid(~ Zona, scales = "free_x", space = "free_x") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Stations ordered along the transect", y = paste(rank, "relative abundance")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), panel.spacing.x = unit(0.4, "lines"))
  utils::write.csv(long, file.path(cfg$tables_dir, "community_composition_phylum.csv"), row.names = FALSE)
  ggplot2::ggsave(file.path(cfg$figures_dir, "community_composition_phylum.pdf"), p, width = 11.69, height = 8.27, units = "in", device = grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(cfg$project_root, "paper", "draft", "figures", "community_composition_phylum.png"), p, width = 11.69, height = 8.27, units = "in", dpi = 300, bg = "white")
  invisible(list(plot = p, data = long, top = avg))
}

if (sys.nframe() == 0) {
  x <- run_community_composition()
  print(round(head(x$top, 10), 4))
}
