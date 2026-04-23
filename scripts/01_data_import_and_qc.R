source("/home/userbio/Sea4Blue/scripts/utils/helpers.R")

prepare_sea4blue_objects <- function(cfg = sea4blue_config()) {
  ensure_result_dirs(cfg)
  objs <- load_phyloseq_atlantic(cfg)
  objs
}

export_modular_seascape_bundle <- function(objs, out_dir = file.path("/home/userbio/Sea4Blue", "results", "derived_data", "Data")) {
  load_packages(c("readr", "tibble", "Biostrings", "ape", "dplyr"))
  dir.create(file.path(out_dir, "Meta_Data", "Prok"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "Count_Data", "Processed", "Prok"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "Count_Data", "Fasta", "Prok"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "Trees", "Prok"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "Taxonomy"), recursive = TRUE, showWarnings = FALSE)

  sample_data(objs$phyloseq_obj) %>%
    as_tibble() %>%
    dplyr::rename(Sample_ID = Nominativo.campione.) %>%
    dplyr::mutate(Sample_ID = as.character(Sample_ID)) %>%
    dplyr::relocate(Sample_ID, .before = 1) %>%
    readr::write_tsv(file.path(out_dir, "Meta_Data", "Prok", "Meta_Data.tsv"))

  otu <- otu_table(objs$phyloseq_obj) %>% as.data.frame() %>% tibble::rownames_to_column("OTU_ID")
  readr::write_tsv(otu, file.path(out_dir, "Count_Data", "Processed", "Prok", "Full_Prok_Count.tsv"))
  Biostrings::writeXStringSet(refseq(objs$phyloseq_obj), file.path(out_dir, "Count_Data", "Fasta", "Prok", "Full_Prok_Sequences.fasta"))
  ape::write.tree(phy_tree(objs$phyloseq_obj), file.path(out_dir, "Trees", "Prok", "Prok_Combined.tree"))
  tax_table(objs$phyloseq_obj) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("OTU_ID") %>%
    readr::write_tsv(file.path(out_dir, "Taxonomy", "Prok_Taxonomy.tsv"))

  invisible(out_dir)
}
