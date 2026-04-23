source("/home/userbio/Sea4Blue/scripts/utils/helpers.R")

get_mechanism_proportions <- function(cfg = sea4blue_config()) {
  load_packages(c("dplyr", "ggplot2", "reshape2", "readr", "magrittr"))
  objs <- load_phyloseq_atlantic(cfg)
  meta <- objs$sample_data
  meta$Sample_ID <- rownames(meta)

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

  pairwise <- reshape2::melt(bNTI_mod) %>%
    cbind(reshape2::melt(RC_BC_mod)) %>%
    magrittr::set_colnames(c("Sample_ID", "To_Sample", "bNTI", "x1", "x2", "RC_BC")) %>%
    dplyr::select(Sample_ID, To_Sample, bNTI, RC_BC) %>%
    dplyr::filter(!is.na(bNTI)) %>%
    dplyr::mutate(Mechanism = classify_mechanism(bNTI, RC_BC))

  summary_tbl <- pairwise %>%
    dplyr::count(Mechanism, name = "Num") %>%
    dplyr::mutate(Num = Num / sum(Num)) %>%
    dplyr::mutate(Mechanism_Group = ordered(ifelse(Mechanism %in% c("Homogeneous Selection", "Heterogeneous Selection"), "Selection", ifelse(Mechanism %in% c("Homogenising Dispersal", "Dispersal Limitation"), "Dispersal", "Drift")), levels = c("Drift", "Dispersal", "Selection")))

  plot <- ggplot(summary_tbl, aes(x = Mechanism_Group, y = Num * 100, fill = Mechanism)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(y = "Proportion of mechanisms (%)", x = "") +
    theme_bw()

  list(pairwise = pairwise, summary = summary_tbl, plot = plot)
}
