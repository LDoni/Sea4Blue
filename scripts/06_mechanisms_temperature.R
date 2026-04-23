source("/home/userbio/Sea4Blue/scripts/utils/helpers.R")

run_mechanisms_temperature <- function(cfg = sea4blue_config()) {
  load_packages(c("dplyr", "ggplot2", "reshape2", "magrittr", "readr"))
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

  detailled <- reshape2::melt(bNTI_mod) %>%
    cbind(reshape2::melt(RC_BC_mod)) %>%
    magrittr::set_colnames(c("From_Sample", "To_Sample", "bNTI", "x1", "x2", "RC_BC")) %>%
    dplyr::select(From_Sample, To_Sample, bNTI, RC_BC) %>%
    dplyr::filter(!is.na(bNTI)) %>%
    dplyr::mutate(Mechanism = classify_mechanism(bNTI, RC_BC)) %>%
    dplyr::left_join(meta %>% dplyr::select(Sample_ID, thetao) %>% dplyr::rename(From_Temp = thetao), by = c("From_Sample" = "Sample_ID")) %>%
    dplyr::left_join(meta %>% dplyr::select(Sample_ID, thetao) %>% dplyr::rename(To_Temp = thetao), by = c("To_Sample" = "Sample_ID")) %>%
    dplyr::mutate(Temp_Diff = abs(as.numeric(From_Temp) - as.numeric(To_Temp)), Temp_Diff_Grp = seq(0, 28, 1)[findInterval(Temp_Diff, seq(0, 28, 1))])

  mech_temp <- detailled %>%
    dplyr::group_by(Temp_Diff_Grp, Mechanism) %>%
    dplyr::summarise(N = n(), .groups = "drop") %>%
    dplyr::group_by(Temp_Diff_Grp) %>%
    dplyr::mutate(Prop = N / sum(N)) %>%
    dplyr::ungroup()

  plot <- ggplot(mech_temp, aes(x = Temp_Diff_Grp, y = Prop * 100, fill = Mechanism)) +
    geom_bar(stat = "identity", colour = "black", linewidth = 0.2) +
    labs(x = "Temperature difference (°C)", y = "Proportion of mechanisms (%)") +
    theme_bw()

  list(detail = detailled, summary = mech_temp, plot = plot)
}
