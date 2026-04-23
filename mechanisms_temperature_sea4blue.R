get_mechanism_prop_with_thetao <- function(datalist, bNTI, RC_BC) {
  
  bNTI_mod <- bNTI %>%
    magrittr::set_colnames(rownames(.)) %>%
    as.matrix()
  bNTI_mod[lower.tri(bNTI_mod)] <- NA
  bNTI_mod[diag(bNTI_mod)] <- NA
  
  RC_BC_mod <- RC_BC %>%
    magrittr::set_rownames(rownames(bNTI_mod)) %>%
    magrittr::set_colnames(colnames(bNTI_mod)) %>%
    as.matrix()
  RC_BC_mod[lower.tri(RC_BC_mod)] <- NA
  RC_BC_mod[diag(RC_BC_mod)] <- NA
  
  # dettaglio edge-level
  detailled <- bNTI_mod %>%
    reshape2::melt() %>%
    with(., cbind(., reshape2::melt(RC_BC_mod))) %>%
    magrittr::set_colnames(c("From_Sample", "To_Sample", "bNTI", "1", "2", "RC_BC")) %>%
    dplyr::select(From_Sample, To_Sample, bNTI, RC_BC) %>%
    filter(!is.na(bNTI)) %>%
    mutate(Mechanism = ifelse(bNTI > 2, "Heterogeneous Selection", 
                              ifelse(bNTI < -2, "Homogeneous Selection",
                                     ifelse(RC_BC < 0.95 & RC_BC > -0.95, "Drift",
                                            ifelse(RC_BC > 0.95, "Dispersal Limitation", "Homogenising Dispersal"))))) 
  datalist$Meta_Data %>% as.data.frame() %>% colnames()
  # aggiungo ΔT usando thetao
  Meta_Data <- datalist$Meta_Data %>% as.data.frame()
  detailled <- detailled %>%
    left_join(Meta_Data %>% dplyr::select(Sample_ID, thetao) %>% dplyr::rename(From_Temp = thetao),
              by = c("From_Sample" = "Sample_ID")) %>%
    left_join(Meta_Data %>% dplyr::select(Sample_ID, thetao) %>% dplyr::rename(To_Temp = thetao),
              by = c("To_Sample" = "Sample_ID")) %>%
    mutate(
      Temp_Diff = abs(as.numeric(From_Temp) - as.numeric(To_Temp)),
      Temp_Diff_Grp = seq(0, 28, 1)[findInterval(Temp_Diff, seq(0, 28, 1))]
    )
  
  return(detailled)
}



detailled_Atlantic <- get_mechanism_prop_with_thetao(
  datalist_Atlantic,
  bNTI = as.matrix(read.csv("output/Community_Mechanisms/Prokaryotes_Atlantic_weighted_bNTI.csv", row.names = 1)),
  RC_BC = as.matrix(read.csv("output/Community_Mechanisms/Raup_Crick_Prok.csv", row.names = 1))
)
mech_temp_Atlantic <- detailled_Atlantic %>%
  dplyr::group_by(Temp_Diff_Grp, Mechanism) %>%
  dplyr::summarise(N = n(), .groups = "drop") %>%
  dplyr::group_by(Temp_Diff_Grp) %>%
  dplyr::mutate(Prop = N / sum(N)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Mechanism = factor(Mechanism,
                            levels = c("Homogeneous Selection", "Heterogeneous Selection",
                                       "Homogenising Dispersal", "Dispersal Limitation", "Drift"),
                            labels = c("Homogeneous selection", "Heterogeneous selection",
                                       "Homogenising dispersal", "Dispersal limitation", "Drift")))

ggplot(mech_temp_Atlantic,
       aes(x = Temp_Diff_Grp, y = Prop*100, fill = Mechanism)) +
  geom_bar(stat = "identity", colour = "black", size = 0.2) +
  scale_fill_manual(values = cbbPalette[c(7,5,3,6,4)]) +
  labs(x = "Temperature difference (°C)", y = "Proportion of mechanisms (%)") +
  theme_bw()

