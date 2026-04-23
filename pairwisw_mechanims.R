get_mechanism_pairwise <- function(bNTI, RC_BC) {
  
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
  
  pairwise_df <- bNTI_mod %>%
    reshape2::melt() %>%
    dplyr::rename(Sample_1 = Var1, Sample_2 = Var2, bNTI = value) %>%
    dplyr::mutate(RC_BC = reshape2::melt(RC_BC_mod)$value) %>%
    dplyr::filter(!is.na(bNTI)) %>%
    dplyr::mutate(
      Mechanism = dplyr::case_when(
        bNTI >  2  ~ "Heterogeneous Selection",
        bNTI < -2  ~ "Homogeneous Selection",
        RC_BC >  0.95  ~ "Dispersal Limitation",
        RC_BC < -0.95  ~ "Homogenising Dispersal",
        TRUE ~ "Drift"
      ),
      Mechanism_Group = factor(
        dplyr::case_when(
          Mechanism %in% c("Homogeneous Selection", "Heterogeneous Selection") ~ "Selection",
          Mechanism %in% c("Dispersal Limitation", "Homogenising Dispersal") ~ "Dispersal",
          TRUE ~ "Drift"
        ),
        levels = c("Drift", "Dispersal", "Selection")
      )
    )
  
  return(pairwise_df)
}


get_mechanism_prop <- function(pairwise_df) {
  
  pairwise_df %>%
    dplyr::count(Mechanism, Mechanism_Group) %>%
    dplyr::mutate(Num = n / sum(n))
}



pairwise_mech_Atlanti <- get_mechanism_pairwise(
  bNTI = read.csv("output/Community_Mechanisms/Prokaryotes_Atlantic_weighted_bNTI.csv",
                  row.names = 1) %>% as.matrix(),
  RC_BC = read.csv("output/Community_Mechanisms/Raup_Crick_Prok.csv",
                   row.names = 1) %>% as.matrix()
)

mechanism_prop_Atlantic <- get_mechanism_prop(pairwise_mech_Atlanti)



pairwise_mech_Atlantic_clean <- pairwise_mech_Atlanti %>%
  dplyr::filter(Sample_1 != Sample_2) %>%
  dplyr::mutate(
    Sample_1 = as.character(Sample_1),
    Sample_2 = as.character(Sample_2),
    site_min = pmin(Sample_1, Sample_2),
    site_max = pmax(Sample_1, Sample_2)
  )




