# Get relative importance of ecological mechanisms from Stegen et al. framework
# Requires bNTI table and RC_BC matrix
# -> See R-scripts folder

library(tidyverse)

source("Modular_Seascape/R/Datalist_Wrangling_Functions.R")
source("Modular_Seascape/R/Import_Data.R")


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

get_mechanism_prop <- function(datalist, bNTI, RC_BC) {
  
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
  
  merged <- bNTI_mod %>%
    reshape2::melt() %>%
    with(., cbind(., reshape2::melt(RC_BC_mod))) %>%
    magrittr::set_colnames(c("Sample_ID", "To_Sample", "bNTI", "1", "2", "RC_BC")) %>%
    dplyr::select(Sample_ID, To_Sample, bNTI, RC_BC) %>%
    filter(!is.na(bNTI)) %>%
    mutate(Mechanism = ifelse(bNTI > 2, "Heterogeneous Selection", 
                              ifelse(bNTI < -2, "Homogeneous Selection",
                                     ifelse(RC_BC < 0.95 & RC_BC > -0.95, "Drift",
                                            ifelse(RC_BC > 0.95, "Dispersal Limitation", "Homogenising Dispersal"))))) %>%
    left_join(., select(datalist$Meta_Data, Sample_ID, Depth_Grp), by = "Sample_ID") %>%
    dplyr::rename("From_Depth_Grp" = "Depth_Grp", "From_Sample" = "Sample_ID", "Sample_ID" = "To_Sample") %>%
    left_join(., select(datalist$Meta_Data, Sample_ID, Depth_Grp), by = "Sample_ID") %>%
    dplyr::rename("To_Depth_Grp" = "Depth_Grp", "To_Sample" = "Sample_ID") %>%
    filter(From_Depth_Grp == To_Depth_Grp) %>%
    group_by(Mechanism, From_Depth_Grp) %>%
    summarize(Num = n()) %>%
    group_by(From_Depth_Grp) %>%
    mutate(Num = Num/sum(Num)) %>%
    ungroup() %>%
    mutate(Mechanism_Group = ordered(ifelse(Mechanism == "Homogeneous Selection" | Mechanism == "Heterogeneous Selection", "Selection",
                                            ifelse(Mechanism == "Homogenising Dispersal" | Mechanism == "Dispersal Limitation", "Dispersal",
                                                   "Drift")),
                                     levels = c("Drift", "Dispersal", "Selection"))) %>%
    mutate(Mechanism = ordered(Mechanism, levels = c("Homogeneous Selection", "Heterogeneous Selection", 
                                                     "Homogenising Dispersal", "Dispersal Limitation", "Drift")))
  
  return(merged)
  
}


### quella che ho usato !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!_--->>>>>
get_mechanism_prop2 <- function(datalist, bNTI, RC_BC) {
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
  
  merged <- bNTI_mod %>%
    reshape2::melt() %>%
    with(., cbind(., reshape2::melt(RC_BC_mod))) %>%
    magrittr::set_colnames(c("Sample_ID", "To_Sample", "bNTI", "1", "2", "RC_BC")) %>%
    dplyr::select(Sample_ID, To_Sample, bNTI, RC_BC) %>%
    filter(!is.na(bNTI)) %>%
    dplyr::mutate(Mechanism = ifelse(bNTI > 2, "Heterogeneous Selection", 
                              ifelse(bNTI < -2, "Homogeneous Selection",
                                     ifelse(RC_BC < 0.95 & RC_BC > -0.95, "Drift",
                                            ifelse(RC_BC > 0.95, "Dispersal Limitation", "Homogenising Dispersal"))))) %>%
    group_by(Mechanism) %>%
    dplyr::summarize(Num = n()) %>%
    ungroup() %>%
    dplyr::mutate(Num = Num/sum(Num)) %>%
    dplyr::mutate(Mechanism_Group = ordered(ifelse(Mechanism == "Homogeneous Selection" | Mechanism == "Heterogeneous Selection", "Selection",
                                            ifelse(Mechanism == "Homogenising Dispersal" | Mechanism == "Dispersal Limitation", "Dispersal",
                                                   "Drift")),
                                     levels = c("Drift", "Dispersal", "Selection"))) %>%
    dplyr::mutate(Mechanism = ordered(Mechanism, levels = c("Homogeneous Selection", "Heterogeneous Selection", 
                                                     "Homogenising Dispersal", "Dispersal Limitation", "Drift")))
  
  return(merged)
}









datalist_Atlantic <- import_data2("Data", kingdom = "Prok", abundance_filter = T, min_counts = 2000) 
 
Meta_Data <- datalist_Atlantic$Meta_Data %>% as.data.frame()

class(Meta_Data)



merged_Atlantic <- bind_rows(
          get_mechanism_prop2(datalist_Atlantic, 
                             read.csv("output/Community_Mechanisms/Prokaryotes_Atlantic_weighted_bNTI.csv", 
                                      row.names = 1)%>% as.matrix(), 
                             read.csv("output/Community_Mechanisms/Raup_Crick_Prok.csv", 
                                      row.names = 1)%>% as.matrix())) 
          
 
# A tibble: 3 × 3
Mechanism               Num Mechanism_Group
<ord>                 <dbl> <ord>          
  1 Dispersal Limitation  0.125 Dispersal      
2 Drift                 0.158 Drift          
3 Homogeneous Selection 0.717 Selection


  ggplot(merged_Atlantic, aes(x = Mechanism_Group, y = Num*100, fill = Mechanism)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(y = "Proportion of Mechanisms (%)", x = "") +
    scale_fill_manual(values = cbbPalette[c(7,5, 3, 6, 4)]) +
    theme_bw() +
    theme(legend.position="right") +
    guides(fill = guide_legend(nrow = 5, title.position = "top")) +
    ylim(c(0, 75))

  
  
  
  library(reshape2)
  library(dplyr)
  library(ggplot2)
  
  # Carica bNTI e RC_BC (li hai già)
  bNTI <- read.csv("output/Community_Mechanisms/Prokaryotes_Atlantic_weighted_bNTI.csv",
                   row.names = 1) %>% as.matrix()
  
  RC_BC <- read.csv("output/Community_Mechanisms/Raup_Crick_Prok.csv",
                    row.names = 1) %>% as.matrix()
  
  # Converti in formato long
  bNTI_df <- melt(bNTI, varnames = c("Sample1","Sample2"), value.name="bNTI")
  RC_df   <- melt(RC_BC, varnames = c("Sample1","Sample2"), value.name="RC_BC")
  
  # Merge
  plot_df <- bNTI_df %>%
    left_join(RC_df, by = c("Sample1","Sample2")) %>%
    filter(Sample1 != Sample2)  # rimuovi la diagonale
  
  ggplot(plot_df, aes(x = bNTI, y = RC_BC)) +
    geom_point(size = 2, alpha = 0.65) +
    geom_vline(xintercept = c(-2, 2), 
               linetype = "dashed", color = "red", size = 1) +
    geom_hline(yintercept = c(-0.95, 0.95), 
               linetype = "dashed", color = "blue", size = 1) +
    theme_bw() +
    labs(
      x = "βNTI",
      y = "RC_bray",
      title = "Stegen Framework: βNTI vs RC_bray",
      subtitle = "Boundary lines: βNTI = ±2, RC_bray = ±0.95"
    )
  
  
  
  ## bootstrap per confidenza
  bNTI <- read.csv(
    "output/Community_Mechanisms/Prokaryotes_Atlantic_weighted_bNTI.csv",
    row.names = 1
  ) |> as.matrix()
  
  RC_BC <- read.csv(
    "output/Community_Mechanisms/Raup_Crick_Prok.csv",
    row.names = 1
  ) |> as.matrix()
  
  
  pairwise_mechanisms <- function(bNTI, RC_BC){
    
    # controlli base
    stopifnot(all(dim(bNTI) == dim(RC_BC)))
    stopifnot(all(rownames(bNTI) == rownames(RC_BC)))
    
    # usa solo triangolo superiore
    idx <- which(upper.tri(bNTI), arr.ind = TRUE)
    
    df <- data.frame(
      Sample1 = rownames(bNTI)[idx[,1]],
      Sample2 = colnames(bNTI)[idx[,2]],
      bNTI    = bNTI[idx],
      RC_BC   = RC_BC[idx],
      stringsAsFactors = FALSE
    )
    
    # classificazione meccanismi
    df$Mechanism <- dplyr::case_when(
      df$bNTI > 2 ~ "Heterogeneous Selection",
      df$bNTI < -2 ~ "Homogeneous Selection",
      df$RC_BC > 0.95 ~ "Dispersal Limitation",
      df$RC_BC < -0.95 ~ "Homogenising Dispersal",
      TRUE ~ "Drift"
    )
    
    return(df)
  }
  pairs <- pairwise_mechanisms(bNTI, RC_BC)
  
  head(pairs)  
  dim(pairs)
  
  
  
  bootstrap_mechanisms <- function(pairs, n_boot = 3000){
    
    mechanisms <- unique(pairs$Mechanism)
    
    boot_mat <- replicate(n_boot, {
      
      samp <- pairs[sample(nrow(pairs), replace = TRUE), ]
      
      tab <- table(factor(samp$Mechanism, levels = mechanisms))
      
      prop <- tab / sum(tab)
      
      return(prop)
      
    })
    
    boot_mat <- t(boot_mat)
    
    data.frame(
      Mechanism = mechanisms,
      Mean = colMeans(boot_mat),
      Lower95 = apply(boot_mat, 2, quantile, 0.025),
      Upper95 = apply(boot_mat, 2, quantile, 0.975)
    )
  }
  
  prop.table(table(pairs$Mechanism))
  table(pairs$Mechanism)
  
  set.seed(123)
  
  boot_results <- bootstrap_mechanisms(pairs, n_boot = 3000)
  
  boot_results
  
 # boot_results
#  Mechanism       Mean    Lower95    Upper95
#  Dispersal Limitation   Dispersal Limitation 0.14353333 0.07619048 0.21904762
#  Drift                                 Drift 0.05009841 0.00952381 0.08571429
#  Homogeneous Selection Homogeneous Selection 0.80847937 0.73333333 0.88571429
  
  
  
  ggplot(boot_results,
         aes(x = Mechanism, y = Mean)) +
    
    geom_col(fill = "steelblue") +
    
    geom_errorbar(
      aes(ymin = Lower95,
          ymax = Upper95),
      width = 0.2
    ) +
    
    theme_bw() +
    ylab("Proportion")
    
  
  
  
  
  ###meccanismi in base a distanza
  
  
  head(coord.df)
  coord.df$id <- tolower(coord.df$id)
  library(geosphere)
  
  coords <- coord.df[, c("lon", "lat")]
  rownames(coords) <- coord.df$id
  
  geo_dist <- distm(coords, fun = distHaversine) / 1000  # km
  rownames(geo_dist) <- tolower(coord.df$id)
  colnames(geo_dist) <- tolower(coord.df$id)
  
  # Triangolo superiore
  idx <- which(upper.tri(geo_dist), arr.ind = TRUE)
  geo_df <- data.frame(
    Sample1 = rownames(geo_dist)[idx[,1]],
    Sample2 = colnames(geo_dist)[idx[,2]],
    geo_dist = geo_dist[idx]
  )
  
  # Assicurati che i nomi in 'pairs' siano minuscoli
  pairs$Sample1 <- tolower(pairs$Sample1)
  pairs$Sample2 <- tolower(pairs$Sample2)
  
  # Unisci distanze geografiche
  pairs_geo <- left_join(pairs, geo_df, by = c("Sample1","Sample2"))
  

  breaks <- seq(0, 6500, length.out = 6)  # 5 bin di ~1300 km ciascuno
  pairs_geo$distance_bin <- cut(pairs_geo$geo_dist,
                                breaks = breaks,
                                include.lowest = TRUE)
  

  prop_df <- pairs_geo %>%
    group_by(distance_bin, Mechanism) %>%
    summarise(n = n(), .groups="drop") %>%
    group_by(distance_bin) %>%
    mutate(prop = n / sum(n))
  
  ggplot(prop_df,
         aes(x = distance_bin,
             y = prop,
             fill = Mechanism)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c(
      "Dispersal Limitation" = "#F0E442",
      "Drift" = "#56B4E9",
      "Homogeneous Selection" = "#D55E00",
      "Heterogeneous Selection" = "#D55E00",
      "Homogenising Dispersal" = "#F0E442"
    )) +
    theme_bw() +
    xlab("Geographic distance (km)") +
    ylab("Proportion of mechanisms") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  
  summary(as.vector(geo_dist))
  
  prop_df %>% 
    print(n = 24)
  