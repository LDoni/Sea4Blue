library(igraph)
library(tidyverse)

source("R/Import_Data.R")
source("R/Datalist_Wrangling_Functions.R")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

datalist_Atlantic <- import_data("data/Atlantic/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso")) %>%
  mutate_meta_datalist(Ocean = "Atlantic") 

datalist_Pacific <- import_data("data/Pacific/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) %>%
  mutate_meta_datalist(Depth_Grp = ifelse(Depth <= DCM, "Epi", "Meso")) %>%
  mutate_meta_datalist(Ocean = "Pacific") 

datalist <- combine_data(mutate_meta_datalist(datalist_Atlantic, Station = as.character(Station)),
                         mutate_meta_datalist(datalist_Pacific, Station = as.character(Station)))

##
#network_pos <- read_graph("output/SparCC_Network_pos.txt", format = "graphml")

### o network_combined do 2_sparcc_analysis
network_pos <-network_combined
results <- tibble(Sample_ID = NA, Modularity = NA, Transitivity = NA, Mean_Dist = NA, Type = NA)

for (i in datalist$Meta_Data$Sample_ID) {
  
  tmp <- datalist %>%
    filter_station_datalist(Sample_ID == !!i)

  sub_network_pos <- subgraph(network_pos, which(names(V(network_pos)) %in% tmp$Count_Data$OTU_ID)) %>%
    delete_vertices(., which(degree(.) == 0))
  
  tmp_modul_pos <- cluster_walktrap(sub_network_pos, weights = NULL) %>%
    modularity()
  
  tmp_transit_pos <- transitivity(sub_network_pos, type = "average")
  
  tmp_dist_pos <- mean_distance(sub_network_pos)
  
  results <- results %>% 
    add_case(Sample_ID = !!i, Modularity = tmp_modul_pos, Transitivity = tmp_transit_pos, 
             Mean_Dist = tmp_dist_pos, Type = "Positive") 
}

data <- results %>%
  dplyr::slice(-1) %>%
  left_join(., datalist$Meta_Data, by = "Sample_ID")

write_csv(data, "output/Network_Properties.csv")

data <- read_csv("output/Network_Properties.csv")

data <- data %>%
  mutate(Longitude.x = as.numeric(Longitude.x))





data %>%
  filter(Type == "Positive") %>%
  ggplot(., aes(x = Longitude.x, y = Modularity, col = as.factor(Zona))) +
  geom_point(size=4) +
  geom_smooth(
    aes(group = 1),
    method = "loess",
    span = 0.5,
    se = FALSE,
    color = "black",
    linewidth = 1
  )+
  scale_color_manual(values = cbbPalette[c(4,7,3)]) +
  scale_x_continuous(
    breaks = c(-80, -60, -40, -20),
    labels = c("80°W", "60°W", "40°W", "20°W"),
    limits = c(-80, -5)
  )




data %>%
  filter(Type == "Positive", Transitivity > 0.4) %>%
  ggplot(aes(x = Longitude.x, y = Transitivity, col = as.factor(Zona))) +
  geom_point(size=4) +
  geom_smooth(
    aes(group = 1),
    method = "loess",
    span = 0.5,
    se = FALSE,
    color = "black",
    linewidth = 1
  )+
  scale_color_manual(values = cbbPalette[c(4,7,3)]) +
  scale_x_continuous(
    breaks = c(-80, -60, -40, -20),
    labels = c("80°W", "60°W", "40°W", "20°W"),
    limits = c(-80, -5)
  ) +
  labs(y = "Clustering coefficient") +
  theme_bw()











#cowplot::plot_grid(p2, p1, legend, nrow = 1, ncol = 3)

#ggsave("figs/network_properties_pos.pdf", width = 16, height = 5, dpi = 300)
