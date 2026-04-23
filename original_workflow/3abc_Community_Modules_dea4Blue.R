library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(viridis)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
source("R/Import_Data.R")
#source("../Jonatan_Project/R/Datalist_Wrangling_Functions.R")
source("Modular_Seascape/R/Import_SparCC_Network.R")
source("Modular_Seascape/R/Similarity_Indices.R")

#datalist <- import_data2("data/Atlantic/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) 


datalist_Atlantic <- import_data2("Data", kingdom = "Prok", abundance_filter = T, min_counts = 2000) 


datalist_cluster <- datalist_Atlantic

 
num_cols <- c("Latitude.x", "Longitude.x", "Latitude.y", "Longitude.y",
              "thetao", "so", "uo", "vo", "zos", "chl", "o2", "no3", "po4", "si", "Counts_Total")

# converti character -> numeric
datalist_cluster$Meta_Data[num_cols] <- lapply(datalist_cluster$Meta_Data[num_cols], as.numeric)


head(readRDS("output/SparCC/cluster.rds"))

datalist_cluster$Count_Data <-  datalist_cluster %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  .$Count_Data %>%
  left_join(., readRDS("output/SparCC/cluster.rds"), by = "OTU_ID") %>%
  filter(!is.na(Cluster)) %>%
  group_by(Cluster) %>%
  summarize_if(is.numeric, sum) %>%
  mutate(Cluster = as.character(Cluster)) %>%
  select(-n, -Degree)

n <- readRDS("output/SparCC/cluster.rds") %>%
  select(Cluster, n) %>%
  distinct() %>%
  group_by(Cluster) %>%
  summarize_all(sum) %>%
  arrange(Cluster) %>%
  .$n







plot_df <- datalist_cluster$Count_Data %>%
  pivot_longer(cols = -Cluster, names_to = "Sample_ID", values_to = "Abundance") %>%
  left_join(datalist_cluster$Meta_Data %>%
              select(Sample_ID,  Longitude = Longitude.x, Latitude=Latitude.x),
            by = "Sample_ID") %>%
  # 2️⃣ Calcola percentuale rispetto al totale del campione
  group_by(Sample_ID) %>%
  mutate(Abundance_pct = Abundance / sum(Abundance)) %>%
  ungroup()





colours <- readRDS("output/SparCC/cluster_w_coloros.rds")%>%
  mutate(Cluster = ordered(Cluster, levels = c(as.character(seq(1, length(unique(Cluster))-1)), "Others"))) %>%
  arrange(Cluster) %>%
  select(Cluster, Colour) %>%
  distinct() %>%
  .$Colour 


# 4️⃣ Plot area lungo la latitudine
ggplot(plot_df, aes(x = Longitude, y = Abundance  , fill = Cluster)) +
  geom_area(col = "black", size = 0.2) +
  scale_fill_manual(values = colours) +
  theme_test() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

#ggsave("figs/Network_Cluster_Abundance_FL_Complete.png", width = 8, height = 7, dpi = 300)

#### Get taxonomic composition of clusters with average abundances in dataset ####

average_asv_abundance <- datalist$Count_Data %>%
  mutate_if(is.numeric, function(x) x/sum(x)*100) %>%
  mutate(average = rowMeans(select_if(., is.numeric))) %>%
  select(OTU_ID, average)

datalist_cluster$Count_Data <- readRDS("output/SparCC/cluster_w_coloros.rds") %>%
  left_join(., average_asv_abundance, by = "OTU_ID") %>%
  select(OTU_ID, Cluster, average) %>%
  reshape2::dcast(OTU_ID~Cluster) %>%
  as_tibble() %>%
  mutate_if(is.numeric, function(x) ifelse(is.na(x), 0, x)) %>%
  left_join(., select_if(datalist$Count_Data, is.character), by = "OTU_ID") 





# 1. Leggi il file di tassonomia
taxonomy <- readr::read_tsv("Data/Taxonomy/Prok_Taxonomy.tsv")

datalist_clusterBackup<-datalist_cluster

datalist_cluster$Count_Data <- datalist_cluster$Count_Data %>%
  dplyr::left_join(taxonomy, by = "OTU_ID")

# Trasforma i moduli in formato long
long_modules <- datalist_cluster$Count_Data %>%
  pivot_longer(
    cols = c(`1`,`2`,`3`,`4`,`5`,`6`,`7`,`8`,`9`,`10`,`11`,`12`,`Others`),
    names_to = "Module",
    values_to = "Abundance"
  ) %>%
  # Combina Class e Family in una sola colonna
  mutate(Class_Family = str_c(Class, Family, sep = " - ")) %>%
  # opzionale: rimuovi "Others" se vuoi
  filter(Module != "Others")




otherThreshold <- 0.01

df_plot <- long_modules %>%
  dplyr::group_by(Class_Family) %>%
  dplyr::summarise(
    Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    rel_abundance = Abundance / sum(Abundance),
    Class_Family_thr = ifelse(rel_abundance < otherThreshold,
                              "Other",
                              Class_Family)
  ) %>%
  dplyr::select(Class_Family, Class_Family_thr) %>%
  right_join(long_modules, by = "Class_Family") %>%
  dplyr::group_by(Module, Class_Family_thr) %>%
  dplyr::summarise(
    Abundance = sum(Abundance),
    .groups = "drop"
  )

library(viridis)
library(forcats)

# livelli
taxa_levels <- levels(df_plot$Class_Family_thr)

# separo Other
taxa_main <- setdiff(taxa_levels, "Other")

# prendo 12 colori ben separati da Turbo
cols_main <- viridis(
  length(taxa_main),
  option = "turbo",
  begin = 0.05,
  end = 0.95
)

# mescolo per rompere il gradiente
set.seed(123)
cols_main <- sample(cols_main)

# palette finale
palette_cf <- c(
  setNames(cols_main, taxa_main),
  "Other" = "grey80"
)

ggplot(df_plot, aes(x = Module, y = Abundance, fill = Class_Family_thr)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette_cf) +
  theme_minimal() +
  labs(
    x = "Module",
    y = "Abundance",
    fill = "Class - Family"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.size = unit(0.6, "cm")
  )



ggplot(df_plot, aes(x = Module, y = Abundance, fill = Class_Family_thr)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d(option = "turbo", end = 0.95,) +
  theme_minimal() +
  labs(
    x = "Module",
    y = "Abundance",
    fill = "Class - Family"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )






### albero dei clusters

Tree <- ape::read.tree("Data/Trees/Prok/Prok_Combined.tree")

Count_Table <- datalist_cluster$Count_Data %>%
  select_if(is.numeric) %>%
  as.matrix() %>%
  magrittr::set_rownames(datalist_cluster$Count_Data$OTU_ID) %>%
  Matrix::Matrix()

my.ps <- phyloseq::phyloseq(phyloseq::otu_table(as.matrix(Count_Table), taxa_are_rows=T), 
                            phyloseq::phy_tree(Tree))

Unifrac_dist <- distance_wrapper(my.ps, method = "UniFrac_weighted", size.thresh = 1, pseudocount = 10^-6, nblocks = 100, 
                                 use.cores = 4, cor.use = "na.or.complete")




Unifrac_dist %>%
  as.matrix() %>%
  .[!rownames(.) %in% "Others",
    !colnames(.) %in% "Others"] %>%
  as.dist() %>%
  hclust(method = "mcquitty") %>%
  plot(hang = -1, axes = FALSE,
       main = "UniFrac clusters")


















### da qui non va

list_datatable <- create_datatable(datalist_cluster, grpBy = Family, upper_grp = Class, 
                                   otherThreshold = 0.01, addColorScheme = T)


list_datatable <- create_datatable(
  datalist_cluster, 
  grpBy = "Family", 
  upper_grp = "Class", 
  otherThreshold = 0.01, 
  addColorScheme = TRUE
)





datatable <- list_datatable$table %>%
  select(1:4) %>%
  filter(Sample_ID != "Others") %>%
  mutate(Sample_ID = as.numeric(Sample_ID)) %>%
  mutate(Group = plyr::mapvalues(Group, levels(Group), str_replace(levels(Group), pattern = ";", replacement = " - ")))

colorvalues <- list_datatable$color

Group_Levels <- str_replace_all(levels(datatable$Group), pattern = "Unknown Marinimicrobia \\(SAR406 clade\\) \\- Unknown Marinimicrobia \\(SAR406 clade\\)",
                replacement = "Unknown Marinimicrobia (SAR406 clade)")

datatable %>%
  mutate(Group = ordered(Group, levels = levels(datatable$Group), labels = Group_Levels)) %>%
  
  ggplot(., aes(x = as.factor(Sample_ID), y = Abundance, fill = Group)) +
    geom_bar(stat = "identity", position = "stack", width = .8, color = "black", size = .3) +
    scale_fill_manual(values = colorvalues) +
    labs(fill = "Class - Family", x = "Cluster", y = "Relative abundance within dataset (%)") +
    theme_bw() +
    theme(legend.position = "right",
          legend.text = element_text(size = 10),
          legend.key.size = unit(.8, "line")) +
    guides(fill = guide_legend(nrow = ceiling(length(unique(datatable$Group))), 
                               title.position = "top"))

ggsave("figs/Network_Cluster_Composition_rel_dataset.png", width = 10, height = 5.5, dpi = 300)
