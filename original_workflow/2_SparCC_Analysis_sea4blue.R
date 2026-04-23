library(igraph)
library(tidyverse)

source("Modular_Seascape/R/Import_Data.R")
source("Modular_Seascape/R/Datalist_Wrangling_Functions.R")
source("Modular_Seascape/R/Import_SparCC_Network.R")
source("Modular_Seascape/R/Similarity_Indices.R")
get_colors_cont <- function(vec, palette, reverse = F, n = 9) {
  
  vec <- ifelse(is.na(vec), 0, vec)
  
  if (reverse) { ramp <- rev(RColorBrewer::brewer.pal(n, palette)) %>% colorRamp(.) } else 
  { ramp <- RColorBrewer::brewer.pal(n, palette) %>%  colorRamp(.) } 
  
  ramp(vegan::decostand(vec, method = "range")) %>%
    rgb(., maxColorValue = 255)
}


library(plyr)
library(doParallel)
library(foreach) 
library(bigmemory)
library(Matrix)
library(biganalytics)
library(BiocParallel)

 
#datalist <- import_data2("Data_net/", kingdom = "Prok", abundance_filter = F, min_counts = 0) 
#usata questa intera:Count_Data
# A tibble: 1,223 × 16
datalist <- import_data2("Data/", kingdom = "Prok", abundance_filter = F, min_counts = 2000) 

source("Modular_Seascape/R/SparCC_Wrapper_sea4blue.R")

#result <- sparCC_wrapper(datalist, n_boot = 999, frac = F)
 
result <- sparCC_wrapper_with_pb(datalist, n_boot = 99, frac = F)

# Salva la matrice di correlazione con nomi di riga
write.csv(result$cor, "output/SparCC/Cor_SparCC_Prok_all.csv")

# Salva la matrice di p-value con nomi di riga
write.csv(result$pVal, "output/SparCC/Pval_SparCC_Prok_all.csv")
  

Max_Count <- datalist %>%
  mutate_count_datalist(function(x) x/sum(x)) %>%
  .$Count_Data %>%
  select_if(is.numeric) %>%
  mutate(Max = apply(., 1, which.max)) %>%
  mutate(OTU_ID = datalist$Count_Data$OTU_ID) %>%
  select(OTU_ID, Max) %>%
  cbind(., datalist$Meta_Data[.$Max,]) %>%
  as_tibble() %>%
  left_join(., select_if(datalist$Count_Data, is.character), by = "OTU_ID")

r_threshold = 0.51

files <- list.files("output/SparCC", pattern = "^Cor_SparCC_Prok_all.*", full.name = T)




data_combined <- map(files, function(x) {
  import_sparcc_network(cor_file = x, 
                        pval_file = str_replace_all(x, pattern = "Cor", replacement = "Pval"),
                        min_r = r_threshold, min_p = 0.05)
}) %>%
  bind_rows() %>%
  filter(Type == "Positive") %>%
  group_by(From, To) %>%
  summarize(weight = max(weight)) 

network_combined <- data_combined %>%
  select(-weight) %>%
  graph_from_data_frame(d = ., directed = F,
                        vertices = slice(Max_Count, match(unique(c(pull(., From), 
                                                                   pull(., To))),
                                                          OTU_ID)) %>%
                          mutate_if(is.factor, as.character) %>%
                          select_if(function(x) is.character(x) | is.numeric(x))) %>%
  set_vertex_attr(graph = ., name = "label", value = NA)

deg <- degree(network_combined, mode = "all")
V_size <- ifelse((log(deg)) < 3 , 2.5, (log(deg)*1.8))

layout_network <- layout_nicely(network_combined)

cluster <- tibble(Cluster = cluster_edge_betweenness(network_combined)$membership,
                  OTU_ID = V(network_combined)$name) %>%
  mutate(n = table(Cluster)[Cluster]) %>%
  mutate(Cluster = ifelse(n < 10, 99, Cluster)) %>%
  mutate(Cluster = c(seq(1, length(unique(Cluster))-1), "Others")[factor(Cluster, levels = unique(Cluster))]) %>%
  mutate(Degree = deg)


library(igraph)

write_graph(
  network_combined,
  file = "output/SparCC_Network_pos.graphml",
  format = "graphml"
)



unique(cluster$Cluster)
head(cluster)

#saveRDS(cluster, file = "output/SparCC/cluster.rds")
#write_csv(as.data.frame(cluster), "output/SparCC_Cluster.csv")#questo non funzia

#### Cluster validation ####
library(vegan)
library(BiocParallel)

count_y <- datalist$Count_Data %>%
  filter(OTU_ID %in% cluster$OTU_ID) %>%
  dplyr::slice(match(cluster$OTU_ID, OTU_ID)) %>%
  select_if(is.numeric)

adonis2(count_y ~ Cluster, data = cluster)
#adonis2(formula = count_y ~ Cluster, data = cluster) ----->>>>  fatto il 27/01/2026
#Df SumOfSqs      R2      F Pr(>F)    
#Model     12   114.85 0.39588 36.697  0.001 ***
### --> Probability that distances between groups defined by modules are obtained by chance: p < 0.001 
### (i.e. random assignment of ASVs to modules)

bplapply(seq_len(Nperm), function(x) {
  
  count_y <- datalist$Count_Data %>%
    filter(OTU_ID %in% cluster$OTU_ID) %>%
    dplyr::slice(match(cluster$OTU_ID, OTU_ID)) %>%
    select_if(is.numeric) %>%
    select_if(colSums(.) >= 8000) %>%
    t() %>%
    rrarefy(., sample = 8000) %>%
    t()
  
  perm <- adonis2(count_y ~ Cluster, data = cluster, by = "onedf")
  
})

modularity1 = cluster_edge_betweenness(network_combined)

obs <- modularity(modularity1)

Nperm = 1000

randomized.modularity <- bplapply(seq_len(Nperm), function(x){
  randomnet <- rewire(network_combined, with=each_edge(0.5))
  return(cluster_edge_betweenness(randomnet) %>% modularity())
})
head(randomized.modularity)

lapply(randomized.modularity, function(x) x > obs) %>%
  unlist() %>% sum()/Nperm
#[1] 0 ------->>>>> fatto il 27/01/2026

### --> Probability that a network with comparable modularity is formed: p < 0.001
### (i.e. the modularity of the graph partitioning)
### (using random rewiring -> 50% chance for an edge to rewire terminal node)



#### Network visualization ####

 
plot(network_combined,
     vertex.size = V_size,
     layout = layout_network * 0.03,
     vertex.color = get_colors_cont(as.numeric(vertex_attr(network_combined, "thetao")),
       "RdBu", reverse = TRUE),
     main = "Temperature")

#par(mfrow = c(2,3))

par(mfrow = c(2,3))

plot(network_combined,
     vertex.size = V_size,
     layout = layout_network * 0.03,
     vertex.color = get_colors_cont(as.numeric(vertex_attr(network_combined, "thetao")), "RdBu", reverse = TRUE),
     main = "Temperature")


plot(network_combined,
     vertex.size = V_size,
     layout = layout_network * 0.03,
     vertex.color = get_colors_cont(as.numeric(vertex_attr(network_combined,"Latitude.x")), "YlOrBr"),
     main = "Latitude")

plot(network_combined,
     vertex.size = V_size,
     layout = layout_network * 0.03,
     vertex.color = get_colors_cont(as.numeric(vertex_attr(network_combined,"Longitude.x")), "YlOrBr"),
     main = "Longitude") 


library(viridisLite)

plot(network_combined,
     vertex.size = V_size,
     layout = layout_network * 0.03,
     vertex.color = viridis(length(V(network_combined)))[
       rank(as.numeric(V(network_combined)$Longitude.x))
     ],
     main = "Longitude")




library(geosphere)

coords <- cbind(
  as.numeric(V(network_combined)$Longitude.x),
  as.numeric(V(network_combined)$Latitude.x)
)

dist_km <- distHaversine(coords,
                         coords[which.min(coords[,1]), ])

plot(network_combined,
     vertex.size = V_size,
     layout = layout_network * 0.03,
     vertex.color = get_colors_cont(dist_km, "YlOrBr"),  # palette valida
     main = "Distance from Florida (km)")




plot(network_combined,
     vertex.size = V_size,
     layout = layout_network * 0.03,
     vertex.color = get_colors_cont(as.numeric(vertex_attr(network_combined,"so")), "RdYlBu", reverse = TRUE),
     main = "Salinity")

plot(network_combined,
     vertex.size = V_size,
     layout = layout_network * 0.03,
     vertex.color =  get_colors_cont(as.numeric(vertex_attr(network_combined,"si")), "Reds"),
     main = "Silicate")

plot(network_combined,
     vertex.size = V_size,
     layout = layout_network * 0.03,
     vertex.color =  get_colors_cont(as.numeric(vertex_attr(network_combined,"chl")), "Greens"),
     main = "Chlorophyll")


plot(network_combined,
     vertex.size = V_size,
     layout = layout_network * 0.03,
     vertex.color = get_colors_cont(as.numeric(vertex_attr(network_combined,"Zona")), "RdYlBu", reverse = TRUE),
     main = "Zona")



plot(network_combined,
     vertex.size = V_size,
     layout = layout_network * 0.03,
     vertex.color = get_colors_cont(as.numeric(vertex_attr(network_combined,"uo")), "RdYlBu", reverse = TRUE),
     main = "uo ")
plot(network_combined,
     vertex.size = V_size,
     layout = layout_network * 0.03,
     vertex.color = get_colors_cont(as.numeric(vertex_attr(network_combined,"vo")), "RdYlBu", reverse = TRUE),
     main = "vo ")

plot(network_combined,
     vertex.size = V_size,
     layout = layout_network * 0.03,
     vertex.color = get_colors_cont(as.numeric(vertex_attr(network_combined,"o2")), "RdYlBu", reverse = TRUE),
     main = "o2")
plot(network_combined,
     vertex.size = V_size,
     layout = layout_network * 0.03,
     vertex.color = get_colors_cont(as.numeric(vertex_attr(network_combined,"po4")), "RdYlBu", reverse = TRUE),
     main = "po4")

plot(network_combined,
     vertex.size = V_size,
     layout = layout_network * 0.03,
     vertex.color = get_colors_cont(as.numeric(vertex_attr(network_combined,"si")), "RdYlBu", reverse = TRUE),
     main = "si")

plot(network_combined,
     vertex.size = V_size,
     layout = layout_network * 0.03,
     vertex.color = get_colors_cont(as.numeric(vertex_attr(network_combined,"Giorno.di.navigazione.")), "RdYlBu", reverse = TRUE),
     main = "Giorno.di.navigazione.")


V(network_combined)$current_speed <- with(
  as.data.frame(vertex.attributes(network_combined)),
  sqrt(as.numeric(uo)^2 + as.numeric(vo)^2)
)



plot(network_combined,
     vertex.size = V_size,
     layout = layout_network * 0.03,
     vertex.color = get_colors_cont(as.numeric(vertex_attr(network_combined,"current_speed")), "RdYlBu", reverse = TRUE),
     main = "Current speed (m/s)")
coords <- layout_network * 0.03



plot(
  network_combined,
  layout = coords,
  vertex.size = V_size,
  vertex.color = "grey80",
  edge.color = "grey70",
  main = "Surface currents"
)

arrows(
  x0 = coords[,1],
  y0 = coords[,2],
  x1 = coords[,1] + as.numeric(V(network_combined)$uo),
  y1 = coords[,2] + as.numeric(V(network_combined)$vo),
  length = 0.05,
  col = "dodgerblue",
  lwd = 1
)










library(igraph)
library(RColorBrewer)

# Estrai l'attributo Zona
zonas <- vertex_attr(network_combined, "Zona")

# Crea una palette discreta (una colore per ogni zona)
palette <- brewer.pal(n = length(unique(zonas)), name = "Set3") # Set3 va bene per categorie
colors <- setNames(palette, unique(zonas))
vertex_colors <- colors[zonas]

# Plot del network
plot(network_combined,
     vertex.size = V_size,
     layout = layout_network * 0.03,
     vertex.color = vertex_colors,
     main = "Zona")

legend("topright", legend = names(colors), col = colors, pch = 19, pt.cex = 1.5, bty = "n")


clusters_unique <- unique(cluster$Cluster)
n_clusters <- length(clusters_unique)

# Palette Set3 (12 colori) + aggiungiamo un colore per 'Others' se necessario
palette <- c(brewer.pal(12, "Set3"), "grey")
cluster_colors <- setNames(palette[1:n_clusters], clusters_unique)

#saveRDS(cluster, file = "output/SparCC/cluster_w_coloros.rds")



# Aggiungi la colonna colore all'oggetto cluster
cluster$Colour <- cluster_colors[cluster$Cluster]



unique()
# ---- 4️⃣ Usa i colori per il network ----
node_colors <- cluster$Colour[match(V(network_combined)$name, cluster$OTU_ID)]

plot(network_combined,
     vertex.size = V_size,
     rescale = TRUE,
     layout = layout_network * 0.03,
     vertex.color = node_colors,
     main = "Cluster")
legend(
  "topleft",
  legend = c(setdiff(clusters_unique, "Others"), "Others"),
  col = cluster_colors[c(setdiff(clusters_unique, "Others"), "Others")],
  pch = 19,
  pt.cex = 1.4,
  cex = 0.8,
  ncol = 2,
  bty = "n",
  title = "Cluster"
)
