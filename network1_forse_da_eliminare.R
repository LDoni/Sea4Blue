# --- 0️⃣ Librerie ---
library(phyloseq)
library(igraph)
library(SpiecEasi)
library(tidyverse)
library(ggraph)

# --- 1️⃣ Filtra taxa abbondanti ---
# Ad esempio genera presenti almeno al 5% in almeno 2 campioni
phy_gen <- tax_glom(phyloseq_obj_css, taxrank = "Genus")
phy_gen1 <- prune_taxa(taxa_sums(phy_gen) > 0.005 * sum(taxa_sums(phy_gen)), phy_gen)

# --- 2️⃣ Estrai matrice OTU normalizzata ---
otu_mat <- as.matrix(otu_table(phy_gen1))
if(taxa_are_rows(phy_gen1)){
  otu_mat <- t(otu_mat)
}

# --- 3️⃣ SPIEC-EASI network inference ---
# Metodo: glasso (sparse inverse covariance)
se <- spiec.easi(otu_mat, method='glasso', lambda.min.ratio=1e-2, nlambda=20, pulsar.params=list(thresh=0.05))

# --- 4️⃣ Estrai rete ---
adj <- as.matrix(getRefit(se))
g <- graph_from_adjacency_matrix(adj, mode="undirected", diag=FALSE)

# Aggiungi nome genera ai nodi
V(g)$name <- taxa_names(phy_gen1)


colnames(adj)  



# --- 6️⃣ Visualizza rete ---
ggraph(g, layout="fr") +
  geom_edge_link() +
  geom_node_point(aes(color=factor(module)), size=5) +
  geom_node_text(aes(label=name), repel=TRUE, size=3) +
  theme_void() +
  ggtitle("Co-occurrence network - Genus level")

# --- 7️⃣ Salva risultati ---
saveRDS(g, "cooccurrence_network_genus.rds")
