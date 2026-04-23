# --- 0️⃣ Librerie ---
library(phyloseq)
library(igraph)
library(SpiecEasi)
library(tidyverse)
library(ggraph)
library(tidygraph)

# --- 1️⃣ Filtra taxa PRESENTI in più campioni ---
phy_gen <- tax_glom(phyloseq_obj_css, taxrank = "Genus")
# Opzione A: Filtra per prevalenza (presenza in % di campioni)
prev_threshold <- 0.1  # 10% dei campioni
phy_filt <- filter_taxa(phy_gen, function(x) sum(x > 0) > (prev_threshold * nsamples(phy_gen)), TRUE)

phy_filt


### da questo phyloseq_obj estraggo i dati che poi verrranno usati per l'analisi dei meccanismi ecologici di comunità
dir.create("Data_net", recursive = TRUE)
dir.create("Data_net/Meta_Data/Prok", recursive = TRUE)
dir.create("Data_net/Count_Data/Processed/Prok", recursive = TRUE)
dir.create("Data_net/Count_Data/Fasta/Prok", recursive = TRUE)
dir.create("Data_net/Trees/Prok", recursive = TRUE)
dir.create("Data_net/Processed/Prok", recursive = TRUE)
dir.create("Data_net/Taxonomy", recursive = TRUE)

#metadata
colnames(sample_data(phy_filt))[1] <- "Sample_ID"


sample_data(phy_filt)  %>%
  as_tibble()  %>%
  mutate(Sample_ID = as.character(Sample_ID)) %>%
  relocate(Sample_ID, .before = 1)%>%
  write_tsv("Data_net/Meta_Data/Prok/Meta_Data.tsv")

# OTU TABLE

otu<-otu_table(phy_filt) %>%
  as.data.frame() %>%
  rownames_to_column("OTU_ID") 
write_tsv(otu,"Data_net/Count_Data/Processed/Prok/Full_Prok_Count.tsv")  
write_tsv(otu,"Data_net/Processed/Prok/Full_Prok_Count.tsv")

#   Sequenze FASTA
writeXStringSet(refseq(phy_filt), "Data_net/Count_Data/Fasta/Prok/Full_Prok_Sequences.fasta")

# 4.4 Albero filogenetico
write.tree(phy_tree(phy_filt), "Data_net/Trees/Prok/Prok_Combined.tree")

# 4.5 Tassonomia
tax_table(phy_filt) %>%
  as.data.frame() %>%
  rownames_to_column("OTU_ID") %>%
  write_tsv("Data_net/Taxonomy/Prok_Taxonomy.tsv")














# Opzione B: Filtra per abbondanza minima in almeno 2 campioni
#phy_filt <- prune_taxa(
#  taxa_sums(phy_gen) > 10,  # o un valore sensato per i tuoi dati
#  phy_gen
#)

cat("Prima del filtraggio:", ntaxa(phy_gen), "generi\n")
cat("Dopo il filtraggio:", ntaxa(phy_filt), "generi\n")

# --- 2️⃣ Prepara matrice per SPIEC-EASI ---
# CSS normalization è già stata fatta (phyloseq_obj_css)
# Estrai e trasponi se necessario
otu_mat <- as(otu_table(phy_filt), "matrix")
if(taxa_are_rows(phy_filt)){
  otu_mat <- t(otu_mat)
}

# Verifica dimensioni
cat("Dimensioni matrice OTU:", dim(otu_mat), "\n")
cat("Numero di campioni:", nrow(otu_mat), "\n")
cat("Numero di generi:", ncol(otu_mat), "\n")

# --- 3️⃣ SPIEC-EASI con parametri adatti per dataset piccolo ---
# IMPORTANTE: Centra i dati prima di SPIEC-EASI
otu_mat_centered <- scale(otu_mat, center = TRUE, scale = FALSE)

se <- spiec.easi(
  otu_mat_centered, 
  method = 'mb',  # 'mb' (Meinshausen-Bühlmann) spesso meglio di 'glasso' per dati biologici
  lambda.min.ratio = 1e-3,  # Più basso per esplorare più lambda
  nlambda = 50,  # Più lambda
  pulsar.params = list(
    thresh = 0.05,
    subsample.ratio = 0.8,  # Importante per pochi campioni
    rep.num = 20  # Più repliche per stabilità
  ),
  icov.select.params = list(rep.num = 50)  # Più stabilità nella selezione
)










# --- 4️⃣ Estrai rete e verifica ---
adj <- getRefit(se)

# Verifica dimensioni
cat("Dimensioni matrice di adiacenza:", dim(adj), "\n")
cat("Numero di nodi:", nrow(adj), "\n")

# Controlla se corrisponde al numero di generi
if(nrow(adj) != ncol(otu_mat)){
  cat("ATTENZIONE: Dimensioni non corrispondono!\n")
  cat("Probabilmente alcuni taxa sono stati rimossi da SPIEC-EASI\n")
  
  # Estrai i taxa mantenuti
  kept_taxa <- colnames(otu_mat)[1:nrow(adj)]  # SPIEC-EASI mantiene l'ordine
  V(g)$name <- kept_taxa
} else {
  V(g)$name <- colnames(otu_mat)
}
head(sample_data(phy_filt))
phy_filt
# Crea grafo
g <- adj2igraph(adj)  # Funzione specifica di SpiecEasi

# Versione minimale per il tuo codice originale
tax_table_df <- as.data.frame(tax_table(phy_filt))

# Crea nomi per i nodi
node_names <- ifelse(
  is.na(tax_table_df$Genus) | tax_table_df$Genus == "",
  tax_table_df$Family,
  tax_table_df$Genus
)

# Assicurati l'ordine corretto
ordered_indices <- match(colnames(otu_mat), rownames(tax_table_df))
V(g)$name <- node_names[ordered_indices]

# Ora visualizza
ggraph(g, layout = "fr") +
  geom_edge_link() +
  geom_node_point(aes(color = factor(module)), size = 5) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3,check_overlap = Inf) +
  theme_void() +
  ggtitle("Co-occurrence network - Genus level")

ggraph(g, layout="fr") +
  geom_edge_link() +
  geom_node_point(aes(color=factor(module)), size=5) +
  geom_node_text(aes(label=name), repel=TRUE, size=3) +
  theme_void() +
  ggtitle("Co-occurrence network - Genus level")












# --- 5️⃣ ANALIZZA la rete ---
# Verifica se la rete è vuota
if(ecount(g) == 0){
  cat("ATTENZIONE: La rete è VUOTA! Nessuna connessione trovata.\n")
  cat("Possibili cause:\n")
  cat("1. Troppi taxa per pochi campioni\n")
  cat("2. Parametri lambda troppo conservativi\n")
  cat("3. Assenza di reali pattern di co-occorrenza\n")
}

# Statistiche della rete
cat("\n=== STATISTICHE RETE ===\n")
cat("Nodi (generi):", vcount(g), "\n")
cat("Connessioni (archi):", ecount(g), "\n")
cat("Densità:", edge_density(g), "\n")
cat("Grado medio:", mean(degree(g)), "\n")

# --- 6️⃣ Modularità (communities) ---
# Solo se ci sono archi
if(ecount(g) > 0){
  # Usa cluster_louvain per modularità
  comm <- cluster_louvain(g)
  V(g)$module <- membership(comm)
  
  cat("Numero di moduli:", length(unique(V(g)$module)), "\n")
  cat("Modularità:", modularity(comm), "\n")
} else {
  V(g)$module <- 1  # Tutti nello stesso modulo se non ci sono archi
}

# --- 7️⃣ Visualizza SOLO se la rete non è vuota ---
if(ecount(g) > 0 && vcount(g) > 1){
  # Versione migliorata
  deg <- degree(g)
  deg_thr <- quantile(deg, 0.75)
  V(g)$degree <- deg
  V(g)$label <- ifelse(deg > deg_thr, V(g)$name, "")
  
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(alpha = 0.3, width = 0.5) +
    geom_node_point(aes(color = as.factor(module)), size = 4) +
    geom_node_text(
      aes(label = label),
      repel = TRUE,
      size = 2.5,
      max.overlaps = 20
    ) +
    theme_void() +
    labs(
      title = paste("Co-occurrence network -", vcount(g), "generi"),
      subtitle = paste(
        ecount(g),
        "connessioni | Modularità:",
        round(modularity(comm), 3)
      )
    )
  
  print(p)
  
  # Salva plot
  ggsave("network_plot.png", p, width = 10, height = 8)
} else {
  cat("\nNon è possibile visualizzare una rete vuota o con un solo nodo.\n")
  
  # Alternative: mostra heatmap delle correlazioni
  cor_mat <- cor(otu_mat_centered)
  heatmap(cor_mat, main = "Matrice di correlazione (alternativa)")
}

# --- 8️⃣ Salva risultati dettagliati ---
results <- list(
  network = g,
  adjacency_matrix = adj,
  spiec_easi_object = se,
  otu_matrix = otu_mat,
  filtered_taxa = taxa_names(phy_filt),
  network_stats = c(
    nodes = vcount(g),
    edges = ecount(g),
    density = graph.density(g)
  )
)

saveRDS(results, "cooccurrence_network_results.rds")

# --- 9️⃣ ESPORTA TABELLA per analisi esterna ---
if(ecount(g) > 0){
  # Crea edge list con nomi dei generi
  edges <- as_data_frame(g, what = "edges")
  
  # Aggiungi tassonomia se disponibile
  tax_info <- as.data.frame(tax_table(phy_filt))
  tax_info$Genus <- rownames(tax_info)
  
  edges <- edges %>%
    left_join(tax_info, by = c("from" = "Genus")) %>%
    left_join(tax_info, by = c("to" = "Genus"), suffix = c("_from", "_to"))
  
  write.csv(edges, "network_edges_with_taxonomy.csv", row.names = FALSE)
  cat("\nEdge list esportata in 'network_edges_with_taxonomy.csv'\n")
}



table(V(g)$module)
length(unique(V(g)$module))


### --->>>>i moduli spiegano la composizione?


### metto i moduli nella tassonomia

tax_df <- as.data.frame(tax_table(phy_filt))
tax_df$module <- V(g)$module
tax_table(phy_filt) <- tax_table(as.matrix(tax_df))


transitivity(g, type = "global")   # clustering coefficient
diameter(g)

comm <- cluster_louvain(g)
modularity(comm)
length(unique(membership(comm)))



#Test randomizzazione rete (fortemente consigliato)

obs_mod <- modularity(comm)

rand_mod <- replicate(999, {
  g_rand <- rewire(g, with = keeping_degseq(niter = ecount(g) * 10))
  modularity(cluster_louvain(g_rand))
})

p_val <- mean(rand_mod >= obs_mod)
p_val







# PERMANOVA: i moduli spiegano la composizione?
otu <- (otu_table(phy_filt))
head(otu)
mod <- tax_table(phy_filt)[, "module"]

adonis2(
  otu ~ mod,
  method = "bray",
  permutations = 9999
)

Permutation test for adonis under reduced model
Permutation: free
Number of permutations: 9999

adonis2(formula = otu ~ mod, permutations = 9999, method = "bray")
Df SumOfSqs      R2      F Pr(>F)    
Model     10   15.700 0.32567 6.6165  1e-04 ***
  Residual 137   32.508 0.67433                  
Total    147   48.208 1.00000     



#Abbondanza per modulo
library(phyloseq)
library(tidyverse)
library(microbiome)
library(ggplot2)

# --- 1️⃣ Trasforma in abbondanza relativa ---
phy_rel <- transform_sample_counts(phy_filt, function(x) x / sum(x))

# --- 2️⃣ Estrai OTU table e taxa ---
otu_rel <- as.data.frame(t(otu_table(phy_rel)))  # ASV come colonne
otu_rel$Sample <- rownames(otu_rel)

# --- 3️⃣ Aggiungi moduli e Genus reale dai tax_table ---
tax_df <- as.data.frame(tax_table(phy_rel))
tax_df$module <- V(g)$module[match(rownames(tax_df), taxa_names(phy_rel))]
tax_df$Genus_real <- tax_df$Genus  # se vuoi il nome reale del Genus

head(tax_df)

# --- 4️⃣ Ottieni info campione ---
sample_df <- meta(phy_rel)
sample_df$Sample <- rownames(sample_df)
# rinomina la colonna senza punto finale se presente
if("Giorno.di.navigazione." %in% colnames(sample_df)){
  sample_df <- sample_df %>%
    rename(Giorno.di.navigazione = Giorno.di.navigazione.)
}
head(sample_df)

tax_df2 <- tax_df %>%
  rownames_to_column(var = "ASV")

# Pivot lungo e join corretto
df_long <- otu_rel %>%
  pivot_longer(
    cols = -Sample,
    names_to = "ASV",
    values_to = "Abundance"
  ) %>%
  # join con tax_table
  left_join(tax_df2[, c("ASV", "module", "Genus_real")], by = "ASV") %>%
  # join con info campione
  left_join(sample_df[, c("Sample", "Giorno.di.navigazione","Zona","Latitude.x", "Longitude.x","thetao")], by = "Sample") %>%
  mutate(Giorno.di.navigazione = as.numeric(Giorno.di.navigazione)) %>%
  arrange(Giorno.di.navigazione)





library(RColorBrewer)

ggplot(df_long, aes(x = factor(Sample ), y = Abundance, fill = factor(module))) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ factor(Zona, levels = c("ANW", "ANC", "ANE")), 
             scales = "free_x") +
  theme_bw() +
  labs(
    x = "Giorno di navigazione",
    y = "Abbondanza relativa",
    fill = "Modulo"
  ) +
  scale_fill_brewer(palette = "Set3") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(1, "lines")
  )



## --->>>>   Moduli ↔ ambiente

library(ggplot2)

# Raggruppiamo per modulo e campione
df_env <- df_long %>%
  group_by(Sample, module, Giorno.di.navigazione) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  left_join(sample_df[, c("Sample", "thetao")], by = "Sample")

ggplot(df_env, aes(x = thetao, y = Abundance, color = factor(module))) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  labs(
    x = "Thetao (°C)",
    y = "Abbondanza relativa del modulo",
    color = "Modulo"
  ) +
  scale_color_viridis(discrete = TRUE, option = "D")


## --->>>>   Moduli ↔ filogenesi
library(ggtree)
library(ape)

# Supponiamo tu abbia l'albero phylo dei tuoi ASV
phylo_tree <- phy_tree(phy_filt)  # phyloseq tree
#phylo_tree <- as.phylo(phy_filt)   # se necessario

# Costruisci un dataframe per mappare i colori dei moduli
tip_data <- data.frame(
  ASV = phylo_tree$tip.label,
  module = tax_df$module[match(phylo_tree$tip.label, rownames(tax_df))],
  Genus_real = tax_df$Genus_real[match(phylo_tree$tip.label, rownames(tax_df))]
)

# Plot albero con colori moduli
p <- ggtree(phylo_tree) %<+% tip_data +
  geom_tippoint(aes(color = factor(module)), size = 3) +
  scale_color_viridis(discrete = TRUE, option = "D") +
  theme_tree2() +
  labs(color = "Modulo")

p






















