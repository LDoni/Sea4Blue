library(umap)
library(ggplot2)
library(ggrepel)
library(dplyr)
### UMAP!!!


bray <- phyloseq::distance(phyloseq_obj_css, method="bray")
umap_res <- umap(as.matrix(bray))

zones <- sample_data(phyloseq_obj_css)$Zone

# assegno colori alle 3 zone
zone_colors <- c("ANW" = "#4280fc", "ANC" = "#ffb452", "ANE" = "#f7170a")

umap_df <- as.data.frame(umap_res$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")

# aggiungi metadata
umap_df$Sample <- rownames(umap_df)
meta_df <- as.data.frame(sample_data(phyloseq_obj_css))
meta_df$Sample <- rownames(meta_df)

umap_df <- umap_df %>%
  left_join(meta_df, by = "Sample")

library(ggplot2)
library(ggrepel)
library(rlang)
# plot con ggplot
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Zona, label = Sample )) +
  geom_point(size = 3, alpha = 0.8) +
  ggrepel::geom_text_repel(max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = c("ANW" = "#4280fc",
                                "ANC" = "#ffb452",
                                "ANE" = "#f7170a")) +
  theme_bw() +
  labs(
    title = "UMAP su distanza Bray-Curtis",
    color = "Zone"
  )


# 1. Calcolo UMAP sui Bray-Curtis
bray <- phyloseq::distance(phyloseq_obj_css, method = "bray")
umap_res <- umap(as.matrix(bray))

# 2. Costruisco un data.frame per ggplot
umap_df <- as.data.frame(umap_res$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")

# aggiungo metadata
umap_df$Sample <- rownames(umap_df)
meta <- as.data.frame(sample_data(phyloseq_obj_css))
meta$Sample <- rownames(meta)

umap_df <- umap_df %>%
  left_join(meta, by = "Sample")

# 3. Plot con ggplot
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Zona)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = Sample), show.legend = FALSE) +
  scale_color_manual(values = c("ANW"="#4280fc", "ANC"="#ffb452", "ANE"="#f7170a")) +
  theme_bw() +
  labs(
    title = "UMAP su distanza Bray-Curtis",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Zone"
  )





library(Rtsne)
library(phyloseq)

bray <- phyloseq::distance(phyloseq_obj_css, method = "bray")

set.seed(123)
tsne_res <- Rtsne(
  as.matrix(bray),
  is_distance = TRUE,
  perplexity = 3,   # prova anche 4
  theta = 0.5
)
library(ggplot2)
library(ggrepel)
library(dplyr)

tsne_df <- data.frame(
  TSNE1 = tsne_res$Y[,1],
  TSNE2 = tsne_res$Y[,2],
  Sample = rownames(as.matrix(bray))
)

# aggiungi metadata
meta <- as.data.frame(sample_data(phyloseq_obj_css))
meta$Sample <- rownames(meta)

tsne_df <- tsne_df %>%
  left_join(meta, by = "Sample")

# colori coerenti con UMAP/PCoA
zone_colors <- c(
  "ANW"="#4280fc",
  "ANC"="#ffb452",
  "ANE"="#f7170a"
)

ggplot(tsne_df, aes(x = TSNE1, y = TSNE2, color = Zona)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(aes(label = Sample), show.legend = FALSE, max.overlaps = Inf) +
  scale_color_manual(values = zone_colors) +
  theme_bw() +
  labs(
    title = "t-SNE (Bray–Curtis distance)",
    x = "t-SNE 1",
    y = "t-SNE 2",
    color = "Zone"
  )
