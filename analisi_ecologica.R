## =========================
## 0) Packages
## =========================
pkgs <- c("tidyverse","vegan","geodist","MicEco","picante","phyloseq","ggrepel")
to_install <- pkgs[!pkgs %in% installed.packages()[,1]]
if(length(to_install)) install.packages(to_install, dependencies = TRUE)
lapply(pkgs, library, character.only = TRUE)

## =========================
## 1) DATA: asv, meta, (tree)
## =========================
# EXPECTED:
# asv: data.frame/matrix, rows = samples, cols = ASVs, values = counts
# meta: data.frame with columns: sample_id, lat, lon, and environmental vars
# tree: (optional) phylo
# ---- Replace with your loads ----
# load("asv.RData"); load("meta.RData"); load("tree.RData")

stopifnot(all(rownames(asv) %in% meta$sample_id))
meta <- meta %>% filter(sample_id %in% rownames(asv))
asv  <- asv[meta$sample_id, , drop = FALSE]

# rimuovi colonne ASV con somma 0
asv <- asv[, colSums(asv) > 0, drop = FALSE]

## =========================
## 2) Trasformazioni standard
## =========================
# relative abundance e trasformazione Hellinger (per metodi lineari)
rel_abund <- sweep(asv, 1, rowSums(asv), FUN = "/")
rel_abund[is.na(rel_abund)] <- 0
hell <- decostand(rel_abund, method = "hellinger")

## =========================
## 3) Distanze
## =========================
# 3.1 Bray–Curtis (composizione)
d_bray <- vegdist(rel_abund, method = "bray")

# 3.2 UniFrac (opzionale, se hai l'albero)
use_unifrac <- exists("tree") && inherits(tree, "phylo")
if (use_unifrac) {
  ps <- phyloseq(otu_table(as.matrix(asv), taxa_are_rows = FALSE),
                 sample_data(meta %>% column_to_rownames("sample_id")),
                 phy_tree(tree))
  d_unifrac <- phyloseq::distance(ps, method = "unifrac", weighted = TRUE)
}

# 3.3 Distanza ambientale (scala variabili, Euclidea)
env_vars <- c("temp","sal","chl","o2","no3","po4","si")  # <-- adatta ai tuoi nomi
env_df <- meta %>% 
  select(all_of(c("sample_id", env_vars))) %>% 
  column_to_rownames("sample_id") %>% 
  mutate(across(everything(), as.numeric)) %>% 
  drop_na()
# match order
common_ids <- intersect(rownames(env_df), rownames(rel_abund))
env_df <- env_df[common_ids, , drop = FALSE]
rel_abund <- rel_abund[common_ids, , drop = FALSE]
hell <- hell[common_ids, , drop = FALSE]
d_bray <- vegdist(rel_abund, method = "bray")

d_env <- dist(scale(env_df), method = "euclidean")

# 3.4 Distanza geografica (km; great-circle)
geo_mat <- geodist(meta %>% select(sample_id, lon, lat),
                   paired = FALSE, measure = "geodesic")
rownames(geo_mat) <- meta$sample_id; colnames(geo_mat) <- meta$sample_id
geo_mat <- geo_mat[common_ids, common_ids]
d_geo <- as.dist(geo_mat/1000)  # km

## =========================
## 4) PERMANOVA (adonis2)
## =========================
# Formula con tutte le variabili ambientali disponibili (senza NA)
env_df2 <- env_df  # già allineato
form <- as.formula(paste("d_bray ~", paste(colnames(env_df2), collapse = " + ")))
set.seed(123)
perm <- adonis2(form, data = as.data.frame(env_df2), permutations = 999, by = "margin")
print(perm)

## =========================
## 5) db-RDA (CAP) e varpart
## =========================
# db-RDA (capscale su Bray-Curtis) con predittori ambientali
cap <- capscale(d_bray ~ ., data = as.data.frame(env_df2), add = TRUE) # add=TRUE per correzione Cailliez
anova_cap <- anova(cap, permutations = 999)
anova_terms <- anova(cap, by = "terms", permutations = 999)
print(anova_cap); print(anova_terms)

# Varianza fra Ambiente e Geografia: usa PCoA della distanza geografica come set di covariate
# (consigliato: pochi assi principali che spiegano >80% varianza)
pcoa_geo <- cmdscale(d_geo, k = max(2, min(5, nrow(env_df2)-1)), eig = TRUE)
G <- as.data.frame(pcoa_geo$points)
colnames(G) <- paste0("G", seq_len(ncol(G)))
E <- as.data.frame(scale(env_df2))
# varpart richiede la matrice (Hellinger) della comunità, non la distanza
vp <- varpart(hell, E, G)
print(vp)
# test frazioni (redundancy analysis)
rda_E <- rda(hell ~ ., data = E)
rda_G <- rda(hell ~ ., data = G)
rda_EG <- rda(hell ~ . + Condition(as.matrix(G)), data = E)
anova(rda_E); anova(rda_G); anova(rda_EG)

## =========================
## 6) Mantel e Partial Mantel
## =========================
set.seed(123)
mantel_env <- mantel(d_bray, d_env, permutations = 999)
mantel_geo <- mantel(d_bray, d_geo, permutations = 999)
mantel_part_env <- mantel.partial(d_bray, d_env, d_geo, permutations = 999)
mantel_part_geo <- mantel.partial(d_bray, d_geo, d_env, permutations = 999)
print(list(mantel_env=mantel_env,
           mantel_geo=mantel_geo,
           partial_env=mantel_part_env,
           partial_geo=mantel_part_geo))

## =========================
## 7) Distance–decay (Bray vs km)
## =========================
to_long <- function(d, name){
  dd <- as.matrix(d)
  tibble(
    i = rep(rownames(dd), times = ncol(dd)),
    j = rep(colnames(dd), each  = nrow(dd)),
    value = as.vector(dd)
  ) %>% 
    filter(i < j) %>% 
    rename(!!name := value)
}

df_bray <- to_long(d_bray, "bray")
df_km   <- to_long(d_geo, "km")
dd_df <- df_bray %>% inner_join(df_km, by = c("i","j"))

gg_distdecay <- ggplot(dd_df, aes(km, bray)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_x_continuous("Geographic distance (km)") +
  scale_y_continuous("Bray–Curtis dissimilarity") +
  ggtitle("Distance–decay of community similarity")
print(gg_distdecay)

# Modello (lineare e log-distanza)
m_lin <- lm(bray ~ km, data = dd_df)
m_log <- lm(bray ~ log1p(km), data = dd_df)
summary(m_lin); summary(m_log)

## =========================
## 8) Modello neutrale di Sloan (MicEco)
## =========================
# Richiede vettore di frequenze d'occurrence e abbondanze relative globali
# MicEco::sncm.fit implementa Sloan Neutral Community Model
set.seed(123)
sncm <- MicEco::sncm.fit(as.matrix(rel_abund))  # taxa in colonne, samples in righe
# Risultati: R2 del fit neutrale, m (immigration), devianza
print(sncm$R2)
print(sncm$fit$par)  # include m
# Classificazione taxa (sopra, dentro, sotto i bounds neutrale 95%)
neu_class <- MicEco::sncm.classify(sncm)
table(neu_class$group)
# Plot diagnostico
plot(sncm)

## =========================
## 9) (Opzionale) UniFrac-based analyses
## =========================
if (use_unifrac) {
  # Mantel con UniFrac
  mantel_unifrac_env <- mantel(as.dist(as.matrix(d_unifrac)), d_env, permutations = 999)
  mantel_unifrac_geo <- mantel(as.dist(as.matrix(d_unifrac)), d_geo, permutations = 999)
  print(list(mantel_unifrac_env=mantel_unifrac_env,
             mantel_unifrac_geo=mantel_unifrac_geo))
  
  # CAP con UniFrac: usa capscale direttamente sulla matrice distanza UniFrac
  cap_u <- capscale(as.dist(as.matrix(d_unifrac)) ~ ., data = as.data.frame(env_df2), add = TRUE)
  print(anova(cap_u, permutations = 999))
}

## =========================
## 10) Piccole utilità di report
## =========================
# Spiega quanta varianza spiega l'ambiente (PERMANOVA) e distance-decay
cat("\n=== SUMMARY HINTS ===\n")
cat(sprintf("- PERMANOVA R2(tot): %.3f (pseudo-F=%.2f, p=%.3f)\n",
            perm$R2[nrow(perm)-1], perm$F[nrow(perm)-1], perm$`Pr(>F)`[nrow(perm)-1]))
cat(sprintf("- Mantel (Bray~Env): r=%.3f, p=%.3f | Mantel (Bray~Geo): r=%.3f, p=%.3f\n",
            mantel_env$statistic, mantel_env$signif,
            mantel_geo$statistic, mantel_geo$signif))
cat(sprintf("- Partial Mantel (Env | Geo): r=%.3f, p=%.3f | (Geo | Env): r=%.3f, p=%.3f\n",
            mantel_part_env$statistic, mantel_part_env$signif,
            mantel_part_geo$statistic, mantel_part_geo$signif))
cat(sprintf("- Neutral model R2: %.3f; immigration m≈%.4f\n",
            sncm$R2, sncm$fit$par["m"]))
