library(phyloseq)
library(dplyr)

filter_abundance_phyloseq <- function(physeq) {
  
  # Estraggo la matrice OTU (taxa x samples)
  otu_mat <- as.data.frame(otu_table(physeq))
  if(taxa_are_rows(physeq) == FALSE) {
    otu_mat <- t(otu_mat)
  }
  
  # Calcolo le proporzioni per ogni campione (colonna)
  prop <- sweep(otu_mat, 2, colSums(otu_mat), FUN = "/")
  
  # 1) Filtro taxa con somma relativa su tutti i campioni > 0.00001
  filter1 <- (rowSums(otu_mat) / sum(otu_mat)) > 0.00001
  
  # 2) Filtro taxa per max proporzione in almeno un campione > 0.01
  filter2 <- apply(prop, 1, max) > 0.01
  
  # 3) Filtro taxa per taxa con proporzione > 0.001 in almeno il 2% dei campioni
  filter3 <- (rowSums(prop > 0.001) / ncol(prop)) > 0.02
  
  # 4) Filtro taxa con presenza (prop > 0) in almeno il 5% dei campioni
  filter4 <- (rowSums(prop > 0) / ncol(prop)) > 0.05
  
  # Unisco i filtri con la stessa logica OR / AND del tuo codice:
  # (filter1 & filter2) | filter3 | filter4
  keep_taxa <- (filter1 & filter2) | filter3 | filter4
  
  # Filtra l'oggetto phyloseq mantenendo solo taxa che passano il filtro
  physeq_filtered <- prune_taxa(keep_taxa, physeq)
  
  return(physeq_filtered)
}


library(furrr)
source("Functions_Similarity_IndicesExCom.R")
source("Similarity_IndicesExCom.R")
ps <- readRDS("ps.rds")
otu_table(ps) <- otu_table(t(otu_table(ps)), taxa_are_rows = TRUE)
head(tax_table(ps))

phyloseq_obj<-ps
grep(pattern = "Mitochondria", tax_table(phyloseq_obj)) 
grep(pattern = "Chloroplast", tax_table(phyloseq_obj)) 
phyloseq_obj <- phyloseq_obj %>% subset_taxa( Family!= "Mitochondria" | is.na(Family) & Class!="Chloroplast" | is.na(Class) ) 

phyloseq_obj <- subset_taxa(phyloseq_obj, (tax_table(phyloseq_obj)[,"Order"]!="Chloroplast") | is.na(tax_table(phyloseq_obj)[,"Order"]))

phyloseq_obj

#filtering singletons
doubleton <- genefilter_sample(phyloseq_obj, filterfun_sample(function(x) x > 1), A=1)
doubleton <- prune_taxa(doubleton, phyloseq_obj) 

doubleton = subset_samples(doubleton,!( Nominativo.campione.=="Positive" | Nominativo.campione.=="Negative"))

sample_names(doubleton)<-sample_data(doubleton)$Nominativo.campione. 
length(sample_names(doubleton))


doubleton_filtered<-filter_abundance_phyloseq(doubleton) 



dist_matrix <-distance_wrapper(
  phyloseq_obj = t(otu_table(doubleton_filtered)),  # Nota il t() per trasporre
  method = "PINA_unweighted",
  use.cores = 4)


shared_asvs <- intersect(
  taxa_names(doubleton_filtered),
  phy_tree(doubleton_filtered)$tip.label
)

length(shared_asvs)

 distance_wrapper( phyloseq_obj = doubleton_filtered,  method = "PINA_weighted",  use.cores = 8)

 
 
 sample_names(doubleton_filtered)
 # Correggi l'orientamento (taxa come righe, campioni come colonne)

 head(otu_table(doubleton_filtered))






















distance_wrapper(phyloseq_obj = doubleton_filtered,
                 method = "bray_curtis",
                 use.cores = 4) 

distance_wrapper(phyloseq_obj = doubleton,
                 method = "TINA_weighted",
                 use.cores = 4)  # puoi aumentare i core se vuoi


distance_wrapper(phyloseq_obj = doubleton,
                 method = "TINA_weighted",
                 size.thresh = 1,  # invece di 5 o 10
                 use.cores = 4,nblocks = 100,pseudocount = 0)


SparCC_wrapper <- function(phyloseq_obj,
                           size.thresh = 5,
                           use.cores = 2,
                           pseudo = 1e-6,
                           nblocks = NULL) {
  
  Count_Table <- as.data.frame(otu_table(phyloseq_obj))
  if (!taxa_are_rows(phyloseq_obj)) {
    Count_Table <- t(Count_Table)
  }
  
  # Filtro
  Count_Table_Filtered <- Count_Table[rowSums(Count_Table) > size.thresh, ]
  n_taxa <- nrow(Count_Table_Filtered)
  
  if (n_taxa < 2) {
    stop("Not enough taxa after filtering. Consider reducing size.thresh.")
  }
  
  # Definisci nblocks
  if (is.null(nblocks)) {
    nblocks <- max(1, min(floor(n_taxa / use.cores), n_taxa))
  }
  if (nblocks > n_taxa) {
    nblocks <- n_taxa
  }
  
  message("Using ", nblocks, " blocks across ", n_taxa, " taxa.")
  
  Count_Table_Filtered <- Count_Table_Filtered + pseudo
  
  cor.est <- sparcc(Count_Table_Filtered,
                    iter = 20,
                    inner.iter = 10,
                    th = 0.1,
                    exclusion.threshold = 0.1,
                    nblocks = nblocks)
  
  return(cor.est$Cor)
}
