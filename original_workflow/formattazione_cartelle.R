### I CHL E MITOK SONO STATI RIMOSSI PRIMA DI FARE TUTTO





## formattazione caretelle
#prima era sbagliato (ancora da correggere nello script)
  mkdir -p Data/Processed/Prok/
  mkdir -p Data/Meta_Data/Prok/
  cp Data/Count_Data/Processed/Prok/Full_Prok_Count.tsv Data/Processed/Prok/
  cp Data/Meta_Data/Prok/Meta_Data.tsv Data/Meta_Data/Prok/  # (already correct)
  
    
    
    ps <- readRDS("ps.rds")
    otu_table(ps) <- otu_table(t(otu_table(ps)), taxa_are_rows = TRUE)
    head(tax_table(ps))
    
    phyloseq_obj<-ps
    grep(pattern = "Mitochondria", tax_table(phyloseq_obj)) 
    grep(pattern = "Chloroplast", tax_table(phyloseq_obj)) 
    phyloseq_obj <- phyloseq_obj %>% subset_taxa( Family!= "Mitochondria" | is.na(Family) & Class!="Chloroplast" | is.na(Class) ) 
    
    phyloseq_obj <- subset_taxa(phyloseq_obj, (tax_table(phyloseq_obj)[,"Order"]!="Chloroplast") | is.na(tax_table(phyloseq_obj)[,"Order"]))
    
    phyloseq_obj  
    
    
  sample_data(ps) %>%
  subset_samples(!(sample.description %in% c("Laposnegative", "Lapospositive"))) %>%
  as_tibble() %>%
  rename(Sample_ID = Nominativo.campione.) %>%
  mutate(Sample_ID = as.character(Sample_ID)) %>%
  relocate(Sample_ID, .before = 1)%>%
  write_tsv("Data/Meta_Data/Prok/Meta_Data.tsv")

library(phyloseq)
library(tidyverse)
library(vegan)
library(picante)

# 1. CREAZIONE STRUTTURA CARTELLE ----
dir.create("Data", recursive = TRUE)
dir.create("Data/Meta_Data/Prok", recursive = TRUE)
dir.create("Data/Count_Data/Processed/Prok", recursive = TRUE)
dir.create("Data/Count_Data/Fasta/Prok", recursive = TRUE)

# 2. ESPORTAZIONE DATI PHYLOSEQ ----
# 2.1 File dei conteggi processati
otu_table(ps) %>%
  as.data.frame() %>%
  rownames_to_column("OTU_ID") %>%
  write_tsv("Data/Count_Data/Processed/Prok/Full_Prok_Count.tsv") %>%
  write_tsv("Data/Processed/Prok/Full_Prok_Count.tsv")
#a volte serve in questa cartella
otu_table(ps) %>%
  as.data.frame() %>%
  rownames_to_column("OTU_ID") %>%
  write_tsv("Data/Processed/Prok/Full_Prok_Count.tsv")


# 2.2 Metadati
sample_data(ps) %>%
  as.data.frame() %>%
  rownames_to_column("Sample_ID") %>%
  write_tsv("Data/Meta_Data/Prok/Meta_Data.tsv")

# 2.3 Sequenze FASTA
writeXStringSet(refseq(ps), "Data/Count_Data/Fasta/Prok/Full_Prok_Sequences.fasta")

write tree... "Data/Trees/Prok/Prok_Combined.tree"

write taxonomy: Data/Taxonomy

