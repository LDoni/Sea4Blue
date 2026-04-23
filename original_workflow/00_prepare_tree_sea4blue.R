#### prepare Fasta for SINA alignment ####

library(tidyverse)

source("Modular_Seascape/R/Import_Data.R")
source("Modular_Seascape/R/Datalist_Wrangling_Functions.R")
source("Modular_Seascape/R/Import_SparCC_Network.R")

datalist_Atlantic <- import_data2("Data/", kingdom = "Prok", abundance_filter = T, min_counts = 2000) 

prok_fasta_Atlantic <- seqinr::read.fasta("Data/Count_Data/Fasta/Prok/Full_Prok_Sequences.fasta")


#https://www.arb-silva.de/aligner/


# 1. Convert all sequences to safe format first
safe_sequences <- lapply(prok_fasta_Atlantic, function(x) {
  if(is.null(x) || !is.character(x)) return(NA)
  paste(x, collapse = "")
})

# 2. Remove any invalid sequences
valid_seqs <- !is.na(safe_sequences)
prok_fasta_clean <- safe_sequences[valid_seqs]

# 3. Verify matching with count data
matched_inds <- match(datalist_Atlantic$Count_Data$OTU_ID, names(prok_fasta_clean))
matched_seqs <- prok_fasta_clean[na.omit(matched_inds)]

# 4. Write in chunks with proper error handling
write_fasta_safely <- function(seq_list, chunk_size = 612, prefix = "Prok_Atlantic") {
  chunks <- split(seq_list, ceiling(seq_along(seq_list)/chunk_size))
  
  for(i in seq_along(chunks)) {
    out_file <- paste0("Data/Count_Data/Fasta/Prok/", prefix, "_part", i, ".fasta")
    
    # Convert to raw character vectors
    seqs_to_write <- unlist(lapply(chunks[[i]], function(x) {
      if(is.list(x)) return(paste(unlist(x), collapse = ""))
      x
    }))
    
    seqinr::write.fasta(
      sequences = as.list(seqs_to_write),  # Must be a list of character vectors
      names = names(chunks[[i]]),
      file.out = out_file,
      as.string = TRUE  # Critical parameter
    )
    message("Successfully wrote ", length(chunks[[i]]), " sequences to ", out_file)
  }
}

# 5. Execute
write_fasta_safely(matched_seqs)


### Then use online-version of SINA aligner in standard settings ####

#### Create Tree with FastTree ####

system(paste0("~/PhD/Statistics/FastTree/FastTree -gtr -nt < ",
              "data/Atlantic/Aligned/Prok/SINA_Aligned_Prok_Atlantic.fasta > ",
              "data/Atlantic/Tree/Prok/Atlantic_Prok_FastTree.tree"))


###  Alla fine ho già il tree fatto con dada2
