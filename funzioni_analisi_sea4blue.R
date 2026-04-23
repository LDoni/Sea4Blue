

#import data
import_data2 <- function(file_ASV, kingdom = "Prok", rare_lim = NULL, 
                         drop_rare = TRUE, abundance_filter = FALSE, min_counts = NULL) {
  
  # 1. Funzione per leggere i dati
  data_select <- function(file_ASV, kingdom) {
    # Percorsi file (corretti)
    count_path <- file.path(file_ASV, "Count_Data", "Processed", kingdom, paste0("Full_", kingdom, "_Count.tsv"))
    meta_path <- file.path(file_ASV, "Meta_Data", kingdom, "Meta_Data.tsv")
    
    # Lettura dati con controllo errori
    Count_Data <- tryCatch({
      suppressMessages(
        readr::read_delim(
          count_path,
          delim = "\t",
          col_types = readr::cols(
            OTU_ID = readr::col_character(),
            .default = readr::col_double()
          ),
          na = c("", "NA", "N/A", "#N/A", "NaN")
        )
      ) %>%
        dplyr::select(-dplyr::matches("[Mm]ock|NC"))
    }, error = function(e) {
      stop("Errore lettura count data: ", e$message)
    })
    
    Meta_Data <- tryCatch({
      suppressMessages(
        readr::read_delim(
          meta_path,
          delim = "\t",
          col_types = readr::cols(.default = readr::col_character())
        )
      )
    }, error = function(e) {
      stop("Errore lettura metadata: ", e$message)
    })
    
    # Verifica corrispondenza campioni (versione semplificata)
    count_samples <- names(Count_Data)[-1]
    meta_samples <- Meta_Data$Sample_ID
    
    if(!all(count_samples %in% meta_samples)) {
      missing_samples <- setdiff(count_samples, meta_samples)
      warning(length(missing_samples), " campioni nei dati di conteggio non presenti nei metadati: ", 
              paste(missing_samples, collapse = ", "))
    }
    
    return(list(Meta_Data = Meta_Data, Count_Data = Count_Data))
  }
  
  # 2. Funzione correzione tassonomia (semplificata)
  correct_ambiguous <- function(datalist) {
    counts <- datalist$Count_Data
    tax_cols <- names(counts)[sapply(counts, is.character)]
    
    patterns <- c('bacterium$', "uncultured", "metagenome", "unidentified", 
                  "Ambiguous_taxa", "unknown", "unidentified marine bacterioplankton")
    
    for(col in tax_cols) {
      for(pattern in patterns) {
        rows <- grep(pattern, counts[[col]])
        if(length(rows) > 0) {
          counts[rows, col] <- paste("Unknown", counts[rows, max(1, which(tax_cols == col)-1)])
        }
      }
    }
    
    datalist$Count_Data <- counts
    return(datalist)
  }
  
  # Pipeline principale (semplificata ma robusta)
  tryCatch({
    data_import <- data_select(file_ASV, kingdom)
    
    # Aggiungi conteggi totali
    data_import$Meta_Data$Counts_Total <- colSums(
      dplyr::select(data_import$Count_Data, where(is.numeric)), 
      na.rm = TRUE
    )[match(data_import$Meta_Data$Sample_ID, 
            names(dplyr::select(data_import$Count_Data, where(is.numeric))))]
    
    # Filtro minimo conteggi
    if(!is.null(min_counts)) {
      valid_samples <- data_import$Meta_Data %>% 
        dplyr::filter(Counts_Total > min_counts) %>% 
        dplyr::pull(Sample_ID)
      
      data_import$Count_Data <- data_import$Count_Data %>%
        dplyr::select(OTU_ID, all_of(valid_samples))
      
      data_import$Meta_Data <- data_import$Meta_Data %>%
        dplyr::filter(Sample_ID %in% valid_samples)
    }
    
    # Rarefazione
    if(!is.null(rare_lim)) {
      numeric_data <- data_import$Count_Data %>%
        dplyr::select(where(is.numeric)) %>%
        as.matrix()
      
      rarefied <- vegan::rrarefy(t(numeric_data), rare_lim) %>% t()
      
      data_import$Count_Data <- data_import$Count_Data %>%
        dplyr::select(!where(is.numeric)) %>%
        dplyr::bind_cols(as.data.frame(rarefied)) %>%
        dplyr::filter(rowSums(select(., where(is.numeric))) > 0)
    }
    
    # Filtro abbondanza
    if(abundance_filter) {
      data_import <- filter_abundance(data_import)
    }
    
    # Correzione tassonomia
    data_import <- correct_ambiguous(data_import)
    
    return(data_import)
    
  }, error = function(e) {
    stop("Errore in import_data: ", e$message)
  })
}



get_mechanism_prop2 <- function(datalist, bNTI, RC_BC) {
  bNTI_mod <- bNTI %>%
    magrittr::set_colnames(rownames(.)) %>%
    as.matrix() 
  
  bNTI_mod[lower.tri(bNTI_mod)] <- NA
  bNTI_mod[diag(bNTI_mod)] <- NA
  
  RC_BC_mod <- RC_BC %>%
    magrittr::set_rownames(rownames(bNTI_mod)) %>%
    magrittr::set_colnames(colnames(bNTI_mod)) %>%
    as.matrix()
  
  RC_BC_mod[lower.tri(RC_BC_mod)] <- NA
  RC_BC_mod[diag(RC_BC_mod)] <- NA
  
  merged <- bNTI_mod %>%
    reshape2::melt() %>%
    with(., cbind(., reshape2::melt(RC_BC_mod))) %>%
    magrittr::set_colnames(c("Sample_ID", "To_Sample", "bNTI", "1", "2", "RC_BC")) %>%
    dplyr::select(Sample_ID, To_Sample, bNTI, RC_BC) %>%
    filter(!is.na(bNTI)) %>%
    dplyr::mutate(Mechanism = ifelse(bNTI > 2, "Heterogeneous Selection", 
                                     ifelse(bNTI < -2, "Homogeneous Selection",
                                            ifelse(RC_BC < 0.95 & RC_BC > -0.95, "Drift",
                                                   ifelse(RC_BC > 0.95, "Dispersal Limitation", "Homogenising Dispersal"))))) %>%
    group_by(Mechanism) %>%
    dplyr::summarize(Num = n()) %>%
    ungroup() %>%
    dplyr::mutate(Num = Num/sum(Num)) %>%
    dplyr::mutate(Mechanism_Group = ordered(ifelse(Mechanism == "Homogeneous Selection" | Mechanism == "Heterogeneous Selection", "Selection",
                                                   ifelse(Mechanism == "Homogenising Dispersal" | Mechanism == "Dispersal Limitation", "Dispersal",
                                                          "Drift")),
                                            levels = c("Drift", "Dispersal", "Selection"))) %>%
    dplyr::mutate(Mechanism = ordered(Mechanism, levels = c("Homogeneous Selection", "Heterogeneous Selection", 
                                                            "Homogenising Dispersal", "Dispersal Limitation", "Drift")))
  
  return(merged)
}