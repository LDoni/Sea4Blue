#### Import Function ####
# Start of every Analysis Pipeline
# 
# Define:
# file_ASV = Location of folder structure (top level of substructure -> see Github ReadMe)
# kingdom = Prokaryotes, Chloroplasts, Eukaryotes (Prok, Chloroplast, Euk)
# rare_lim = Integer defining the rarefying level. If NULL no rarefying will be done
# drop_rare = Logical defining if samples below rare_lim will be dropped from analysis
# abundance_filter = Logical defining if abundance filter after Milici et al. 2016 should be applied
# min_count = Integer defining the minimum count number a sample should have. NULL for no filtering

import_data <- function(file_ASV, kingdom = "Prok", rare_lim = NULL, drop_rare = T, 
                        abundance_filter = F, min_counts = NULL) {
  
  data_select <- function(file_ASV, kingdom = "Prok") {
    
    Count_Data <- suppressMessages(read_delim(paste(file_ASV, "Processed/", kingdom, "/Full_", kingdom,"_Count.tsv", sep = ""), 
                                              del = "\t")) %>%
      select(-grep("[Mm]ock", names(.))) %>%
      select(-grep("NC", names(.))) 
    
    Meta_Data <- suppressMessages(read_delim(paste(file_ASV, "Meta_Data/", kingdom, "/Meta_Data.tsv", sep = ""), 
                                             del = "\t"))
    
    return(list(Meta_Data = Meta_Data, 
                Count_Data = Count_Data))
    
  }
  
  filter_abundance <- function(datalist) {
    
    counts <- datalist$Count_Data
    prop <- datalist$Count_Data %>%
      mutate_if(is.numeric,
                function(x) x/sum(x))
    
    counts_filtered <- counts %>%
      filter(((rowSums(select_if(counts, is.numeric))/sum(rowSums(select_if(counts, is.numeric)))) > 0.00001) &
               (apply(select_if(prop, is.numeric),1,max) > 0.01) |
               ((apply(select_if(prop, is.numeric) > 0.001, 1, sum) / ncol(select_if(prop, is.numeric))) > 0.02) |
               ((apply(select_if(prop, is.numeric) > 0, 1, sum) / ncol(select_if(prop, is.numeric))) > 0.05))
    
    datalist$Count_Data <- counts_filtered
    
    return(datalist)
    
  }
  
  rarefy_datalist <- function(datalist, rare_lim, drop = F) {
    
    count_rared <- datalist$Count_Data %>%
      select_if(is.numeric) %>%
      select_if(colSums(.) >= ifelse(drop, rare_lim, 0)) %>%
      t() %>% 
      vegan::rrarefy(., rare_lim) %>%
      t() %>%
      as_tibble() %>%
      bind_cols(select_if(datalist$Count_Data, is.character), .) %>%
      filter(rowSums(select_if(., is.numeric)) > 0)
    
    meta_subset <- datalist$Meta_Data %>%
      dplyr::slice(match(names(select_if(count_rared, is.numeric)), datalist$Meta_Data$Sample_ID))
    
    datalist$Count_Data <- count_rared
    datalist$Meta_Data <- meta_subset
    
    return(datalist)
    
  }
  
  correct_ambiguous <- function(datalist, fromTaxLvl = 8) {
    
    replacer <- function(Count_Data, taxLvl, replaceLvl, pattern) {
      Count_Data[grep(pattern, x = as_vector(Count_Data[, taxLvl])), taxLvl] <- paste0("Unknown ", as_vector(Count_Data[grep(pattern, x = as_vector(Count_Data[, taxLvl])), replaceLvl]))
      return(Count_Data)
    }
    
    for (taxLvl in fromTaxLvl:1) {
      
      for (i in (taxLvl-1):1) {
        
        datalist$Count_Data <- replacer(datalist$Count_Data, taxLvl, replaceLvl = i, pattern = 'bacterium$') %>%
          replacer(., taxLvl, replaceLvl = i, pattern = "uncultured") %>%
          replacer(., taxLvl, replaceLvl = i, pattern = "metagenome") %>%
          replacer(., taxLvl, replaceLvl = i, pattern = "unidentified") %>%
          replacer(., taxLvl, replaceLvl = i, pattern = "Ambiguous_taxa") %>%
          replacer(., taxLvl, replaceLvl = i, pattern = "unknown") %>%
          replacer(., taxLvl, replaceLvl = i, pattern = "unidentified marine bacterioplankton") %>%
          replacer(., taxLvl, replaceLvl = i, pattern = "^Unknown.*bacterium$") %>%
          replacer(., taxLvl, replaceLvl = i, pattern = "^Uncultured.*bacterium$") 
        
      }
      
    }  
    
    return(datalist)
    
  }
  
  data_import <- data_select(file_ASV, kingdom = kingdom) %>%
    mutate_meta_datalist(Counts_Total = colSums(select_if(.$Count_Data, is.numeric))) %>%
    with(., if (!is.null(min_counts)) filter_station_datalist(., Counts_Total > !!min_counts) else .) %>%
    with(., if (!is.null(rare_lim)) rarefy_datalist(., rare_lim, drop_rare) else .) %>%
    with(., if (abundance_filter) filter_abundance(.) else .) %>%
    correct_ambiguous()
  
  return(data_import)
  
}




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

