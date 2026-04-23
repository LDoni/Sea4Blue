##quello fatto da me 
sparCC_wrapper_with_pb <- function(datalist, envir_filter = NULL, n_boot = 99, frac = FALSE) {
  
  library(progress)
  
  shuffle <- function(count) {
    tmp <- count
    for (i in 1:ncol(tmp)) {
      tmp[,i] <- tmp[sample(nrow(tmp), replace = T), i]
    }
    return(tmp)
  }
  
  devtools::source_url("https://raw.githubusercontent.com/huayingfang/CCLasso/master/R/SparCC.R")
  
  envir_filter <- enquo(envir_filter)
  
  OTU_ID <- datalist %>%
    with(., if (!rlang::quo_is_null(envir_filter)) filter_station_datalist(., !!envir_filter) else .) %>%
    .$Count_Data %>%
    .$OTU_ID
  
  count <- datalist %>%
    with(., if (!rlang::quo_is_null(envir_filter)) filter_station_datalist(., !!envir_filter) else .) %>%
    .$Count_Data %>%
    select_if(is.numeric) %>%
    t() %>%
    as_tibble()
  
  cat("=== ANALISI SparCC ===\n")
  cat("Inizio:", format(Sys.time()), "\n")
  cat("OTU:", length(OTU_ID), "| Campioni:", nrow(count), "| Bootstrap:", n_boot, "\n")
  cat("-----------------------------\n")
  
  # 1. Correlazioni reali
  cat("1. Calcolo correlazioni reali...\n")
  pb1 <- progress_bar$new(total = 1, format = "[:bar] :percent")
  pb1$tick()
  
  if (frac) {
    SparCC_true_cor <- count %>% SparCC.frac() %>% .$cor.w  
  } else {
    SparCC_true_cor <- count %>% SparCC.count() %>% .$cor.w  
  }
  
  # 2. Bootstrap
  cat("2. Bootstrap shuffle...\n")
  pb2 <- progress_bar$new(total = n_boot, format = "[:bar] :percent (:eta rimasti)")
  
  boot_list <- list()
  for (i in 1:n_boot) {
    boot_list[[i]] <- shuffle(count)
    pb2$tick()
  }
  
  # 3. Calcolo p-value
  cat("3. Calcolo p-value...\n")
  pb3 <- progress_bar$new(total = n_boot, format = "[:bar] :percent (:eta rimasti)")
  
  SparCC_boot_cor <- list()
  for (i in 1:n_boot) {
    if (frac) {
      cor_mat <- SparCC.frac(boot_list[[i]])$cor.w
    } else {
      cor_mat <- SparCC.count(boot_list[[i]])$cor.w
    }
    
    SparCC_boot_cor[[i]] <- abs(cor_mat) >= abs(SparCC_true_cor) & 
      sign(cor_mat) == sign(SparCC_true_cor)
    pb3$tick()
  }
  
  # Somma i risultati
  SparCC_boot_cor <- Reduce(`+`, SparCC_boot_cor) / n_boot
  
  colnames(SparCC_true_cor) <- OTU_ID
  rownames(SparCC_true_cor) <- OTU_ID
  
  cat("-----------------------------\n")
  cat("Fine:", format(Sys.time()), "\n")
  cat("Correlazioni significative (p<0.05):", 
      sum(SparCC_boot_cor < 0.05, na.rm = TRUE), "\n")
  
  return(list(cor = SparCC_true_cor,
              pVal = SparCC_boot_cor))
}






###quello vero
sparCC_wrapper <- function(datalist, envir_filter = NULL, n_boot = 99, frac = FALSE) {
  
  shuffle <- function(count) {
    
    tmp <- count
    
    for (i in 1:ncol(tmp)) {
      tmp[,i] <- tmp[sample(nrow(tmp), replace = T), i]
    }
    return(tmp)
  }
  
  devtools::source_url("https://raw.githubusercontent.com/huayingfang/CCLasso/master/R/SparCC.R")
  
  envir_filter <- enquo(envir_filter)
  
  OTU_ID <- datalist %>%
    with(., if (!rlang::quo_is_null(envir_filter)) filter_station_datalist(., !!envir_filter) else .) %>%
    .$Count_Data %>%
    .$OTU_ID
  
  count <- datalist %>%
    with(., if (!rlang::quo_is_null(envir_filter)) filter_station_datalist(., !!envir_filter) else .) %>%
    .$Count_Data %>%
    select_if(is.numeric) %>%
    t() %>%
    as_tibble()
  
  if (frac) {
    SparCC_true_cor <- count %>%
      SparCC.frac() %>%
      .$cor.w  
  } else {
    SparCC_true_cor <- count %>%
      SparCC.count() %>%
      .$cor.w  
  }
  
  boot_list <- rep(list(count), n_boot) %>%
    bplapply(., shuffle)
  
  SparCC_boot_cor <- boot_list %>%
    bplapply(., function(x) {
      if (frac) SparCC.frac(x) else SparCC.count(x) }) %>%
    map(., "cor.w") %>%
    map(., function(x) abs(x) >= abs(SparCC_true_cor) & sign(x) == sign(SparCC_true_cor)) %>%
    reduce(`+`)/n_boot
  
  colnames(SparCC_true_cor) <- OTU_ID
  rownames(SparCC_true_cor) <- OTU_ID
  
  return(list(cor = SparCC_true_cor,
              pVal = SparCC_boot_cor))
  
}