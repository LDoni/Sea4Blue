library(reticulate)
library(ncdf4)
library(dplyr)
library(purrr)

# Inizializza Copernicus Marine
cm <- import("copernicusmarine")

# I tuoi campioni
df <- data.frame(
  Nominativo.campione. = c(
    "eDNA1", "eDNA18", "eDNA19", "eDNA25", "eDNA26", "eDNA27", "eDNA28",
    "eDNA3", "eDNA4", "eDNA6", "eDNA7", "eDNA9", "eDNA15", "eDNA16", "eDNA17"
  ),
  Latitude = c(
    29.30070, 35.79043, 36.99627, 35.89713, 36.13633, 36.25035, 36.07880,
    30.65897, 31.18923, 31.75563, 31.82482, 31.99925, 34.21340, 34.55028, 35.44867
  ),
  Longitude = c(
    -77.936917, -34.595667, -32.055683, -17.887783, -14.540200, -11.687733, -7.874533,
    -75.508383, -74.859950, -69.329683, -67.287833, -62.128333, -44.076200, -40.685317, -37.395217
  ),
  data = as.Date(c(
    "2022-05-23", "2022-06-09", "2022-06-10", "2022-06-23", "2022-06-24",
    "2022-06-25", "2022-06-26", "2022-05-25", "2022-06-26", "2022-05-28",
    "2022-05-29", "2022-05-31", "2022-06-06", "2022-06-07", "2022-06-08"
  ))
)

# Dataset CORRETTI da usare
datasets <- list(
  physics = list(
    id = "cmems_mod_glo_phy_myint_0.083deg_P1D-m",
    vars = c("thetao", "so", "uo", "vo", "zos")
  ),
  bgc_model = list(
    id = "cmems_mod_glo_bgc_my_0.25deg_P1D-m",
    vars = c("chl", "o2", "no3", "po4", "si")
  ))

# Funzione migliorata per estrarre valori
extract_nc_values <- function(ncfile, vars) {
  nc <- nc_open(ncfile)
  out <- list()
  for (v in vars) {
    if (v %in% names(nc$var)) {
      val <- ncvar_get(nc, v)
      out[[v]] <- as.numeric(val) # perché ora c’è un solo valore
    } else {
      out[[v]] <- NA
    }
  }
  nc_close(nc)
  return(out)
}

# Funzione principale corretta
get_sample_data <- function(sample, dataset_info) {
  result <- as.list(sample)
  
  for (ds_name in names(dataset_info)) {
    ds <- dataset_info[[ds_name]]
    message("Scarico ", ds_name, " per ", sample[["Nominativo.campione."]])
    
    tmpfile <- tempfile(fileext = ".nc")
    
    tryCatch({
      # Chiamata corretta a subset
      cm$subset(
        dataset_id = ds$id,
        variables = ds$vars,
        minimum_longitude = as.numeric(sample["Longitude"]),
        maximum_longitude = as.numeric(sample["Longitude"]),
        minimum_latitude = as.numeric(sample["Latitude"]),
        maximum_latitude = as.numeric(sample["Latitude"]),
        start_datetime = format(as.Date(sample[["data"]]), "%Y-%m-%d"),
        end_datetime = format(as.Date(sample[["data"]]), "%Y-%m-%d"),
        minimum_depth = 0.493,
        maximum_depth = 0.5058,
        coordinates_selection_method = "nearest",
        output_filename = tmpfile
      )
      
      if (file.exists(tmpfile) && file.size(tmpfile) > 0) {
        vals <- extract_nc_values(tmpfile, ds$vars)
        result <- c(result, vals)
      } else {
        warning("File non scaricato per ", sample[["Nominativo.campione."]], " - ", ds_name)
        # Aggiungi NA per tutte le variabili
        for (v in ds$vars) {
          result[[v]] <- NA
        }
      }
      
    }, error = function(e) {
      message("Errore per ", sample[["Nominativo.campione."]], " - ", ds_name, ": ", e$message)
      for (v in ds$vars) {
        result[[v]] <- NA
      }
    })
    
    # Pulisci file temporaneo
    if (file.exists(tmpfile)) file.remove(tmpfile)
  }
  
  return(as.data.frame(result))
}

# Esegui l'estrazione con gestione errori
all_data <- map_dfr(1:nrow(df), function(i) {
  sample <- df[i, ]
  tryCatch({
    get_sample_data(sample, datasets)
  }, error = function(e) {
    message("Errore completo per campione ", sample[["Nominativo.campione."]], ": ", e$message)
    # Restituisce riga con tutti NA per le variabili
    result <- as.list(sample)
    for (ds in datasets) {
      for (v in ds$vars) {
        result[[v]] <- NA
      }
    }
    return(as.data.frame(result))
  })
})

# Visualizza risultato
print(head(all_data))



#all_data_forsesbagliato<-all_data
 




cm$subset()










