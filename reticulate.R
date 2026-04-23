


py_run_string("
import copernicusmarine
copernicusmarine.login(
    username = 'ldoni',  # Sostituisci con il tuo username
    password = 'Carambola999'  # Sostituisci con la tua password reale
)
")





library(tidyverse)
 
library(reticulate)
reticulate::virtualenv_list()



#virtualenv_remove("something", confirm = FALSE)  # confirm = FALSE evita la richiesta di conferma
install_python() 



virtualenv_create(envname = "something")
virtualenv_install("something", packages = c("copernicusmarine"))
use_virtualenv("something", required = TRUE)
 

py_install("copernicusmarine")

cm <- import("copernicusmarine")
 
 cm$CopernicusMarineDataset

 

result <- cm$subset(
  dataset_id = "cmems_mod_ibi_phy_anfc_0.027deg-2D_PT1H-m",
  start_datetime = "2022-05-23",
  end_datetime = "2022-05-23",
  variables = list("thetao"),  # Temperatura superficiale
  minimum_longitude = -9.64671906997131,
  maximum_longitude = 7.364793639893833,
  minimum_latitude = 36.7143279919028,
  maximum_latitude = 41.884404655537004,
  output_filename = "dati_temperatura.nc"
)
















library(reticulate)
library(dplyr)
library(readr)
library(tidync)
library(purrr)
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

cm <- import("copernicusmarine")

# 1. Caricamento dati campioni --------------------------------------------
  # oppure crea manualmente il tuo dataframe

df <- df %>% 
  mutate(data = as.Date(data, format = "%Y-%m-%d"))  # Assicurati del formato corretto

# 2. Definizione datasets Copernicus --------------------------------------
datasets <- list(
  FISICO = list(
    id = "cmems_mod_glo_phy_anfc_0.083deg_P1D-m",
    vars = c("thetao", "so", "uo", "vo")
  ),
  BIO = list(
    id = "cmems_mod_glo_bgc_anfc_0.25deg_P1D-m",
    vars = c("chl", "no3", "po4", "si")
  )
)

# 3. Funzione per scaricare un singolo punto ------------------------------
download_single_point <- function(sample_id, date, lon, lat, dataset_id, variables) {
  message(sprintf("Processing %s (%s) | Dataset: %s", sample_id, date, dataset_id))
  
  output_file <- tempfile(fileext = ".nc")
  
  tryCatch({
    result <- cm$subset(
      dataset_id = dataset_id,
      variables = variables,
      start_datetime = as.character(date),
      end_datetime = as.character(date),
      minimum_longitude = lon - 0.05,
      maximum_longitude = lon + 0.05,
      minimum_latitude = lat - 0.05,
      maximum_latitude = lat + 0.05,
      minimum_depth = 0,
      maximum_depth = 5,
      output_filename = output_file
    )
    
    if (file.exists(output_file)) {
      dati <- tidync::tidync(output_file) %>%
        hyper_tibble() %>%
        mutate(
          Sample = sample_id,
          Date = as.Date(date),
          Dataset = dataset_id
        )
      file.remove(output_file)
      return(dati)
    } else {
      warning(paste("Fallito per", sample_id, "| file non trovato"))
      return(NULL)
    }
  }, error = function(e) {
    warning(paste("Errore critico per", sample_id, ":", e$message))
    return(NULL)
  })
}


# 4. Download multiplo -----------------------------------------------------
risultati <- list()

for (ds_name in names(datasets)) {
  ds <- datasets[[ds_name]]
  message("Inizio elaborazione dataset: ", ds_name)
  
  dati_ds <- list()
  
  for (i in seq_len(nrow(df))) {
    riga <- df[i, ]
    dati <- download_single_point(
      sample_id = riga$Nominativo.campione.,
      date = riga$data,
      lon = riga$Longitude,
      lat = riga$Latitude,
      dataset_id = ds$id,
      variables = ds$vars
    )
    
    if (!is.null(dati)) {
      dati_ds[[i]] <- dati
    }
    
    if (i %% 5 == 0) Sys.sleep(1)
  }
  
  risultati[[ds_name]] <- bind_rows(dati_ds)
  message(sprintf(
    "Completato %s: %d/%d campioni riusciti",
    ds_name,
    sum(!sapply(dati_ds, is.null)),
    nrow(df)
  ))
}

# 5. Unione e salvataggio finale -------------------------------------------
dati_completi <- risultati %>% 
  reduce(full_join, by = c("longitude", "latitude", "depth", "time", "Sample", "Date")) %>%
  mutate(
    current_speed = sqrt(uo^2 + vo^2),
    current_dir_deg = (atan2(vo, uo) * 180/pi) %% 360
  )

# 6. Salva in CSV e RDS ----------------------------------------------------
write_csv(dati_completi, "dati_oceanici_completi.csv")
saveRDS(dati_completi, "dati_oceanici.rds")

# 7. Report finale ---------------------------------------------------------
success_rates <- sapply(risultati, function(x) round(100 * nrow(x)/nrow(df), 1))
message("\nRiepilogo finale:")
message(paste("-", names(success_rates), ":", success_rates, "% successo", collapse = "\n"))













