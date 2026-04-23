# Installa e carica i pacchetti necessari

library(CopernicusMarine)
library(tidync)
library(dplyr)
library(lubridate)

 

df <- ps %>%
  subset_samples(!(sample.description %in% c("Laposnegative", "Lapospositive"))) %>%
  meta() %>%
  dplyr::select(Nominativo.campione., Latitude, Longitude,data)
df$data <- as.Date(df$data, format = "%d/%m/%Y")  

# Scarica i dati per ogni punto del tuo dataframe




 
library(reticulate)



virtualenv_create(envname = "something")
virtualenv_install("something", packages = c("copernicusmarine"))
use_virtualenv("something", required = TRUE)
py_install("copernicusmarine")

py_run_string("
import copernicusmarine
copernicusmarine.login(
    username = 'ldoni',  
    password = 'Carambola999'  
)
")



cm <- import("copernicusmarine")



result <- cm$subset(
  dataset_id = "cmems_mod_ibi_phy_anfc_0.027deg-2D_PT1H-m",
  start_datetime = "2022-11-23",  # Data di inizio valida
  end_datetime = "2022-11-23",
  variables = list("thetao"),
  minimum_longitude = -9.6,       # Entro i limiti [-19.08, 5.08]
  maximum_longitude = 5.0,
  minimum_latitude = 36.71,
  maximum_latitude = 41.88,
  output_filename = "dati_corretti.nc"
)




library(reticulate)
library(dplyr)
library(purrr)
library(readr)

# 1. Configurazione Python -------------------------------------------------
# Assicurati di aver fatto login prima via terminale:
# python -c "import copernicusmarine; copernicusmarine.login()"

cm <- import("copernicusmarine")

# 2. Caricamento dati -----------------------------------------------------
df <- df %>% 
  mutate(data = as.Date(data, format = "%d/%m/%Y"))

# 3. Definizione dataset --------------------------------------------------
datasets <- list(
  FISICO = list(
    id = "cmems_mod_glo_phy_anfc_0.083deg_P1D-m",
    vars = c("thetao", "so", "uo", "vo")
  ),
  BIO = list(
    id = "cmems_mod_glo_bgc_anfc_0.25deg_P1D-m",
    vars = c("chl")
  )
)

# 4. Funzione di download con logging --------------------------------------
download_single_point <- function(sample_id, date, lon, lat, dataset_id, variables) {
  message(sprintf("Processing %s (%s) | Dataset: %s", sample_id, date, dataset_id))
  
  output_file <- tempfile(fileext = ".nc")
  
  tryCatch({
    # Chiamata Python diretta con gestione errori
    py_run_string(sprintf(
      "
      try:
          result = copernicusmarine.subset(
              dataset_id='%s',
              variables=%s,
              start_datetime='%s',
              end_datetime='%s',
              minimum_longitude=%f,
              maximum_longitude=%f,
              minimum_latitude=%f,
              maximum_latitude=%f,
              minimum_depth=0,
              maximum_depth=5,
              output_filename='%s',
              force_download=True
          )
      except Exception as e:
          result = str(e)
      ",
      dataset_id,
      reticulate::r_to_py(variables),
      date, date,
      lon - 0.05, lon + 0.05,
      lat - 0.05, lat + 0.05,
      output_file
    ))
    
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
      warning(paste("Fallito per", sample_id, "| Errore:", py$result))
      return(NULL)
    }
  }, error = function(e) {
    warning(paste("Errore critico per", sample_id, ":", e$message))
    return(NULL)
  })
}

# 5. Loop principale ------------------------------------------------------
risultati <- list()

for (ds_name in names(datasets)) {
  ds <- datasets[[ds_name]]
  message("Inizio elaborazione dataset: ", ds_name)
  
  dati_ds <- list()
  
  for (i in seq_len(nrow(df))) {
    riga <- df[i, ]
    dati <- download_single_point(
      riga$Nominativo.campione.,
      riga$data,
      riga$Longitude,
      riga$Latitude,
      ds$id,
      ds$vars
    )
    
    if (!is.null(dati)) {
      dati_ds[[i]] <- dati
    }
    
    # Pausa per evitare rate limiting
    if (i %% 5 == 0) Sys.sleep(1)
  }
  
  risultati[[ds_name]] <- bind_rows(dati_ds)
  message(sprintf(
    "Completato %s: %d/%d campioni successo",
    ds_name,
    sum(!sapply(dati_ds, is.null)),
    nrow(df)
  )
}

# 6. Elaborazione risultati -----------------------------------------------
dati_completi <- risultati %>% 
  reduce(full_join, by = c("longitude", "latitude", "depth", "time", "Sample", "Date")) %>%
  mutate(
    current_speed = sqrt(uo^2 + vo^2),
    current_dir_deg = (atan2(vo, uo) * 180/pi) %% 360
  )

# 7. Salvataggio ---------------------------------------------------------
write_csv(dati_completi, "dati_oceanici_completi.csv")
saveRDS(dati_completi, "dati_oceanici.rds")

# Report finale
success_rates <- sapply(risultati, function(x) round(100 * nrow(x)/nrow(df), 1)
                        message("\nRiepilogo finale:")
                        message(paste("-", names(success_rates), ":", success_rates, "% successo", collapse = "\n"))



                        library(reticulate)
                        library(dplyr)
                        library(readr)
                        library(tidync)
                        library(purrr)
                        
                        copernicusmarine <- import("copernicusmarine")
                        
                        
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
                            py_run_string(sprintf(
                              "
      try:
          result = copernicusmarine.subset(
              dataset_id='%s',
              variables=%s,
              start_datetime='%s',
              end_datetime='%s',
              minimum_longitude=%f,
              maximum_longitude=%f,
              minimum_latitude=%f,
              maximum_latitude=%f,
              minimum_depth=0,
              maximum_depth=5,
              output_filename='%s',
              force_download=True
          )
      except Exception as e:
          result = str(e)
      ",
                              dataset_id,
                              reticulate::r_to_py(variables),
                              date, date,
                              lon - 0.05, lon + 0.05,
                              lat - 0.05, lat + 0.05,
                              output_file
                            ))
                            
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
                              warning(paste("Fallito per", sample_id, "| Errore:", py$result))
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
                        

                        py_run_string("
from datetime import datetime
info = copernicusmarine.get_dataset_info('cmems_mod_glo_phy_anfc_0.083deg_P1D-m')
print(f'Temporange: {info.time_coverage_start} to {info.time_coverage_end}')
")

test_point <- df[1, ]

test_download <- function() {
  output <- tempfile(fileext = ".nc")
  
  py_run_string(sprintf("
try:
    result = copernicusmarine.subset(
        dataset_id='cmems_mod_glo_phy_anfc_0.083deg_P1D-m',
        variables=['thetao'],
        start_datetime='%s',
        end_datetime='%s',
        minimum_longitude=%f,
        maximum_longitude=%f,
        minimum_latitude=%f,
        maximum_latitude=%f,
        minimum_depth=0,
        maximum_depth=5,
        output_filename='%s',
        force_download=True
    )
    print('Success! File:', result)
except Exception as e:
    print('Error:', str(e))
", 
  test_point$data, test_point$data,
  test_point$Longitude - 0.05, test_point$Longitude + 0.05,
  test_point$Latitude - 0.05, test_point$Latitude + 0.05,
  output))
  
  if (file.exists(output)) {
    print("File scaricato con successo!")
    print(list.files(tempdir()))
    file.remove(output)
  }
}

test_download()               
