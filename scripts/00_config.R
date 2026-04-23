sea4blue_config <- function(
  source_root = "/media/shared1/atlantic_sea4blue/16S_Illumina/Lapo_18May",
  project_root = "/home/userbio/Sea4Blue"
) {
  list(
    source_root = normalizePath(source_root, mustWork = FALSE),
    project_root = normalizePath(project_root, mustWork = FALSE),
    figures_dir = file.path(project_root, "results", "figures"),
    tables_dir = file.path(project_root, "results", "tables"),
    ocean_dir = file.path(project_root, "results", "oceanography"),
    paper_dir = file.path(source_root, "paper"),
    quarto_dir = file.path(project_root, "supplement"),
    paths = list(
      phyloseq_rds = file.path(source_root, "ps.rds"),
      env_csv = file.path(source_root, "dati_scaricati_copernicus.csv"),
      metadata_csv = file.path(source_root, "metadata_sea4blue.csv"),
      sst_nc = file.path(source_root, "AQUA_MODIS.20220101_20221231.L3m.YR.SST.sst.4km.nc"),
      longhurst_shp = file.path(source_root, "longhurst_v4_2010", "Longhurst_world_v4_2010.shp"),
      btnti_csv = file.path(source_root, "output", "Community_Mechanisms", "Prokaryotes_Atlantic_weighted_bNTI.csv"),
      rc_bray_csv = file.path(source_root, "output", "Community_Mechanisms", "Raup_Crick_Prok.csv"),
      network_properties_csv = file.path(source_root, "output", "Network_Properties.csv"),
      sparcc_cor_csv = file.path(source_root, "output", "SparCC", "Cor_SparCC_Prok_all.csv"),
      sparcc_pval_csv = file.path(source_root, "output", "SparCC", "Pval_SparCC_Prok_all.csv"),
      sparcc_cluster_rds = file.path(source_root, "output", "SparCC", "cluster.rds"),
      graphml = file.path(source_root, "output", "SparCC_Network_pos.graphml"),
      stations_csv = file.path(source_root, "parcels", "stations.csv"),
      currents_nc = file.path(source_root, "parcels", "data_large", "cmems_mod_glo_phy_my_0.083deg_P1D-m_uo-vo_90.00W-0.00E_25.00N-45.00N_0.49m_2022-01-01-2022-12-29.nc")
    )
  )
}

ensure_result_dirs <- function(cfg) {
  dir.create(cfg$figures_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cfg$tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cfg$ocean_dir, recursive = TRUE, showWarnings = FALSE)
  invisible(cfg)
}
