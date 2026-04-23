plot_mechanism_vs_env <- function(
    datalist,
    bNTI,
    RC_BC,
    env_var,
    breaks,
    env_label,
    palette,
    transform_fun = NULL
) {
  
  ## ---------- STEP 1: calcolo edge + Δenv ----------
  detailled <- get_mechanism_prop_with_env(
    datalist = datalist,
    bNTI = bNTI,
    RC_BC = RC_BC,
    env_var = env_var,
    breaks = breaks,
    transform_fun = transform_fun
  )
  
  diff_grp_col <- paste0(env_var, "_Diff_Grp")
  
  ## ---------- STEP 2: aggregazione ----------
  mech_env <- detailled %>%
    dplyr::group_by(.data[[diff_grp_col]], Mechanism) %>%
    dplyr::summarise(N = n(), .groups = "drop") %>%
    dplyr::group_by(.data[[diff_grp_col]]) %>%
    dplyr::mutate(Prop = N / sum(N)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Mechanism = factor(
        Mechanism,
        levels = c("Homogeneous Selection", "Heterogeneous Selection",
                   "Homogenising Dispersal", "Dispersal Limitation", "Drift"),
        labels = c("Homogeneous selection", "Heterogeneous selection",
                   "Homogenising dispersal", "Dispersal limitation", "Drift")
      )
    )
  
  ## ---------- STEP 3: plot ----------
  p <- ggplot(
    mech_env,
    aes(x = .data[[diff_grp_col]], y = Prop * 100, fill = Mechanism)
  ) +
    geom_bar(stat = "identity", colour = "black", size = 0.2) +
    scale_fill_manual(values = palette) +
    labs(
      x = env_label,
      y = "Proportion of mechanisms (%)",
      fill = "Assembly mechanism"
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(colour = "black"),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(colour = "black"),
      legend.position = "right"
    )
  
  ## ---------- OUTPUT ----------
  return(list(
    data = mech_env,
    plot = p
  ))
}

datalist_Atlantic$Meta_Data <- datalist_Atlantic$Meta_Data %>%
  mutate(
    uo = as.numeric(uo),
    vo = as.numeric(vo),
    current_speed = sqrt(uo^2 + vo^2)
  )

temp_out <- plot_mechanism_vs_env(
  datalist = datalist_Atlantic,
  bNTI = as.matrix(read.csv("output/Community_Mechanisms/Prokaryotes_Atlantic_weighted_bNTI.csv", row.names = 1)),
  RC_BC = as.matrix(read.csv("output/Community_Mechanisms/Raup_Crick_Prok.csv", row.names = 1)),
  env_var = "thetao",
  breaks = seq(0, 28, 1),
  env_label = "Temperature difference (°C)",
  palette = cbbPalette[c(7,5,3,6,4)]
)

temp_out$plot
sal_out <- plot_mechanism_vs_env(
  datalist_Atlantic,
  bNTI,
  RC_BC,
  env_var = "so",
  breaks = seq(0, 5, 0.2),
  env_label = "Salinity difference (PSU)",
  palette = cbbPalette[c(7,5,3,6,4)]
)

sal_out$plot


chl_out <- plot_mechanism_vs_env(
  datalist_Atlantic,
  bNTI,
  RC_BC,
  env_var = "chl",
  breaks = seq(0, 3, 0.2),
  env_label = expression(log[10]~"Chlorophyll-a difference"),
  palette = cbbPalette[c(7,5,3,6,4)],
  transform_fun = function(x) log10(x + 1e-6)
)

chl_out$plot
current_out <- plot_mechanism_vs_env(
  datalist_Atlantic,
  bNTI,
  RC_BC,
  env_var = "current_speed",
  breaks = seq(0, 1.5, 0.05),
  env_label = "Current speed difference (m s⁻¹)",
  palette = cbbPalette[c(7,5,3,6,4)]
)

current_out$plot



ggplot(temp_out$data,
       aes(x = thetao_Diff_Grp, y = Prop,
           colour = Mechanism)) +
  geom_point(size = 3) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 4),
    se = TRUE
  ) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = expression(Delta*"Temperature (°C)"),
    y = "Proportion of assembly mechanisms"
  ) +
  theme_classic()


