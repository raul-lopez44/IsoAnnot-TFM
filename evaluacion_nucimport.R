# Evaluacion del modulo NLS (NucImport): Tabla 2 (grid search), Fig. 13-18 y Tabla 3.
# Requiere: ggplot2, dplyr, tidyr.
# Entrada: CSV con columnas Protein_ID, Set, Model, ImpProb, cNLSProb, Match.

library(ggplot2)
library(dplyr)
library(tidyr)

# --- Configuracion (ajustar rutas antes de ejecutar) ----------------------------
INPUT_CSV <- "tfm/_test_figures/df.csv"
OUT_DIR   <- "tfm/_test_figures/figuras_nls"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

UMBRAL_IMP  <- 0.60
UMBRAL_CNLS <- 0.35

TFM_GEOM <- 4.5
FIG_DPI  <- 300

# --- Carga de datos -------------------------------------------------------------
df <- read.csv(INPUT_CSV, stringsAsFactors = FALSE)
df$ImpProb  <- as.numeric(df$ImpProb)
df$cNLSProb <- as.numeric(df$cNLSProb)

clasificar <- function(datos, u_imp, u_cnls) {
  datos %>%
    mutate(
      Pasa_Imp  = ImpProb >= u_imp,
      Pasa_cNLS = cNLSProb >= u_cnls,
      Status = case_when(
        Set == "positivo" & Pasa_Imp & Pasa_cNLS            ~ "TP",
        Set == "positivo" & Pasa_Imp & !Pasa_cNLS           ~ "FN*",
        Set == "positivo" & !Pasa_Imp                       ~ "FN_L",
        Set == "negativo" & (!Pasa_Imp | !Pasa_cNLS)       ~ "TN",
        Set == "negativo" & Pasa_Imp & Pasa_cNLS            ~ "FP",
        TRUE ~ "Error"
      )
    )
}

# --- Tabla 2: grid search con restriccion posicional (Match == YES) ---------------
grid <- expand.grid(
  Umbral_Imp  = seq(0.1, 0.9, by = 0.05),
  Umbral_cNLS = seq(0.1, 0.9, by = 0.05)
)

evaluar_estricto <- function(u_imp, u_cnls, datos) {
  pred <- datos$ImpProb >= u_imp & datos$cNLSProb >= u_cnls
  tp   <- sum(datos$Set == "positivo" & pred & datos$Match == "YES")
  fn   <- sum(datos$Set == "positivo") - tp
  tn   <- sum(datos$Set == "negativo" & !pred)
  fp   <- sum(datos$Set == "negativo" & pred)

  sens <- if (tp + fn > 0) tp / (tp + fn) * 100 else 0
  spec <- if (tn + fp > 0) tn / (tn + fp) * 100 else 0
  prec <- if (tp + fp > 0) tp / (tp + fp) * 100 else 0
  f1   <- if (prec + sens > 0) 2 * prec * sens / (prec + sens) else 0
  youden <- sens + spec - 100

  c(Sens = sens, Spec = spec, F1 = f1, Youden = youden)
}

resultados_grid <- data.frame()
for (modelo in unique(df$Model)) {
  sub <- df %>% filter(Model == modelo)
  for (i in seq_len(nrow(grid))) {
    m <- evaluar_estricto(grid$Umbral_Imp[i], grid$Umbral_cNLS[i], sub)
    resultados_grid <- rbind(resultados_grid, data.frame(
      Modelo = modelo,
      Umbral_Imp = grid$Umbral_Imp[i],
      Umbral_cNLS = grid$Umbral_cNLS[i],
      Sensibilidad_Estricta = round(m["Sens"], 2),
      Especificidad = round(m["Spec"], 2),
      F1_Score = round(m["F1"], 2),
      Indice_Youden = round(m["Youden"], 2)
    ))
  }
}

write.csv(resultados_grid,
          file.path(OUT_DIR, "tabla_2_grid_search.csv"),
          row.names = FALSE)

resumen_optimos <- resultados_grid %>%
  group_by(Modelo) %>%
  slice_max(Indice_Youden, n = 1) %>%
  arrange(desc(Indice_Youden))
cat("\n--- Tabla 2: punto optimo por modelo (maximo Youden) ---\n")
print(resumen_optimos, row.names = FALSE)

# --- Fig. 13-15: barras con Sensibilidad, Especificidad y Sensibilidad Estricta ---
graficar_modelos <- function(datos, u_imp, u_cnls, subtitulo, archivo) {
  d <- clasificar(datos, u_imp, u_cnls)
  metricas <- d %>%
    group_by(Model) %>%
    summarise(
      Sensibilidad   = sum(Status == "TP") / sum(Set == "positivo") * 100,
      Especificidad  = sum(Status == "TN") / sum(Set == "negativo") * 100,
      Match_Correcto = sum(Status == "TP" & Match == "YES") / sum(Set == "positivo") * 100,
      .groups = "drop"
    ) %>%
    pivot_longer(c(Sensibilidad, Especificidad), names_to = "Metrica", values_to = "Valor")

  p <- ggplot(metricas, aes(x = Model, y = Valor, fill = Metrica)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    geom_text(aes(label = paste0(round(Valor, 1), "%")),
              position = position_dodge(0.9), vjust = -0.5, size = TFM_GEOM) +
    geom_bar(
      data = filter(metricas, Metrica == "Sensibilidad"),
      aes(x = Model, y = Match_Correcto, fill = "Sensibilidad Estricta"),
      stat = "identity", width = 0.4, position = position_nudge(x = 0.225)
    ) +
    geom_text(
      data = filter(metricas, Metrica == "Sensibilidad"),
      aes(x = Model, y = Match_Correcto, label = paste0(round(Match_Correcto, 1), "%")),
      position = position_nudge(x = 0.225),
      vjust = 1.5, size = TFM_GEOM, color = "white", fontface = "italic"
    ) +
    scale_fill_manual(
      values = c(
        "Especificidad" = "#F8766D",
        "Sensibilidad" = "#00BFC4",
        "Sensibilidad Estricta" = "#007376"
      ),
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    labs(
      title = "Comparativa de Modelos: Sensibilidad y Especificidad",
      subtitle = subtitulo,
      x = "Modelo", y = "Porcentaje (%)", fill = "Metrica"
    ) +
    theme_minimal(base_size = 14)

  ggsave(archivo, p, width = 12, height = 7, dpi = FIG_DPI, bg = "white")
  message("Figura guardada: ", archivo)
}

graficar_modelos(df, 0.60, 0.35,
  "Import Prob >= 0.60; cNLS Prob >= 0.35",
  file.path(OUT_DIR, "figura_13_modelos_0.60_0.35.png"))
graficar_modelos(df, 0.60, 0.60,
  "Import Prob >= 0.60; cNLS Prob >= 0.60",
  file.path(OUT_DIR, "figura_14_modelos_0.60_0.60.png"))
graficar_modelos(df, 0.90, 0.35,
  "Import Prob >= 0.90; cNLS Prob >= 0.35",
  file.path(OUT_DIR, "figura_15_modelos_0.90_0.35.png"))

# --- Fig. 16: capacidad discriminativa (Modelo 1) --------------------------------
df_m1 <- df %>% filter(Model == "M1")

p_disc <- ggplot(df_m1, aes(x = ImpProb, y = cNLSProb, color = Set)) +
  geom_point(alpha = 0.35, size = 1.8) +
  geom_vline(xintercept = UMBRAL_IMP, linetype = "dashed", color = "darkred") +
  geom_hline(yintercept = UMBRAL_CNLS, linetype = "dashed", color = "darkred") +
  scale_color_manual(values = c("positivo" = "#31a354", "negativo" = "#756bb1")) +
  labs(title = "Capacidad Discriminativa (Modelo 1)",
       x = "Import Probability", y = "cNLS Probability", color = "Conjunto") +
  theme_bw(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6))
ggsave(file.path(OUT_DIR, "figura_16_discriminativa_M1.png"),
       p_disc, width = 10, height = 7, dpi = FIG_DPI, bg = "white")

# --- Fig. 17: set positivo y acierto del motivo ------------------------------------
p_pos <- df_m1 %>%
  filter(Set == "positivo") %>%
  ggplot(aes(x = ImpProb, y = cNLSProb, color = Match)) +
  geom_point(alpha = 0.5, size = 2.2) +
  geom_vline(xintercept = UMBRAL_IMP, linetype = "dashed", color = "darkred") +
  geom_hline(yintercept = UMBRAL_CNLS, linetype = "dashed", color = "darkred") +
  scale_color_manual(values = c("YES" = "#2c7bb6", "NO" = "#d7191c")) +
  labs(title = "Distribucion del Set Positivo (Modelo 1)",
       x = "Import Probability", y = "cNLS Probability", color = "Match UniProt") +
  theme_bw(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6))
ggsave(file.path(OUT_DIR, "figura_17_set_positivo_M1.png"),
       p_pos, width = 10, height = 7, dpi = FIG_DPI, bg = "white")

# --- Fig. 18: diagnostico TP / FN* / FN_L ----------------------------------------
df_diag <- clasificar(df, UMBRAL_IMP, UMBRAL_CNLS) %>%
  filter(Set == "positivo", Status %in% c("TP", "FN*", "FN_L")) %>%
  mutate(Status = factor(Status, levels = c("TP", "FN*", "FN_L"))) %>%
  group_by(Model, Status, Match) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Model, Status) %>%
  mutate(Porcentaje = n / sum(n) * 100)

p_diag <- ggplot(df_diag, aes(x = Status, y = Porcentaje, fill = Match)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = paste0(round(Porcentaje, 1), "%")),
            position = position_stack(vjust = 0.5),
            color = "white", size = TFM_GEOM, fontface = "bold") +
  facet_wrap(~Model) +
  scale_fill_manual(values = c("YES" = "#4393c3", "NO" = "#f4a582")) +
  labs(title = "Localizacion de Motivos: Analisis de TP y Falsos Negativos",
       y = "Porcentaje dentro de la categoria (%)",
       x = "Categoria", fill = "Match UniProt") +
  theme_light(base_size = 14) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "gray90"))
ggsave(file.path(OUT_DIR, "figura_18_diagnostico_localizacion.png"),
       p_diag, width = 14, height = 8, dpi = FIG_DPI, bg = "white")

# --- Tabla 3: ensemble M1 vs M5 (union e interseccion) ---------------------------
get_metrics_strict <- function(call_col, match_col, data) {
  tp <- sum(data$Set == "positivo" & call_col & match_col == 1)
  fn <- sum(data$Set == "positivo") - tp
  tn <- sum(data$Set == "negativo" & !call_col)
  fp <- sum(data$Set == "negativo" & call_col)

  sens <- if (tp + fn > 0) tp / (tp + fn) * 100 else 0
  spec <- if (tn + fp > 0) tn / (tn + fp) * 100 else 0
  prec <- if (tp + fp > 0) tp / (tp + fp) * 100 else 0
  f1   <- if (prec + sens > 0) 2 * prec * sens / (prec + sens) else 0
  youden <- sens + spec - 100

  c(Sens_Estricta = sens, Spec = spec, F1_Score = f1, Youden = youden)
}

# Une predicciones de dos modelos por proteina (pivot_wider + sym)
analyze_ensemble_pair <- function(model_a, model_b, datos, u_imp = 0.6, u_cnls = 0.35) {
  df_pair <- datos %>%
    filter(Model %in% c(model_a, model_b)) %>%
    select(Protein_ID, Set, Model, ImpProb, cNLSProb, Match) %>%
    pivot_wider(
      names_from = Model,
      values_from = c(ImpProb, cNLSProb, Match),
      names_glue = "{.value}_{Model}"
    )

  col_imp_a   <- paste0("ImpProb_", model_a)
  col_imp_b   <- paste0("ImpProb_", model_b)
  col_cnls_a  <- paste0("cNLSProb_", model_a)
  col_cnls_b  <- paste0("cNLSProb_", model_b)
  col_match_a <- paste0("Match_", model_a)
  col_match_b <- paste0("Match_", model_b)

  df_ens <- df_pair %>%
    mutate(
      Call_A = (!!sym(col_imp_a) >= u_imp & !!sym(col_cnls_a) >= u_cnls),
      Call_B = (!!sym(col_imp_b) >= u_imp & !!sym(col_cnls_b) >= u_cnls),
      Match_A_bin = as.integer(!!sym(col_match_a) == "YES"),
      Match_B_bin = as.integer(!!sym(col_match_b) == "YES"),
      Union_Call = Call_A | Call_B,
      Union_Match_bin = as.integer(Match_A_bin == 1 | Match_B_bin == 1),
      Inter_Call = Call_A & Call_B,
      Inter_Match_bin = as.integer(Match_A_bin == 1 & Match_B_bin == 1)
    )

  combos <- list(
    list(nombre = model_a,        call = df_ens$Call_A,      match = df_ens$Match_A_bin),
    list(nombre = model_b,        call = df_ens$Call_B,      match = df_ens$Match_B_bin),
    list(nombre = "Union",        call = df_ens$Union_Call,  match = df_ens$Union_Match_bin),
    list(nombre = "Interseccion", call = df_ens$Inter_Call,  match = df_ens$Inter_Match_bin)
  )

  resultado <- data.frame()
  for (combo in combos) {
    m <- get_metrics_strict(combo$call, combo$match, df_ens)
    resultado <- rbind(resultado, data.frame(
      Modelo = combo$nombre,
      Umbral_Imp = u_imp,
      Umbral_cNLS = u_cnls,
      Sensibilidad_Estricta = round(m["Sens_Estricta"], 2),
      Especificidad = round(m["Spec"], 2),
      F1_Score = round(m["F1_Score"], 2),
      Indice_Youden = round(m["Youden"], 2)
    ))
  }
  resultado %>% arrange(desc(Indice_Youden))
}

tabla_3 <- analyze_ensemble_pair("M1", "M5", df)
write.csv(tabla_3, file.path(OUT_DIR, "tabla_3_ensemble_M1_M5.csv"), row.names = FALSE)
cat("\n--- Tabla 3: ensemble M1 vs M5 ---\n")
print(tabla_3, row.names = FALSE)

cat("\nAnalisis NLS completado. Salidas en:", OUT_DIR, "\n")
