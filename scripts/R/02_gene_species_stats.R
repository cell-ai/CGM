#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: 02_gene_species_stats.R
# Descrição:
#   Calcula estatísticas descritivas das distâncias genéticas por gene e por
#   espécie. Remove genes com menos de 40 observações válidas, computa
#   mediana/média/desvio padrão por gene e mediana global por espécie.
#   Identifica genes nos extremos da distribuição (percentis 15 e 85),
#   que serão usados como entrada para o enriquecimento funcional (Script 03).
#
# Input:
#   - data/processed/adjustedDistances.csv  (gerado pelo Script 01)
#
# Output (em data/processed/):
#   - medianas_genes.csv
#       Mediana da distância genética para cada gene (entre todas as espécies).
#       Utilizado no Script 04 para análise de vias Reactome.
#   - stats_genes.csv
#       Estatísticas completas por gene: mediana, média e desvio padrão.
#   - medianas_especies.csv
#       Mediana global da distância genética por espécie (Homo sapiens vs. cada
#       outra espécie), com número de genes avaliados.
#   - genes_p15.csv
#       Genes com mediana de distância ≤ percentil 15 (genes mais conservados).
#   - genes_p85.csv
#       Genes com mediana de distância ≥ percentil 85 (genes mais divergentes).
# ------------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

# ------------------------------ Paths -----------------------------------------

base_dir <- "~/Documents/materials"
data_dir <- file.path(base_dir, "data", "processed")

input_adjusted_path <- file.path(data_dir, "adjustedDistances.csv")

output_median_species_path <- file.path(data_dir, "medianas_especies.csv")
output_median_genes_path   <- file.path(data_dir, "medianas_genes.csv")
output_stats_genes_path    <- file.path(data_dir, "stats_genes.csv")

# ------------------------ Leitura dos dados -----------------------------------

filtered_distancias <- read.csv(input_adjusted_path)

# --------------------- Filtro: genes com >= 40 observações --------------------

# Separar parte de espécies e parte de genes
species_cols <- filtered_distancias[, 1:2]
gene_cols    <- filtered_distancias[, 3:ncol(filtered_distancias)]

# Manter apenas genes com pelo menos 40 valores não-NA
keep_genes <- colSums(!is.na(gene_cols)) >= 40
gene_cols_filtered <- gene_cols[, keep_genes, drop = FALSE]

filtered_data <- cbind(species_cols, gene_cols_filtered)

# ---------------------- Estatísticas por gene ---------------------------------

medianas_todos_genes <- sapply(filtered_data[, -(1:2), drop = FALSE],
                               median, na.rm = TRUE)
medias_todos_genes   <- sapply(filtered_data[, -(1:2), drop = FALSE],
                               mean,   na.rm = TRUE)
desvios_todos_genes  <- sapply(filtered_data[, -(1:2), drop = FALSE],
                               sd,     na.rm = TRUE)

# Data frame com mediana por gene
medianas_df <- data.frame(
  Symbol          = names(medianas_todos_genes),
  Median_Distance = as.numeric(medianas_todos_genes),
  row.names       = NULL
)

# Data frame com estatísticas completas por gene
stats_df <- data.frame(
  Symbol          = names(medianas_todos_genes),
  Median_Distance = as.numeric(medianas_todos_genes),
  Mean_Distance   = as.numeric(medias_todos_genes),
  Std_Deviation   = as.numeric(desvios_todos_genes),
  row.names       = NULL
)

# ---------------------- Mediana por espécie -----------------------------------

# Aqui cada linha é um par (Species1, Species2, gene)
# Vamos "derreter" os dados para formato longo e agrupar por Species1
medianas_por_especie <- filtered_data %>%
  gather(key = "gene", value = "distance", -Species1, -Species2) %>%
  group_by(Species1) %>%
  summarise(
    median_distance = median(distance, na.rm = TRUE),
    count_genes     = n_distinct(gene[!is.na(distance)])
  ) %>%
  ungroup()

# Deixar nomes de espécies mais legíveis (substitui "_" por espaço e capitaliza)
medianas_por_especie <- medianas_por_especie %>%
  mutate(
    Species1 = gsub("_", " ", Species1),
    Species1 = sapply(
      Species1,
      function(name) paste0(toupper(substring(name, 1, 1)), substring(name, 2))
    )
  ) %>%
  arrange(median_distance)

colnames(medianas_por_especie) <- c("Species", "Median_Genetic_Distance", "Number_of_Genes")

# -------------------------- Salvar saídas -------------------------------------

write.csv(medianas_por_especie, output_median_species_path, row.names = FALSE)
write.csv(medianas_df,          output_median_genes_path,   row.names = FALSE)
write.csv(stats_df,             output_stats_genes_path,    row.names = FALSE)

# ------------------ Estatísticas descritivas globais --------------------------

summary_stats <- medianas_df %>%
  summarise(
    Media        = mean(Median_Distance,   na.rm = TRUE),
    Mediana      = median(Median_Distance, na.rm = TRUE),
    DesvioPadrao = sd(Median_Distance,     na.rm = TRUE),
    IQR          = IQR(Median_Distance,    na.rm = TRUE),
    Minimo       = min(Median_Distance,    na.rm = TRUE),
    Maximo       = max(Median_Distance,    na.rm = TRUE)
  )

cat("\nResumo das medianas de distância por gene:\n")
print(summary_stats)

# ----------------- Percentis 15 e 85 da distribuição --------------------------

p15 <- quantile(medianas_df$Median_Distance, 0.15, na.rm = TRUE)
p85 <- quantile(medianas_df$Median_Distance, 0.85, na.rm = TRUE)

cat("\n15º percentil:", p15, "\n")
cat("85º percentil:", p85, "\n")

# --------------------- Histograma das medianas --------------------------------

hist_plot <- ggplot(medianas_df, aes(x = Median_Distance)) +
  geom_histogram(bins = 60, fill = "lightblue", color = "black") +
  geom_vline(xintercept = p15, color = "darkblue", linetype = "dashed", size = 0.6) +
  geom_vline(xintercept = p85, color = "red",      linetype = "dashed", size = 0.6) +
  labs(
    title = "",
    x     = "Genetic Distance (Median per Gene)",
    y     = "Frequency"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.background   = element_rect(fill = "white", colour = "white"),
    legend.position    = "none"
  )

print(hist_plot)

# ------------------ Genes nos extremos da distribuição ------------------------

genes_p15 <- medianas_df %>%
  filter(Median_Distance <= p15)

genes_p85 <- medianas_df %>%
  filter(Median_Distance >= p85)

cat("\nNúmero de genes no <= 15º percentil:", nrow(genes_p15), "\n")
cat("Número de genes no >= 85º percentil:",  nrow(genes_p85), "\n")

# write.csv(genes_p15, file.path(data_dir, "genes_p15.csv"), row.names = FALSE)
# write.csv(genes_p85, file.path(data_dir, "genes_p85.csv"), row.names = FALSE)

# --------------------- QQplot das medianas ------------------------------------

# Remover NAs para o QQplot
medianas_sem_na <- medianas_df %>%
  filter(!is.na(Median_Distance))

qq_plot <- ggplot(medianas_sem_na, aes(sample = Median_Distance)) +
  stat_qq(size = 1) +
  stat_qq_line(colour = "red", linewidth = 0.5) +
  labs(
    title = "",
    x = "Theoretical Quantiles (Normal)",
    y = "Sample Quantiles (Median Genetic Distance)"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = "white")
  )

print(qq_plot)
