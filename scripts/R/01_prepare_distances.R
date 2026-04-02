#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: 01_prepare_distances.R
# Descrição:
#   Prepara a matriz de distâncias genéticas brutas para as análises seguintes.
#   Filtra pares envolvendo Homo sapiens, remove espécies com baixa cobertura
#   de genes, exclui genes com alinhamentos problemáticos (MUC16, KIAA1958)
#   e aplica um limiar de distância (valores > threshold → NA).
#
# Input:
#   - data/raw/allSpecies_resultsDistances.csv
#       Tabela de distâncias genéticas entre pares de espécies (gene a gene)
#
# Output (em data/processed/):
#   - adjustedDistances.csv
#       Matriz de distâncias gene a gene filtrada e ajustada, com pares
#       Homo sapiens vs. demais espécies. Valores acima do limiar substituídos
#       por NA. Utilizada como entrada para os scripts 02–05.
#   - changeLog.csv
#       Log de todas as substituições realizadas pelo limiar, incluindo
#       espécie, gene e valor original da distância.
# ------------------------------------------------------------------------------

library(dplyr)

# ------------------------------ Paths -----------------------------------------

# Diretório base do projeto (ajuste se necessário)
base_dir <- "~/Documents/materials"

input_distances_path <- file.path(base_dir, "data", "raw", "allSpecies_resultsDistances.csv")

output_dir <- file.path(base_dir, "data", "processed")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

output_adjusted_path <- file.path(output_dir, "adjustedDistances.csv")
output_log_path      <- file.path(output_dir, "changeLog.csv")

# -------------------------- Leitura dos dados ---------------------------------

distancias <- read.csv(input_distances_path)

# Manter apenas pares envolvendo Homo sapiens
homo_distancias <- distancias %>%
  filter(Species1 == "homo_sapiens" | Species2 == "homo_sapiens")

# ---------------------- Remoção de espécies  -----------------------------

# Espécies com poucos genes / baixa cobertura
species_to_remove <- c(
  "choloepus_hoffmanni", "sorex_araneus",  "tupaia_belangeri",
  "erinaceus_europaeus", "echinops_telfairi", "vicugna_pacos",
  "notamacropus_eugenii", "procavia_capensis", "ochotona_princeps"
)

filtered_distancias <- homo_distancias %>%
  filter(!(Species1 %in% species_to_remove | Species2 %in% species_to_remove))

# ---------------------- Remoção de genes outliers -----------------------------

# Removendo colunas de genes com alinhamentos problemáticos
filtered_distancias <- filtered_distancias %>%
  select(-MUC16, -KIAA1958)

# ---------------- Função: aplicar limiar e registrar log ----------------------

apply_threshold_and_log <- function(data, threshold = 10) {
  # Data frame para armazenar log de alterações
  log_entries <- data.frame(
    Species1 = character(),
    Species2 = character(),
    Gene     = character(),
    Distance = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Identificar colunas numéricas (genes)
  numeric_cols <- vapply(data, is.numeric, logical(1L))
  numeric_col_names <- names(data)[numeric_cols]
  
  # Aplicar limiar e registrar mudanças
  for (gene_name in numeric_col_names) {
    high_values_idx <- which(data[[gene_name]] > threshold)
    
    if (length(high_values_idx) > 0) {
      log_entries <- rbind(
        log_entries,
        data.frame(
          Species1 = data$Species1[high_values_idx],
          Species2 = data$Species2[high_values_idx],
          Gene     = gene_name,
          Distance = data[[gene_name]][high_values_idx],
          stringsAsFactors = FALSE
        )
      )
      
      # Substituir valores acima do limiar por NA
      data[[gene_name]][high_values_idx] <- NA
    }
  }
  
  # Retornar lista com dados modificados e log
  list(data = data, log = log_entries)
}

# ------------------------- Aplicar função e salvar ----------------------------

result <- apply_threshold_and_log(filtered_distancias, threshold = 10)

# Salvar dados ajustados e log de alterações
write.csv(result$data, output_adjusted_path, row.names = FALSE)
write.csv(result$log,  output_log_path,      row.names = FALSE)

# Mostrar um resumo do log (opcional)
print(head(result$log))
cat("\nTotal de entradas modificadas:", nrow(result$log), "\n")
