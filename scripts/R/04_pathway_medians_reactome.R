#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: 04_pathway_medians_reactome.R
# Descrição:
#   Avalia se vias imunológicas do Reactome possuem distâncias genéticas
#   significativamente diferentes do restante do genoma. Para cada via,
#   calcula a mediana das distâncias dos genes-membro, intervalo de
#   confiança por bootstrap e teste de Wilcoxon (via vs. resto).
#
# Input:
#   - source/Immune_pathways_Reactome.csv
#       Lista de vias imunológicas Reactome com seus genes associados.
#   - data/processed/medianas_genes.csv  (gerado pelo Script 02)
#       Mediana de distância genética por gene.
#
# Output (em data/processed/PathwayAnalysis/):
#   - reactome_genesets_with_medians.csv
#       Tabela com uma linha por via contendo: mediana da distância,
#       intervalo de confiança (bootstrap), p-valor (Wilcoxon),
#       p-valor ajustado (BH), cobertura de genes e lista de genes
#       presentes na análise.
# ------------------------------------------------------------------------------

library(data.table)

base_dir <- "~/Documents/materials"
data_dir <- file.path(base_dir, "data", "processed")

genesets_path <- file.path(base_dir, "source", "Immune_pathways_Reactome.csv")
medianas_path <- file.path(data_dir, "medianas_genes.csv")
output_path   <- file.path(data_dir, "PathwayAnalysis", "reactome_genesets_with_medians.csv")

dir.create(file.path(data_dir, "PathwayAnalysis"), recursive = TRUE, showWarnings = FALSE)

# -------------------- IC da mediana (bootstrap) ------------------------------

calcular_ic_mediana <- function(x, n_boot = 1000, conf_level = 0.95) {
  x <- x[is.finite(x) & !is.na(x)]
  n <- length(x)
  if (n < 2) return(c(NA_real_, NA_real_))
  boot_meds <- replicate(n_boot, median(sample(x, size = n, replace = TRUE)))
  alpha <- (1 - conf_level) / 2
  quantile(boot_meds, probs = c(alpha, 1 - alpha), na.rm = TRUE)
}

# ------------------------ process_genes (ajustada) ---------------------------

process_genes <- function(x) {
  # x: string do tipo "['GENE1', 'GENE2']" ou "GENE1, GENE2"
  x <- gsub("\\[|\\]|'|\"", "", x)              # remove colchetes/aspas
  parts <- unlist(strsplit(x, "[,;]"))         # quebra por vírgula ou ponto-e-vírgula
  parts <- trimws(parts)                       # tira espaços
  parts <- parts[parts != ""]                  # remove vazios
  return(parts)
}

# --------------------------- ler dados ---------------------------------------

genesets_df <- fread(genesets_path)
medianas_df <- fread(medianas_path)

# aplicar processamento nos genes dos genesets
genesets_df$Genes <- lapply(genesets_df$Genes, process_genes)

# maiúsculas + trim
genesets_df$Genes <- lapply(genesets_df$Genes, toupper)
medianas_df$Symbol <- toupper(medianas_df$Symbol)

genesets_df$Genes <- lapply(genesets_df$Genes, trimws)
medianas_df$Symbol <- trimws(medianas_df$Symbol)

# tirar NAs de mediana
medianas_df <- medianas_df[!is.na(Median_Distance)]

# ---------------------- função calcular_resultados ---------------------------

calcular_resultados <- function(genesets_df, medianas_df) {
  resultados <- data.table(
    Pathway        = character(),
    Mediana_Via    = numeric(),
    IC_Lower       = numeric(),
    IC_Upper       = numeric(),
    PValue         = numeric(),
    Cobertura      = numeric(),
    Genes_Presente = list()
  )
  
  medianas_df <- medianas_df[!is.na(Median_Distance)]
  
  for (i in 1:nrow(genesets_df)) {
    pathway       <- genesets_df$Pathway[i]
    pathway_genes <- genesets_df$Genes[[i]]
    
    total_genes_via     <- length(pathway_genes)
    Genes_Presente      <- intersect(pathway_genes, medianas_df$Symbol)
    num_genes_disponiveis <- length(Genes_Presente)
    cobertura           <- if (total_genes_via > 0) num_genes_disponiveis / total_genes_via else NA_real_
    
    if (num_genes_disponiveis >= 2) {
      distancias_via    <- medianas_df[Symbol %in% Genes_Presente, Median_Distance]
      distancias_outros <- medianas_df[!Symbol %in% Genes_Presente, Median_Distance]
      
      wilcox_result <- wilcox.test(distancias_via, distancias_outros, alternative = "two.sided")
      p_value       <- wilcox_result$p.value
      mediana_via   <- median(distancias_via, na.rm = TRUE)
      ic_mediana    <- calcular_ic_mediana(distancias_via)
      
      resultados <- rbind(
        resultados,
        data.table(
          Pathway        = pathway,
          Mediana_Via    = mediana_via,
          IC_Lower       = ic_mediana[1],
          IC_Upper       = ic_mediana[2],
          PValue         = p_value,
          Cobertura      = cobertura,
          Genes_Presente = list(Genes_Presente)
        ),
        use.names = TRUE,
        fill      = TRUE
      )
    }
  }
  
  if (nrow(resultados) > 0) {
    resultados[, PAdjusted := p.adjust(PValue, method = "BH")]
  } else {
    resultados[, PAdjusted := numeric(0)]
  }
  
  return(resultados)
}

# -------------------------- rodar e salvar -----------------------------------

resultados_all_genes <- calcular_resultados(genesets_df, medianas_df)
cat("Número de vias com >= 2 genes mapeados:", nrow(resultados_all_genes), "\n")

fwrite(resultados_all_genes, output_path)
cat("Resultados salvos em:\n  ", output_path, "\n")
