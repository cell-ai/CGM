#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: 05_gsea_species_reactome.R
# Descrição:
#   Executa Gene Set Enrichment Analysis (GSEA) por espécie usando as
#   distâncias genéticas como ranking. Para cada espécie, ordena os genes
#   pela distância genética (maior → menor) e roda fgseaMultilevel contra
#   o banco de vias Reactome 2024. Permite identificar vias imunológicas
#   enriquecidas entre os genes mais divergentes de cada espécie.
#
# Input:
#   - data/processed/adjustedDistances.csv  (gerado pelo Script 01)
#       Matriz de distâncias genéticas Homo sapiens vs. cada espécie.
#   - source/Reactome2024.txt
#       Arquivo GMT com gene sets do Reactome (vias biológicas).
#
# Output (em data/processed/PathwayAnalysis/):
#   - species_rankings/{espécie}_ranks.txt
#       Ranking de genes por distância genética (TSV, 2 colunas:
#       gene_id e distance), um arquivo por espécie.
#   - {espécie}/fgseaResults_{espécie}_ranks.txt.tsv
#       Resultados completos da GSEA por espécie (pathway, pval, padj,
#       NES, size, leadingEdge).
#   - {espécie}/fgseaResults_{espécie}_ranks.txt_signif.tsv
#       Resultados filtrados por significância (padj ≤ 0.05).
# ------------------------------------------------------------------------------

library(dplyr)
library(data.table)
library(fgsea)

# ------------------------------ Paths -----------------------------------------

base_dir   <- "~/Documents/materials"
data_dir   <- file.path(base_dir, "data", "processed")

# Tabela de distâncias ajustadas (saída do Script 1)
distances_path <- file.path(data_dir, "adjustedDistances.csv")

# Diretório para rankings e resultados por espécie
species_base_dir <- file.path(data_dir, "PathwayAnalysis")
ranks_dir        <- file.path(species_base_dir, "species_rankings")

# Arquivo GMT com vias de interesse (Reactome, 2024)
gmt_file <- file.path(base_dir, "source", "Reactome2024.txt")

# Garantir diretórios
dir.create(species_base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ranks_dir,        recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Construir rankings de genes por espécie
# ------------------------------------------------------------------------------

cat("Lendo distâncias ajustadas em:\n  ", distances_path, "\n")

all_distancias <- read.csv(distances_path)

# OBS:
# - Aqui você está usando Species1 como "espécie de interesse" por linha.
# - Para cada linha (espécie1 vs espécie2), ordena os genes da maior para a menor distância.

cat("Número de linhas (pares de espécies) em 'all_distancias': ", nrow(all_distancias), "\n")
cat("Criando rankings por linha / espécie (baseado em Species1)...\n")

for (i in seq_len(nrow(all_distancias))) {
  species_name <- all_distancias$Species1[i]
  species_data <- all_distancias[i, -(1:2)]  # remove colunas Species1 e Species2
  
  # transpor para ficar um gene por linha
  species_data_t <- t(species_data)
  
  gene_distances <- data.frame(
    gene_id  = rownames(species_data_t),
    distance = as.numeric(species_data_t[, 1]),
    row.names = NULL
  )
  
  gene_distances_clean  <- gene_distances %>% filter(!is.na(distance))
  gene_distances_sorted <- gene_distances_clean %>% arrange(desc(distance))
  
  species_file_name <- file.path(ranks_dir, paste0(species_name, "_ranks.txt"))
  
  write.table(
    gene_distances_sorted,
    file      = species_file_name,
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE
  )
  
  cat("  Processed species:", species_name, "->", species_file_name, "\n")
}

cat("Criação dos arquivos de ranking concluída.\n\n")

# ------------------------------------------------------------------------------
# 2. Função de GSEA (fgseaMultilevel) por espécie
# ------------------------------------------------------------------------------

ssGSEA <- function(gmtfile,
                   fileranks,
                   Ptype       = "padj",
                   pval_cutoff = 0.05,
                   output_dir  = "Species_Outputs") {
  
  # Função interna para ler GMT simples
  read.gmt <- function(fname) {
    gmt.lines <- readLines(fname)
    gmt.list  <- lapply(gmt.lines, function(x) unlist(strsplit(x, split = "\t")))
    gmt.names <- sapply(gmt.list, `[`, 1)
    gmt.genes <- lapply(gmt.list, function(x) x[3:length(x)])
    names(gmt.genes) <- gmt.names
    return(gmt.genes)
  }
  
  # Ler gene sets
  gmt <- read.gmt(gmtfile)
  
  # Ler ranks: espera 2 colunas (gene_id, distance)
  ranks_df <- read.delim(fileranks, header = TRUE, stringsAsFactors = FALSE)
  
  # Remover NA
  ranks_df <- ranks_df[complete.cases(ranks_df), ]
  
  # Vetor nomeado: valores = distance, nomes = gene_id
  rank_vals <- ranks_df$distance
  names(rank_vals) <- ranks_df$gene_id
  
  # Rodar fgsea
  fgseaRes <- fgseaMultilevel(
    pathways  = gmt,
    stats     = rank_vals,
    minSize   = 15,
    maxSize   = 2000,
    nproc     = 1,
    eps       = 0,
    scoreType = "pos"
  )
  
  fgseaRes <- as.data.frame(fgseaRes)
  
  # leadingEdge vira string única separada por vírgula
  fgseaRes$leadingEdge <- vapply(
    fgseaRes$leadingEdge,
    FUN.VALUE = character(1),
    FUN       = function(x) paste(x, collapse = ",")
  )
  
  # Criar output_dir se não existir
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Nome base do arquivo de saída (usa o nome do arquivo de ranks)
  base_out <- file.path(output_dir, paste0("fgseaResults_", basename(fileranks)))
  
  # Resultados completos
  write.table(
    fgseaRes,
    file      = paste0(base_out, ".tsv"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  
  # Filtrar significativos pelo Ptype escolhido (padj por padrão)
  if (!Ptype %in% colnames(fgseaRes)) {
    warning("Coluna ", Ptype, " não encontrada em fgseaRes. Mantendo todos os resultados.")
    fgseaRes_signif <- fgseaRes
  } else {
    fgseaRes_signif <- subset(fgseaRes, fgseaRes[[Ptype]] <= pval_cutoff)
  }
  
  write.table(
    fgseaRes_signif,
    file      = paste0(base_out, "_signif.tsv"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  
  cat("FGSEA finalizado para", fileranks, "→ resultados em:", output_dir, "\n")
  invisible(fgseaRes)
}

# ------------------------------------------------------------------------------
# 3. Loop sobre todas as espécies (arquivos de ranking)
# ------------------------------------------------------------------------------

cat("Rodando GSEA para cada arquivo de ranking em:\n  ", ranks_dir, "\n")

dir_outputs <- species_base_dir
rank_files  <- list.files(ranks_dir, pattern = "_ranks.txt$", full.names = TRUE)

for (f in rank_files) {
  species_name <- gsub("_ranks.txt", "", basename(f))
  
  # Subdiretório por espécie
  out_sp_dir <- file.path(dir_outputs, species_name)
  
  ssGSEA(
    gmtfile     = gmt_file,
    fileranks   = f,
    Ptype       = "padj",
    pval_cutoff = 0.05,
    output_dir  = out_sp_dir
  )
  
  cat("Completed GSEA for species:", species_name, "\n\n")
}
