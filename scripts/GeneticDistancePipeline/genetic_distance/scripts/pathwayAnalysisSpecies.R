# Gene Set Enrichment Analysis using Genetic Distance information for each species of mammals to humans. 
# GSEA using fgsea 

library(tidyr)
library(stringr)
library(dplyr)
library(data.table)
library(clusterProfiler)
library(fgsea)
library(ggplot2)

# Primeiro: criar uma lista de ranks dos genes para cada espécie, ordenada do maior para o menor

all_distancias <- read.csv("/home/maiara/Documents/Mestrado/Novas_analises/adjustedDistances.csv")

# lembrete: aqui, ao contrário do que foi feito no cálculo das medianas de cada gene, não foi feita a remoção de colunas de genes que tenham menos que 40 distâncias

species_list <- all_distancias$Species1
dir_ranks <- "/home/maiara/Documents/Mestrado/Novas_analises/pathwayAnalysisSpecies/species_rankings"


for (i in 1:nrow(all_distancias)) {
  species_name <- all_distancias$Species1[i]
  species_data <- all_distancias[i, -(1:2)]
  species_data_t <- t(species_data)
  gene_distances <- data.frame(
    gene_id = rownames(species_data_t),
    distance = as.numeric(species_data_t[, 1])
  )
  
  gene_distances_clean <- gene_distances %>% filter(!is.na(distance))
  gene_distances_sorted <- gene_distances_clean %>% arrange(desc(distance))
  
  species_file_name <- paste0("/home/maiara/Documents/Mestrado/Novas_analises/pathwayAnalysisSpecies/species_rankings/",species_name,"_ranks.txt")
  
  write.table(
    gene_distances_sorted,
    file = species_file_name, 
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  cat("Processed species:", species_name,"\n")
}


################################################################
################################################################

ssGSEA <- function(gmtfile, fileranks, Ptype = "padj", pval_cutoff = 0.05, output_dir = "Species_Outputs") {
  
  # Função interna para ler GMT
  read.gmt <- function(fname){
    gmt.lines <- readLines(fname)
    gmt.list  <- lapply(gmt.lines, function(x) unlist(strsplit(x, split="\t")))
    gmt.names <- sapply(gmt.list, '[', 1)
    gmt.genes <- lapply(gmt.list, function(x) x[3:length(x)])
    names(gmt.genes) <- gmt.names
    return(gmt.genes)
  }
  
  # Ler gene sets
  gmt <- read.gmt(gmtfile)
  
  # Ler ranks => supostamente 2 colunas: geneid e distance
  ranks_df <- read.delim(fileranks, header = TRUE, stringsAsFactors = FALSE)
  
  # Remover NA
  ranks_df <- ranks_df[complete.cases(ranks_df), ]
  
  # Transformar em vetor named
  rank_vals <- ranks_df$distance
  names(rank_vals) <- ranks_df$gene_id
  
  # Rodar fgsea
  fgseaRes <- fgseaMultilevel(pathways = gmt,
                              stats    = rank_vals,
                              minSize  = 15,
                              maxSize  = 2000,
                              nproc    = 1,
                              eps      = 0, 
                              scoreType="pos")
  
  fgseaRes <- as.data.frame(fgseaRes)
  
  fgseaRes$leadingEdge <- vapply(
    fgseaRes$leadingEdge, 
    FUN = function(x) paste(x, collapse=","),
    FUN.VALUE = character(1)
  )
  
  # Criar output_dir se não existir
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Nome base do arquivo
  base_out <- file.path(output_dir, paste0("fgseaResults_", basename(fileranks)))
  
  write.table(
    fgseaRes,
    file      = paste0(base_out, ".tsv"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  
  # Filtrar por pval_cutoff
  fgseaRes_signif <- subset(fgseaRes, fgseaRes[[Ptype]] <= pval_cutoff)
  
  # significativos
  write.table(
    fgseaRes_signif,
    file      = paste0(base_out, "_signif.tsv"),
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )
  
  cat("FGSEA finalizado para", fileranks, " resultados em:", output_dir, "\n")
  return(invisible(fgseaRes))
}

##################################################
##################################################
##################################################

dir_outputs <- "/home/maiara/Documents/Mestrado/Novas_analises/pathwayAnalysisSpecies/"

rank_files <- list.files(dir_ranks, pattern="_ranks.txt$", full.names=TRUE)

for (f in rank_files) {
  species_name <- gsub("_ranks.txt","", basename(f))
  
  # criar subdir
  out_sp_dir <- file.path(dir_outputs, species_name)
  
  # rodar
  ssGSEA(
    gmtfile     = "/home/maiara/Documents/Mestrado/Novas_analises/pathwayAnalysisSpecies/Reactome2024.txt",
    fileranks   = f,
    Ptype       = "padj",
    pval_cutoff = 0.05,
    output_dir  = out_sp_dir
  )
  
  cat("Completed GSEA for species:", species_name, "\n\n")
}
