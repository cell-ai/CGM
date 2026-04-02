#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: 03_enrichment_go_kegg.R
# Descrição:
#   Realiza análise de enriquecimento funcional (ORA) nos conjuntos de genes
#   conservados (≤ p15) e divergentes (≥ p85) identificados no Script 02.
#   Utiliza Gene Ontology (BP, CC, MF) via clusterProfiler e KEGG pathways.
#
# Input (em data/processed/, gerados pelo Script 02):
#   - genes_p15.csv   — genes com distância ≤ percentil 15 (conservados)
#   - genes_p85.csv   — genes com distância ≥ percentil 85 (divergentes)
#
# Output (em data/processed/enrichment_results/):
#   Tabelas de enriquecimento (CSV):
#   - go_bp_low.csv / go_bp_high.csv   — GO Biological Process
#   - go_cc_low.csv / go_cc_high.csv   — GO Cellular Component
#   - go_mf_low.csv / go_mf_high.csv   — GO Molecular Function
#   - kegg_enrichment_low.csv / kegg_enrichment_high.csv — KEGG pathways
#
#   Figuras (em data/processed/enrichment_results/figures/):
#   - go_enrichment_low.png / go_enrichment_high.png   — barplot facetado GO
#   - kegg_barplot_low.png / kegg_barplot_high.png     — barplot KEGG
#   - kegg_dotplot_low.png / kegg_dotplot_high.png     — dotplot KEGG
# ------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)

# ------------------------------ Paths -----------------------------------------

base_dir <- "~/Documents/materials"
data_dir <- file.path(base_dir, "data", "processed")

# Arquivos de entrada (gerados no Script 2)
genes_p15_path  <- file.path(data_dir, "genes_p15.csv")   # genes com mediana <= p15
genes_p85_path  <- file.path(data_dir, "genes_p85.csv")   # genes com mediana >= p85

# Diretório para salvar tabelas e figuras
results_dir     <- file.path(data_dir, "enrichment_results")
figures_dir     <- file.path(results_dir, "figures")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------- Leitura dos genes ---------------------------------

genes_low_df  <- read.csv(genes_p15_path)
genes_high_df <- read.csv(genes_p85_path)

# Assumindo coluna "Symbol" com os símbolos dos genes
gene_lists <- list(
  low  = unique(genes_low_df$Symbol),
  high = unique(genes_high_df$Symbol)
)

# ------------------- Função auxiliar para GO enrichment -----------------------

run_go_enrichment <- function(gene_list, set_name) {
  message("\nRodando GO enrichment para conjunto: ", set_name,
          " (n = ", length(gene_list), " genes)")
  
  # GO: Biological Process
  egoBP <- enrichGO(
    gene          = gene_list,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.1,
    qvalueCutoff  = 0.05
  )
  
  # GO: Cellular Component
  egoCC <- enrichGO(
    gene          = gene_list,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.1,
    qvalueCutoff  = 0.05
  )
  
  # GO: Molecular Function
  egoMF <- enrichGO(
    gene          = gene_list,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.1,
    qvalueCutoff  = 0.05
  )
  
  # Converter para data.frame (tratando caso de não haver termos significativos)
  df_BP <- if (!is.null(egoBP) && nrow(as.data.frame(egoBP)) > 0) {
    as.data.frame(egoBP)
  } else {
    data.frame()
  }
  
  df_CC <- if (!is.null(egoCC) && nrow(as.data.frame(egoCC)) > 0) {
    as.data.frame(egoCC)
  } else {
    data.frame()
  }
  
  df_MF <- if (!is.null(egoMF) && nrow(as.data.frame(egoMF)) > 0) {
    as.data.frame(egoMF)
  } else {
    data.frame()
  }
  
  # Salvar tabelas completas
  if (nrow(df_BP) > 0) {
    write.csv(df_BP,
              file.path(results_dir, paste0("go_bp_", set_name, ".csv")),
              row.names = FALSE)
  }
  if (nrow(df_CC) > 0) {
    write.csv(df_CC,
              file.path(results_dir, paste0("go_cc_", set_name, ".csv")),
              row.names = FALSE)
  }
  if (nrow(df_MF) > 0) {
    write.csv(df_MF,
              file.path(results_dir, paste0("go_mf_", set_name, ".csv")),
              row.names = FALSE)
  }
  
  # Se não tiver nada significativo, já retorna aqui
  if (nrow(df_BP) == 0 && nrow(df_CC) == 0 && nrow(df_MF) == 0) {
    message("Nenhum termo GO significativo para o conjunto: ", set_name)
    return(invisible(NULL))
  }
  
  # Adicionar coluna de categoria para combinar
  if (nrow(df_BP) > 0) df_BP$Category <- "BP"
  if (nrow(df_CC) > 0) df_CC$Category <- "CC"
  if (nrow(df_MF) > 0) df_MF$Category <- "MF"
  
  combined_df <- dplyr::bind_rows(df_BP, df_CC, df_MF)
  
  # Ordenar por p.adjust dentro de cada categoria
  combined_df <- combined_df %>%
    group_by(Category) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    ungroup()
  
  # Opcional: selecionar top N por categoria para o gráfico
  combined_df <- combined_df %>%
    group_by(Category) %>%
    slice_head(n = 20) %>%
    ungroup()
  
  # Fator ordenado para descrição
  combined_df$Description <- factor(
    combined_df$Description,
    levels = rev(unique(combined_df$Description))
  )
  
  # Gráfico de barras facetado por categoria
  go_plot <- ggplot(combined_df,
                    aes(x = Description, y = -log10(p.adjust), fill = Category)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    facet_wrap(~Category, scales = "free", ncol = 1) +
    labs(
      x     = "GO Terms",
      y     = "-log10(p.adjust)",
      title = paste("GO Enrichment Analysis -", set_name)
    ) +
    theme_minimal() +
    theme(
      axis.text.y  = element_text(size = 7),
      plot.title   = element_text(hjust = 0.5)
    )
  
  print(go_plot)
  
  ggsave(
    filename = file.path(figures_dir, paste0("go_enrichment_", set_name, ".png")),
    plot     = go_plot,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
  
  invisible(list(egoBP = egoBP, egoCC = egoCC, egoMF = egoMF))
}

# ------------------- Função auxiliar para KEGG enrichment ---------------------

run_kegg_enrichment <- function(gene_list, set_name) {
  message("\nRodando KEGG enrichment para conjunto: ", set_name,
          " (n = ", length(gene_list), " genes)")
  
  # Converter símbolos para ENTREZID
  gene_df <- bitr(
    gene_list,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
  
  if (nrow(gene_df) == 0) {
    message("Nenhum gene convertido para ENTREZID no conjunto: ", set_name)
    return(invisible(NULL))
  }
  
  # Análise KEGG (mesma lógica do seu script)
  kegg_enrich <- enrichKEGG(
    gene          = gene_df$ENTREZID,
    organism      = "hsa",
    pvalueCutoff  = 0.01,         # igual ao seu original
    pAdjustMethod = "BH"
  )
  
  if (is.null(kegg_enrich) || nrow(as.data.frame(kegg_enrich)) == 0) {
    message("Nenhuma via KEGG significativa para o conjunto: ", set_name)
    return(invisible(NULL))
  }
  
  # Manter a mesma lógica: criar log_p_adjust DENTRO do objeto
  sorted_kegg <- kegg_enrich
  sorted_kegg@result <- sorted_kegg@result %>%
    arrange(p.adjust) %>%
    mutate(log_p_adjust = -log10(p.adjust))
  
  # (Opcional) salvar tabela completa
  kegg_df <- as.data.frame(sorted_kegg)
  write.csv(
    kegg_df,
    file.path(results_dir, paste0("kegg_enrichment_", set_name, ".csv")),
    row.names = FALSE
  )
  
  # Barplot ordenado pelo p.adjust (top 20 categorias)
  kegg_bar <- barplot(
    sorted_kegg,
    showCategory = 20,
    title = paste("KEGG Pathway Enrichment -", set_name)
  )
  
  ggsave(
    filename = file.path(figures_dir, paste0("kegg_barplot_", set_name, ".png")),
    plot     = kegg_bar,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
  
  # Dotplot usando log_p_adjust criado em @result
  kegg_dot <- dotplot(sorted_kegg, x = "log_p_adjust") +
    labs(
      x     = "-log10(p.adjust)",
      y     = "KEGG Pathways",
      title = paste("KEGG Pathway Enrichment (dotplot) -", set_name)
    )
  
  ggsave(
    filename = file.path(figures_dir, paste0("kegg_dotplot_", set_name, ".png")),
    plot     = kegg_dot,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
  
  print(kegg_bar)
  print(kegg_dot)
  
  invisible(sorted_kegg)
}


# ---------------------- Loop pelos dois conjuntos -----------------------------

for (set_name in names(gene_lists)) {
  gene_list <- gene_lists[[set_name]]
  
  # Remover NAs e strings vazias
  gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
  gene_list <- unique(gene_list)
  
  if (length(gene_list) < 2) {
    message("\nConjunto muito pequeno para enriquecimento: ", set_name)
    next
  }
  
  # GO
  run_go_enrichment(gene_list, set_name)
  
  # KEGG
  run_kegg_enrichment(gene_list, set_name)
}

message("\nAnálises de enriquecimento (GO + KEGG) concluídas para conjuntos low e high.\n")
