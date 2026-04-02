##############################################################################
# Genetic Distance Project 
# Script for checking quality of pruned trees 
# 
# This script loads a single pruned tree and its corresponding alignment
# cleans up the alignment sequence names in the same way they were cleaned
# during pruning, and compares to the tip labels of the tree
#
# 02/2025 - Maiara M  
##############################################################################

# load required packages
library(ape)
library(seqinr)

# base directory
base_dir <- "/home/maiara/Documents/02-2025-mestrado"

# path for pruned tree
output_tree_dir <- file.path(base_dir, "data_processing", "trees", "pruned_trees")
# ex: checking gene CXCL1
tree <- file.path(output_tree_dir, "CXCL13_pruned.nwk")

# path for alignment file 
fasta_dir <- file.path(base_dir, "data_processing", "fastas", "aligned_fasta", "teste")
alignment <- file.path(fasta_dir, "filtered_CXCL13-aligned-trim.fasta")

##############################################################################
# Step 1: read pruned tree and alignment file 

test_tree <- read.tree(tree)
aln <- read.fasta(alignment)

cat("Tree contains", length(test_tree$tip.label), "species.\n")
cat("Alignment has", length(aln), "sequences.\n\n")

##############################################################################
# Step 2: adjusting alignment names

aln_names <- sapply(aln, function(x) attr(x, "name"))

clean_name <- function(original_name) {
  
  # split in |
  parts <- strsplit(original_name, "\\|")[[1]]
  if (length(parts) == 2) {
    nm <- parts[2]
  } else {
    nm <- original_name
  }
  
  # adjusting names 
  nm <- tolower(nm)
  nm <- sub("^([^_]+_[^_]+).*", "\\1", nm)
  
  # Replacing for synonyms
  if (nm == "notamacropus_eugenii")     nm <- "macropus_eugenii"
  if (nm == "physeter_catodon")         nm <- "physeter_macrocephalus"
  if (nm == "cricetulus_griseus")       nm <- "cricetulus_barabensis"
  if (nm == "carlito_syrichta")         nm <- "tarsius_syrichta"
  if (nm == "equus_asinus")             nm <- "equus_africanus"
  if (nm == "nannospalax_galili")       nm <- "spalax_ehrenbergi"
  if (nm == "cebus_imitator")           nm <- "cebus_capucinus"
  if (nm == "cervus_hanglu")            nm <- "cervus_elaphus"
  
  return(nm)
}

aln_names_clean <- sapply(aln_names, clean_name)

cat("Some cleaned names: \n")
print(head(aln_names_clean))

###############################################################################
# Step 3: Comparing names in the alignment and tree 

tree_tips <- test_tree$tip.label

cat("\nAlignmend names (sorted):\n")
print(sort(unique(aln_names_clean)))

cat("\nNames in tree (sorted):\n")
print(sort(unique(tree_tips)))

common_names <- intersect(aln_names_clean, tree_tips)
cat("Names in common:", length(common_names), "\n")

missing_in_tree <- setdiff(aln_names_clean, tree_tips)
if (length(missing_in_tree) == 0) {
  cat("Any species is missing in the tree.\n")
} else {
  cat("Species present in the alignment but missing in the tree:\n")
  cat(paste0(missing_in_tree, collapse = ", "), "\n")
}

# Quem está na árvore, mas não no alinhamento
missing_in_aln <- setdiff(tree_tips, aln_names_clean)
if (length(missing_in_aln) == 0) {
  cat("Any species is missing in the alignment.\n")
} else {
  cat("Species present in the tree but missing in the alignment:\n")
  cat(paste0(missing_in_aln, collapse = ", "), "\n")
}

###############################################################################
# Overview 

cat("\n=== 5) Resumo final ===\n")
cat(" - Número de sequências no alinhamento (bruto):", length(aln), "\n")
cat(" - Número de tips na árvore podada:", length(test_tree$tip.label), "\n")
cat(" - Nomes em comum (após limpeza):", length(common_names), "\n")

if (length(common_names) == length(aln_names_clean) && length(common_names) == length(test_tree$tip.label)) {
  cat("Tudo certo: todos os nomes do alinhamento batem com a árvore!\n")
} else {
  cat("Atenção: há diferenças entre o alinhamento e a árvore.\n")
  cat("Verifique as listas de faltantes acima.\n")
}

cat("\n=== Fim da checagem ===\n")