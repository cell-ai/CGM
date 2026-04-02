# Genetic Distance Project #

# load libraries
library(seqinr)
library(ape)
library(tools)

# base directory
base_dir <- "/home/maiara/Documents/02-2025-mestrado"

# path for reference tree
ref_tree_path <- file.path(base_dir, "data_processing", "trees", "reference_tree", "reference_tree_96.nwk")

# path for FASTAs files
fasta_dir <- file.path(base_dir, "data_processing", "fastas", "aligned_fasta", "teste")

# path for trees 
out_tree_dir <- file.path(base_dir, "data_processing", "trees", "pruned_trees")

# read reference tree
ref_tree <- read.tree(ref_tree_path)

# list FASTAs files in directory
fasta_files <- list.files(
  path = fasta_dir,
  pattern = "\\.fasta$",
  full.names = TRUE
)

# loop for each FASTA
for (ff in fasta_files) {
  
  # read alignments
  aln_seqs <- read.fasta(ff)
  seq_names <- sapply(aln_seqs, function(x) attr(x, "name"))
  
  # remove ID before "|"  
  seq_names_clean <- sapply(seq_names, function(x) {
    parts <- strsplit(x, "\\|")[[1]]
    if (length(parts) == 2) {
      return(parts[2])
    } else {
      return(x)
    }
  })
  
  # adjusting names to match tree species names
  seq_names_clean <- tolower(seq_names_clean)
  seq_names_clean <- sub("^([^_]+_[^_]+).*", "\\1", seq_names_clean)
  
  # replacing synonyms 
  seq_names_clean[seq_names_clean == "notamacropus_eugenii"] <- "macropus_eugenii"
  seq_names_clean[seq_names_clean == "physeter_catodon"]     <- "physeter_macrocephalus"
  seq_names_clean[seq_names_clean == "cricetulus_griseus"]   <- "cricetulus_barabensis"
  seq_names_clean[seq_names_clean == "carlito_syrichta"]   <- "tarsius_syrichta"
  seq_names_clean[seq_names_clean == "equus_asinus"]   <- "equus_africanus"
  seq_names_clean[seq_names_clean == "nannospalax_galili"]   <- "spalax_ehrenbergi"
  seq_names_clean[seq_names_clean == "cebus_imitator"]   <- "cebus_capucinus"
  seq_names_clean[seq_names_clean == "cervus_hanglu"]   <- "cervus_elaphus"
  
  # intersection
  common_tips <- intersect(ref_tree$tip.label, seq_names_clean)
  
  # creating pruned tree
  pruned_tree <- keep.tip(ref_tree, common_tips)

  ### adjusting the file name 
  
  # file path 
  file_stem <- file_path_sans_ext(basename(ff))
  # ex: filtered_CXCL1-aligned-trim
  
  file_stem <- sub("^filtered_", "", file_stem)
  # ex: CXCL1-aligned-trim
  
  file_stem <- sub("-aligned-trim$", "", file_stem)
  # ex: CXCL1
  
  outname <- paste0(file_stem, "_pruned.nwk")
  outpath <- file.path(out_tree_dir, outname)
  
  # saving tree
  write.tree(pruned_tree, outpath)
  
}