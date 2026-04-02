#### Prunning Consensus Tree ######

install.packages("ape")
install.packages("seqinr")

library(ape)

# árvore filogenética consenso de mamíferos  
reference_tree_file <- "/home/maiara/Documents/02-2025-mestrado/data_processing/trees/reference_tree/MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre"

# espécies que ficarão na árvore referência  
distancias = read.csv('/home/maiara/Documents/Mestrado/distances_tables/allSpecies_resultsDistances.csv')
homo_distancias = distancias %>% filter(Species1 == "homo_sapiens" | Species2 == "homo_sapiens")
all_species <- tolower(homo_distancias$Species1)
all_species <- sub("^([^_]+_[^_]+).*", "\\1", all_species)

# corrige sinônimos
all_species[all_species == "notamacropus_eugenii"] <- "macropus_eugenii"
all_species[all_species == "physeter_catodon"]     <- "physeter_macrocephalus"
all_species[all_species == "cricetulus_griseus"]   <- "cricetulus_barabensis"
all_species[all_species == "carlito_syrichta"]   <- "tarsius_syrichta"
all_species[all_species == "equus_asinus"]   <- "equus_africanus"
all_species[all_species == "nannospalax_galili"]   <- "spalax_ehrenbergi"
all_species[all_species == "cebus_imitator"]   <- "cebus_capucinus"
all_species[all_species == "cervus_hanglu"]   <- "cervus_elaphus"
all_species["homo_sapiens"] <- "homo_sapiens"


# lendo a árvore e ajustando o nome para ficar igual ao nome das espécies 
big_tree <- read.nexus(reference_tree_file)

big_tree$tip.label[1:90]

# 1) Remover leading underscore (caso exista)
big_tree$tip.label <- gsub("^_+", "", big_tree$tip.label)

# 2) Converter tudo para minúsculas
big_tree$tip.label <- tolower(big_tree$tip.label)

# 3) Manter só genus_species (corta depois do segundo sublinhado)
# Esta expressão regular captura a parte antes do 3º underscore
big_tree$tip.label <- sub("^([^_]+_[^_]+).*", "\\1", big_tree$tip.label)

# checando 
common_sp <- intersect(big_tree$tip.label, all_species)

length(common_sp) 

ref_tree <- keep.tip(big_tree, common_sp)

write.tree(ref_tree, file="/home/maiara/Documents/02-2025-mestrado/data_processing/trees/reference_tree/reference_tree_96.nwk")

ref_tree$tip.label

missing_sp <- setdiff(all_species, big_tree$tip.label)

print(missing_sp)

