#!/usr/bin/env bash

################################################################################
#                   Genetic Distance Project
#
# Author: Maiara Moro
# 
# This script cleans a FASTA file to match species name in a phylogenetic tree
# - Removes prefixes (ENSEMBL IDs)
# - Applies synonym substitutions
# - Removes redundant sub-species names
# - Ensures only "genus_species" format (two names)
# - Removes canis_lupus_dingo to avoid duplication in RAxML
# 
# Usage: bash clean_fasta.sh <input_fasta> <output_clean_fasta>
#
################################################################################

set -euo pipefail

INPUT_FASTA="$1"
OUTPUT_CLEAN_FASTA="$2"

# Remove prefixes (ENSEMBL IDs)
# Substitute synonyms
# Remove redundant sub-species names

awk '
    BEGIN { OFS="\n"; keep=1 }   # Garante que a saída mantenha formatação correta
    /^>/ {
        gsub(/^>/, "", $0);                      # Remove ">"
        split($0, arr, "|");                     # Divide no "|"
        species = arr[length(arr)];              # Pega apenas o nome final
        
        # Aplicar substituições de sinônimos
        gsub("notamacropus_eugenii", "macropus_eugenii", species);
        gsub("physeter_catodon", "physeter_macrocephalus", species);
        gsub("cricetulus_griseus", "cricetulus_barabensis", species);
        gsub("carlito_syrichta", "tarsius_syrichta", species);
        gsub("equus_asinus", "equus_africanus", species);
        gsub("nannospalax_galili", "spalax_ehrenbergi", species);
        gsub("cebus_imitator", "cebus_capucinus", species);
        gsub("cervus_hanglu_yarkandensis", "cervus_hanglu", species);
        gsub("cervus_hanglu", "cervus_elaphus", species);

        if (species == "canis_lupus_dingo") {
            keep=0; next;
        } 

        split(species, parts, "_");
        if (length(parts) > 2) {
            species = parts[1] "_" parts[2];   # Mantém apenas genus_species
        }

        print ">" species;
        keep=1;
        next;
    } 
    
    keep { print $0 }   # Mantém apenas sequências que não foram removidas

' "$INPUT_FASTA" > "$OUTPUT_CLEAN_FASTA"

echo "FASTA cleaned: $OUTPUT_CLEAN_FASTA"