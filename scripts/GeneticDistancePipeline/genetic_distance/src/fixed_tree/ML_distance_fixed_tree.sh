#!/usr/bin/env bash

######################################################################################
############## Genetic Distance Project - Distance Fixed Tree Script #################
#
# Author: Maiara Moro 
#
# This script is used to calculate the genetic distance between two sequences using a fixed tree
# The input file is a fasta file containing the sequences and a nwk file containing the tree
# The output file will contain the genetic distance between all sequences 
# 
# Program for pair-wise ML distances: raxml (https://cme.h-its.org/exelixis/resource/download/NewManual.pdf)
# ML model parameters will be estimated on an starting species tree (nwk file)
# The algorithm will compute branch lengths between the sequences of an alignment using a ML estimate 
#
# # distanceFixedTree.sh
# Usage: bash ML_distance_fixex_tree.sh <pruned_trees_dir> <trimmed_alignments_dir> <output_dir>
# bash ML_distance_fixed_tree.sh ../../data_processing/trees/pruned_trees ../../data_processing/fastas/aligned_fasta/teste ../../data_processing/raxml_output_fixed_tree
######################################################################################

set -euo pipefail


PRUNED_TREES_DIR="$1"
TRIMMED_ALIGNMENTS_DIR="$2"
OUTPUT_DIR="$(realpath "$3")"

mkdir -p "$OUTPUT_DIR"

echo "===================================================="
echo "    RAxML Distance Calculation with Fixed Trees     "
echo "                                                    "
echo " Pruned trees dir: $PRUNED_TREES_DIR                "
echo " Trimmed alignments dir: $TRIMMED_ALIGNMENTS_DIR    " 
echo " Output dir: $OUTPUT_DIR                            "
echo "                                                    "
echo "===================================================="

# loop over all pruned trees
for TREE_FILE in "$PRUNED_TREES_DIR"/*_pruned.nwk; do

    [ -e "$TREE_FILE" ] || continue

    # get the tree name
    BASENAME=$(basename "$TREE_FILE" "_pruned.nwk")
    GENE="$BASENAME"

    # get the trimmed alignment
    ALN_FILE="$TRIMMED_ALIGNMENTS_DIR/filtered_${GENE}-aligned-trim.fasta"

    if [[ ! -f "$ALN_FILE" ]]; then
        echo "Alignment file not found: $ALN_FILE"
        continue
    fi

    echo "Processing gene $GENE"
    echo "  Tree file: $TREE_FILE"
    echo "  Alignment file: $ALN_FILE"

    # verify if the raxml file already exists
    EXPECTED_OUTPUT="$OUTPUT_DIR/RAxML_distances.${GENE}_matrix"

    # if yes, then skip this gene
    if [[ -f "$EXPECTED_OUTPUT" ]]; then
        echo "Output file already exists for gene $GENE. Skipping..."
        continue
    fi


    # adjust the alignment (making appropriate for the tree) and create a temporary file to store the cleaned alignment (appropiate for the tree) using the script clean_fasta.sh
    CLEAN_FASTA="${OUTPUT_DIR}/${GENE}_clean.fasta"
    bash clean_fasta.sh "$ALN_FILE" "$CLEAN_FASTA"

    JOB_NAME="${GENE}_matrix"

    # RaXML command:
    # -f x: compute pair-wise ML distances
    # -t: fixed tree
    # -m PROTGAMMAAUTO: select between protein models using the ML Criterion
    # -p 12345: random seed
    # -s: alignment file
    # -n: output file
    # -w: output directory

    # run raxml
    raxmlHPC -f x -p 12345 -s "$CLEAN_FASTA" -t "$TREE_FILE" -m PROTGAMMAAUTO -n "$JOB_NAME" -w "$OUTPUT_DIR"

    rm -f "$CLEAN_FASTA"

echo "RaxML finished successfully for gene $GENE"

done

echo "All done!"

