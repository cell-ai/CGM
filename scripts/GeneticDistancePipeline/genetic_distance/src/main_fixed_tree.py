#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#####################################################################################################
#                                   Genetic Distance Project                                        #
# Author: Maiara Moro                                                                               #    
#                                                                                                   #                                                      
# This pipeline automates the process of genetic data analysis. It is divided into 7 steps:         #
#                                                                                                   #                                           
# Based on a initial gene list of interest, the pipeline will:                                      #
#                                                                                                   #
# 1 - Fetch canonical orthologous sequences of mammals from ENSEMBL database using API REST         #
# 2 - Filter to keep only one-to-one orthologs                                                      # 
# 3 - Multiple sequence alignment using MAFFT                                                       #
# 4 - Alignment trimming using Trimal                                                               #
# 5 - Prune species phylogenetic tree based on the species in the alignment                         #
# 6 - Calculate pairwise genetic distances using ML estimates - RAxML                               #    
# 7 - Summarize results in a CSV table                                                              #  
#                                                                                                   #    
#####################################################################################################

# import libraries
import os
import logging
import subprocess
import sys

# import local modules
from fixed_tree.table_distances2 import create_all_distances_table

# initialize logger
logging.basicConfig(level=logging.INFO)

# directory paths
BASE_DIR = "/home/maiara/Documents/02-2025-mestrado" 
PRUNED_TREES_DIR = os.path.join(BASE_DIR, "data_processing/trees/pruned_trees")
ALIGNED_FASTA_DIR = os.path.join(BASE_DIR, "data_processing/fastas/aligned_fasta/teste")
RAXML_OUTPUT_DIR = os.path.join(BASE_DIR, "data_processing/raxml_output_fixed_tree")
DISTANCES_CSV_PATH = os.path.join(BASE_DIR, "distance_results/geneticDistance2/Fixed_tree_allSpecies_resultsDistances.csv")

# create directories
os.makedirs(RAXML_OUTPUT_DIR, exist_ok=True)

def run_raxml_fixed_tree():

    """
    Runs RAxML using pruned species trees and aligned FASTA files to compute genetic distances.

    """
    logging.info("Starting RAxML pairwise genetic distance calculation...")

    raxml_script = os.path.join(BASE_DIR, "src/fixed_tree/ML_distance_fixed_tree.sh")

    try:
        subprocess.run(
            ["bash", raxml_script, PRUNED_TREES_DIR, ALIGNED_FASTA_DIR, RAXML_OUTPUT_DIR],
            check=True,
        )
        logging.info("RAxML processing completed successfully.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running RAxML: {e}")
        sys.exit(1)


def generate_distance_table():

    """
    Summarizes genetic distances from RAxML output into a single CSV file.

    """
    logging.info("Generating summary table of genetic distances...")
    
    create_all_distances_table(RAXML_OUTPUT_DIR, DISTANCES_CSV_PATH)


def main():
    """
    Executes the full pipeline: tree pruning, RAxML distance calculation, and summary table creation.
    """
    logging.info("Starting Genetic Distance Project pipeline...")

    # Step 1: Run RAxML to calculate genetic distances
    # run_raxml_fixed_tree()

    # Step 2: Generate CSV summary of distances
    generate_distance_table()

    logging.info("All tasks completed successfully.")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception("An unexpected error occurred: %s", e)
        sys.exit(1)