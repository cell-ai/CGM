#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
main.py: Main pipeline for bioinformatics project 

This pipeline consists of 4 scripts and has the main function of searching for mammalian orthologous sequences for each gene and creating a table that summarizes the genetic distances between each species for each gene in the list.
It receives as input a csv file containing a list of genes (name and ensembl id) and generates as main output a csv file with the genetic distances between pairs of sequences for each gene.


Esta pipeline é composta por 4 scripts e tem como função principal buscar sequências ortólogas de mamíferos para cada gene e criar uma tabela que resume as distâncias genéticas entre as espécies para cada gene da lista.
Recebe como entrada um arquivo csv contendo uma lista de genes (nome e ensembl id) e gera como saída principal um arquivo csv com as distâncias genéticas entre pares de sequências para cada gene.
"""

# Import standard libraries
import sys
import logging
import subprocess
import os

# Import local modules
from common.importOrthologsFasta import import_orthologs_fasta
from common.filter_fasta import entry_point
from variable_tree.table_distances import create_all_distances_table

# Initialize logging
logging.basicConfig(level=logging.INFO)


def fetchOrthologs(CSV_INPUT_FILE_PATH, OUTPUT_FOLDER_PATH_1):
    """
    Starting point of the script: fetch the orthologous sequences from a gene list.
    """
    logging.info("Starting the first part of pipeline...")
    import_orthologs_fasta(CSV_INPUT_FILE_PATH, OUTPUT_FOLDER_PATH_1)


def removeDuplicates(INPUT_CSV_PATH, OUTPUT_FOLDER_PATH_2):
    """
    Remove duplicated sequences by selecting those with higher %id.
    """
    logging.info("Filtering FASTA files...")
    entry_point(INPUT_CSV_PATH, OUTPUT_FOLDER_PATH_2)


def alignDistances(OUTPUT_FOLDER_PATH_2, OUTPUT_FOLDER_PATH_3):
    """
    Multiple Sequence Alignment and other tasks.
    """
    logging.info("Aligning the sequences and calculating pairwise genetic distances...")
    subprocess.run(
        ["bash", "src/alignDistance.sh", OUTPUT_FOLDER_PATH_2, OUTPUT_FOLDER_PATH_3]
    )


def createDataframeDistances(OUTPUT_FOLDER_PATH_3, CSV_PATH):
    """
    Summarize genetic distances into a dataframe.
    """
    logging.info(f"Updating results table in {OUTPUT_FOLDER_PATH_3} directory...")
    create_all_distances_table(OUTPUT_FOLDER_PATH_3, CSV_PATH)


def main():
    # Paths for input and output files and directories
    CSV_INPUT_FILE_PATH = (
        "/home/maiara/pipeline/Genes_biomart/mart_export.txt"
    )
    OUTPUT_FOLDER_PATH_1 = "/home/maiara/pipeline/outputFasta/orthologs_selectedSpecies"  # remember to change the name of the folders for the gene name!!
    OUTPUT_FOLDER_PATH_2 = (
        "/home/maiara/pipeline/outputFasta/mouse_strains/mouse_filtered"
    )
    OUTPUT_FOLDER_PATH_3 = (
        "/home/maiara/pipeline/outputAlignmentDistancesStrains"
    )
    OUTPUT_FOLDER_PATH_4 = (
        "/home/maiara/pipeline/distances"
    )
    CSV_RESULT_DISTANCES_PATH = OUTPUT_FOLDER_PATH_4 + "/strains/strains_resultsDistances.csv"

    # Create folders if they don't exist
    os.makedirs(OUTPUT_FOLDER_PATH_1, exist_ok=True)
    os.makedirs(OUTPUT_FOLDER_PATH_2, exist_ok=True)
    os.makedirs(OUTPUT_FOLDER_PATH_3, exist_ok=True)
    os.makedirs(OUTPUT_FOLDER_PATH_4, exist_ok=True)

    # Run each step of the pipeline
    #fetchOrthologs(CSV_INPUT_FILE_PATH, OUTPUT_FOLDER_PATH_1)
    #removeDuplicates(OUTPUT_FOLDER_PATH_1, OUTPUT_FOLDER_PATH_2)
    #alignDistances(OUTPUT_FOLDER_PATH_2, OUTPUT_FOLDER_PATH_3)
    create_all_distances_table(OUTPUT_FOLDER_PATH_3, CSV_RESULT_DISTANCES_PATH)

    logging.info("All tasks completed successfully.")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception("An error occurred: %s", e)
        sys.exit(1)
