"""
This script creates a single table that summarizes the genetic distances between pairs of species for each gene analyzed. It receives RAxML output files as input and generates a csv table.

Esse script cria uma única tabela que resume as distâncias genéticas entre pares de espécies para cada gene analisado. Recebe como input arquivos de saída do RAxML e gera uma tabela csv. 

"""

from collections import defaultdict
import csv
import os
import locale


def read_raxml_distances_matrix(file_path, gene_name):
    distances = defaultdict(dict)
    try:
        with open(file_path, "r") as f:
            for line in f:
                seq_ids, distance = line.strip().split("\t")
                species1 = seq_ids.split(" ")[0].split("|")[1]
                species2 = seq_ids.split(" ")[1].split("|")[1]
                species_pair = (species1, species2)
                distances[species_pair][gene_name] = float(distance)
    except Exception as e:
        return str(e)
    return distances


def update_csv_table(csv_file_path, distances, gene_name):
    existing_data = defaultdict(dict)
    try:
        with open(csv_file_path, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                species_pair = (row["Species1"], row["Species2"])
                existing_data[species_pair] = row
    except FileNotFoundError:
        pass
    except Exception as e:
        return False, str(e)

    for species_pair, gene_distances in distances.items():
        species1, species2 = species_pair
        if species_pair not in existing_data:
            existing_data[species_pair] = {"Species1": species1, "Species2": species2}
        existing_data[species_pair][gene_name] = gene_distances.get(gene_name, "N/A")

    # Update the fieldnames to include the new gene if it's not already included
    existing_fieldnames = (
        ["Species1", "Species2"]
        if not existing_data
        else list(next(iter(existing_data.values())).keys())
    )
    if gene_name not in existing_fieldnames:
        existing_fieldnames.append(gene_name)

    try:
        locale.setlocale(locale.LC_NUMERIC, "C")

        with open(csv_file_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=existing_fieldnames)
            writer.writeheader()
            for row in existing_data.values():
                writer.writerow(row)
    except Exception as e:
        return False, str(e)
    finally:
        locale.resetlocale(locale.LC_NUMERIC)
    return True, "Table successfully updated"


# directory_path = "/home/mra/Downloads/Analises-masters/Pipeline/teste/"

# csv_file_path = '/home/mra/Downloads/Analises-masters/Pipeline/teste/resultsDistances.csv'


def create_all_distances_table(directory_path, csv_file_path):
    for filename in os.listdir(directory_path):
        if filename.startswith("RAxML_distances.") and filename.endswith("_matrix"):
            raxml_matrix_file_path = os.path.join(directory_path, filename)

            # Automatically determine the gene name from the file name
            gene_name = filename.replace("RAxML_distances.", "").split("_matrix")[0]

            # Read the genetic distances from the RAxML_distances matrix file
            distances = read_raxml_distances_matrix(raxml_matrix_file_path, gene_name)

            # Update the existing CSV table with the new genetic distances
            update_status, message = update_csv_table(
                csv_file_path, distances, gene_name
            )

            # Print the update status and message
            print(f"Processing {filename}: {update_status}, {message}")
