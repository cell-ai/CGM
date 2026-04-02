# This script filters FASTA files so that it leaves only one sequence for each species, choosing the one with the highest % identity.
# Receives FASTA and CSV files as input with information about %id for each sequence.

# O script realiza a filtragem de arquivos FASTA com sequências ortólogas, de modo que deixe apenas uma sequência para cada espécie,
# escolhendo a que tiver a maior % identidade. Recebe como input arquivos FASTA e CSV com informações sobre % id de cada sequência.

import csv
import os
from Bio import SeqIO


desired_species = [
    "homo_sapiens",
    "macaca_mulatta",
    "pan_troglodytes",
    "capra_hircus",
    "oryctolagus_cuniculus",
    "equus_caballus",
    "canis_lupus_familiaris",
    "cavia_porcellus",
    "felis_catus",
    "sus_scrofa",
    "ovis_aries_rambouillet",
    "mustela_putorius_furo",
    "mesocricetus_auratus",
    "rattus_norvegicus",
    "mus_musculus",
    "bos_taurus",
    "monodelphis_domestica",
]


def read_csv_and_filter(csv_path):
    species_dict = {}
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ortholog_species = row["Ortholog_Species"]
            percentage_id = float(row["Percentage_ID"])
            if ortholog_species in desired_species:
                if ortholog_species in species_dict:
                    if percentage_id > species_dict[ortholog_species]["Percentage_ID"]:
                        species_dict[ortholog_species] = {
                            "ID": row["Ortholog_Protein_ID"],
                            "Percentage_ID": percentage_id,
                        }
                else:
                    species_dict[ortholog_species] = {
                        "ID": row["Ortholog_Protein_ID"],
                        "Percentage_ID": percentage_id,
                    }

    return species_dict


def filter_fasta(fasta_path, species_dict, output_path):
    with open(output_path, "w") as out_f:
        for record in SeqIO.parse(fasta_path, "fasta"):
            description_parts = record.description.split("|")
            if len(description_parts) > 1:
                species = description_parts[1].split(" ")[0]
            else:
                print(
                    f"Skipping record with unexpected description format:{record.description}"
                )
                continue
            print(record.description)
            if "homo_sapiens" in record.description:
                SeqIO.write(record, out_f, "fasta")
            elif species in species_dict:
                if species_dict[species]["ID"] in record.description:
                    SeqIO.write(record, out_f, "fasta")


def entry_point(input_folder_with_csv_and_fasta, output_folder):
    for filename in os.listdir(input_folder_with_csv_and_fasta):
        print(filename)
        if filename.endswith(".csv"):
            base_name = os.path.splitext(filename)[0]
            corresponding_fasta = f"{input_folder_with_csv_and_fasta}/{base_name}.fasta"
            output_fasta = f"{output_folder}/filtered_{base_name}.fasta"

            if os.path.exists(corresponding_fasta):
                species_dict = read_csv_and_filter(
                    f"{input_folder_with_csv_and_fasta}/{filename}"
                )
                filter_fasta(corresponding_fasta, species_dict, output_fasta)
