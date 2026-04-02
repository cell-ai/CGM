"""
This script receives as input a CSV file containing a list of genes (name and ensemble id). From this list, a search is made for orthologous (mammalian) sequences for each gene with the
Ensembl rest API. FASTA and CSV files are generated as output for each gene. CSV files contain information for each extracted sequence. The script also includes the homo sapiens sequence for each gene in the FASTA file with a new API query.

"""

import csv
import requests
import time
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
import concurrent.futures


def requests_retry_session(
    retries=3, backoff_factor=0.3, status_forcelist=(500, 502, 504), session=None
):
    session = session or requests.Session()
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    return session


# function to fetch mammalian orthologs sequences using Ensembl API REST
def fetch_orthologs(ensembl_id):
    server = "https://rest.ensembl.org"
    mammalian_target_taxon = "40674"
    type_of_homology = "orthologues"
    ext = f"/homology/id/{ensembl_id}?target_taxon={mammalian_target_taxon};type={type_of_homology}"
    headers = {"Content-Type": "application/json"}

    try:
        response = requests_retry_session().get(server + ext, headers=headers)

        if not response.ok:
            print("Request gone wrong!!!")
            response.raise_for_status()
            return []

        decoded = response.json()
        mammalian_orthologs = []

        data = decoded.get("data", [])
        if not data:
            print(f"No data found for {ensembl_id}")
            return mammalian_orthologs

        for homology in data[0].get("homologies", []):
            target = homology.get("target", {})
            mammalian_orthologs.append(
                {
                    "Gene_Ensembl_ID": ensembl_id,
                    "Ortholog_Species": target.get("species", ""),
                    "Ortholog_Protein_ID": target.get("protein_id"),
                    "Percentage_ID": target.get("perc_id"),
                    "Percentage_Pos": target.get("perc_pos"),
                }
            )

        return mammalian_orthologs

    except requests.exceptions.RequestException as e:
        print(f"Failed to connect to Ensembl API for orthologs: {e}")
        return []


# function to fetch sequence in FASTA format
def fetch_sequence_fasta(protein_id, species_name):
    server = "https://rest.ensembl.org"
    ext = f"/sequence/id/{protein_id}?type=protein;"
    headers = {"Content-type": "text/x-fasta"}

    try:
        response = requests_retry_session().get(server + ext, headers=headers)

        if not response.ok:
            print(f"failed to fetch sequence for {protein_id}")
            return None

        # Add the species name to the FASTA header
        fasta_text = response.text
        fasta_header, fasta_sequence = fasta_text.split("\n", 1)
        modified_fasta_header = f"{fasta_header}|{species_name}"
        modified_fasta_text = f"{modified_fasta_header}\n{fasta_sequence}"

        return modified_fasta_text

    except requests.exceptions.RequestException as e:
        print(f"Failed to connect to Ensembl API for sequence: {e}")
        return None


# function to fetch human protein sequence from a gene id
def fetch_human_protein_sequence(ensembl_id):
    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{ensembl_id}?expand=1"
    headers = {"Content-Type": "application/json"}

    try:
        response = requests_retry_session().get(server + ext, headers=headers)

        if not response.ok:
            print(f"Failed to fetch data for {ensembl_id}")
            return None

        decoded = response.json()
        if "Transcript" not in decoded:
            print(f"No transcripts found for {ensembl_id}")
            return None

        protein_id = None

        # Get the protein ID from the canonical transcript of homo sapiens
        for transcript in decoded["Transcript"]:
            if "Translation" in transcript and transcript.get("is_canonical"):
                protein_id = transcript["Translation"]["id"]
                break

        # se não tiver sequencia canônica para o gene, pegar o primeiro transcrito da lista
        if protein_id is None:
            for transcript in decoded["Transcript"]:
                if "Translation" in transcript:
                    protein_id = transcript["Translation"]["id"]
                    break

        if protein_id is None:
            print(f"No protein ID found for {ensembl_id}")
            return None

        # coletando a sequencia de proteina canonica
        return fetch_sequence_fasta(protein_id, "homo_sapiens")

    except requests.exceptions.RequestException as e:
        print(f"Failed to connect to Ensembl API for human protein: {ensembl_id}")
        return None


def import_orthologs_fasta(csvInputFilePath: str, outputFolderPath: str, max_workers=8): #obs definir aqui a qntdd de threads que deseja usar - cuidado para não atingir o limite da API 
    # Carrega a lista de genes do arquivo CSV
    gene_list = []
    with open(csvInputFilePath, "r") as readfile:
        reader = csv.DictReader(readfile, delimiter="\t")
        for row in reader:
            if row["Ensembl Canonical"] == "1" and row["Gene name"]:
                gene_list.append((row["Gene stable ID"], row["Gene name"]))

    # Processamento para cada gene na lista
    for ensembl_id, name in gene_list:
        print(f"Processando Ensembl ID: {ensembl_id}, Name: {name}")

        # Buscar a sequência da proteína humana
        human_protein_sequence = fetch_human_protein_sequence(ensembl_id)
        if human_protein_sequence is None:
            print(
                f"Não foi possível obter a sequência da proteína humana para {ensembl_id}. Pulando."
            )
            continue

        # Buscar ortólogos
        orthologs = fetch_orthologs(ensembl_id)
        if not orthologs:
            print(f"Nenhum ortólogo encontrado para {ensembl_id}. Pulando...")
            continue

        # Dicionário para armazenar as sequências dos ortólogos
        ortholog_sequences = {}

        # Busca as sequências em paralelo
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_id = {
                executor.submit(
                    fetch_sequence_fasta,
                    ortholog["Ortholog_Protein_ID"],
                    ortholog["Ortholog_Species"],
                ): ortholog
                for ortholog in orthologs
            }

            for future in concurrent.futures.as_completed(future_to_id):
                ortholog = future_to_id[future]
                try:
                    ortholog_sequence = future.result()
                    if ortholog_sequence:
                        ortholog_sequences[
                            ortholog["Ortholog_Protein_ID"]
                        ] = ortholog_sequence
                except Exception as exc:
                    print(f"Erro durante a busca da sequência: {exc}")

        # Escrevendo no arquivo CSV
        with open(f"{outputFolderPath}/{name}.csv", "w", newline="") as write_file:
            fieldnames = [
                "Gene_Ensembl_ID",
                "Ortholog_Species",
                "Ortholog_Protein_ID",
                "Percentage_ID",
                "Percentage_Pos",
            ]
            writer = csv.DictWriter(write_file, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(orthologs)

        # Escrevendo no arquivo FASTA
        with open(f"{outputFolderPath}/{name}.fasta", "w") as fasta_file:
            if human_protein_sequence:
                fasta_file.write(human_protein_sequence + "\n")
            for ortholog in orthologs:
                protein_id = ortholog["Ortholog_Protein_ID"]
                if protein_id in ortholog_sequences:
                    fasta_file.write(ortholog_sequences[protein_id] + "\n")
