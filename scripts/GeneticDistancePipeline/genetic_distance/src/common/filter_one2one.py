import csv
import requests
import time
from Bio import SeqIO
import concurrent.futures


# Função para obter IDs de proteínas ortólogas one2one
def get_one2one_orthologs(gene_id, tentativas=3, intervalo_inicial=1):
    url = f"https://rest.ensembl.org/homology/id/{gene_id}?type=orthologues"
    headers = {"Content-Type": "application/json"}
    intervalo = intervalo_inicial

    for tentativa in range(tentativas):
        response = requests.get(url, headers=headers)
        if response.ok:
            data = response.json()
            one2one_protein_ids = [homology["target"]["protein_id"] for homology in data["data"][0]["homologies"] if homology["type"] == "ortholog_one2one"]
            return one2one_protein_ids
        else:
            print(f"Tentativa {tentativa + 1} falhou para o gene {gene_id}: {response.status_code}")
            time.sleep(intervalo)
            intervalo *= 2  # Aumenta o intervalo para o próximo backoff

    print(f"Não foi possível obter ortólogos após {tentativas} tentativas para o gene {gene_id}")
    return []


# Função para filtrar arquivo FASTA
def filter_fasta_by_orthologs(fasta_path, one2one_protein_ids, gene_name, output_dir, log_file):
    try:
        sequences = SeqIO.parse(fasta_path, "fasta")
    except FileNotFoundError:
        log_file.write(f"Arquivo FASTA não encontrado para {gene_name}\n")
        print(f"Arquivo FASTA não encontrado para {gene_name}\n")
        return
    except Exception as e:
        log_file.write(f"Erro ao ler o arquivo FASTA para {gene_name}: {e}\n")
        print(f"Erro ao ler o arquivo FASTA para {gene_name}: {e}\n")
        return

    filtered_sequences = []
    homo_sapiens_sequence = None

    for seq in sequences:
        if 'homo_sapiens' in seq.description.lower():
            homo_sapiens_sequence = seq
        else:
            seq_id_parts = seq.id.split(".")[0].split("|")[0]
            if seq_id_parts in one2one_protein_ids:
                filtered_sequences.append(seq)

    if filtered_sequences and homo_sapiens_sequence:
        filtered_sequences.append(homo_sapiens_sequence)
        output_path = f"{output_dir}/{gene_name}_filtered.fasta"
        SeqIO.write(filtered_sequences, output_path, "fasta")
        print(f"Arquivo FASTA filtrado criado para {gene_name}")
    elif not filtered_sequences:
        print(f"Nenhuma correspondência ortóloga one2one encontrada para {gene_name}")

def process_genes_from_tsv_paralelo(tsv_file_path, fasta_input_dir, fasta_output_dir, log_file_path, max_workers=5):
    with open(tsv_file_path, mode="r", newline="") as file, open(log_file_path, mode="w", newline="") as log_file:
        reader = csv.reader(file, delimiter="\t")
        gene_list = [row for row in reader if len(row) >= 2]

        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_gene = {
                executor.submit(process_gene, gene_name, ensembl_id, fasta_input_dir, fasta_output_dir, log_file): (gene_name, ensembl_id)
                for gene_name, ensembl_id in gene_list
            }

            for future in concurrent.futures.as_completed(future_to_gene):
                gene_name, ensembl_id = future_to_gene[future]
                try:
                    future.result()  # Espera a execução e processa o resultado
                except Exception as exc:
                    log_file.write(f"Erro ao processar {gene_name} ({ensembl_id}): {exc}\n")

def process_gene(gene_name, ensembl_id, fasta_input_dir, fasta_output_dir, log_file):
    print(f"Lendo {gene_name} ({ensembl_id}) do TSV")
    one2one_protein_ids = get_one2one_orthologs(ensembl_id)
    if not one2one_protein_ids:
        log_file.write(f"Nenhum ortólogo one2one encontrado para {gene_name} ({ensembl_id})\n")
        return

    fasta_path = f"{fasta_input_dir}/{gene_name}.fasta"
    print(f"Processando arquivo FASTA: {fasta_path}")
    filter_fasta_by_orthologs(fasta_path, one2one_protein_ids, gene_name, fasta_output_dir, log_file)

# Caminhos e chamada à função principal
tsv_file_path = "/home/maiara/pipeline/outputFasta/resto.tsv"
fasta_input_dir = "/home/maiara/pipeline/outputFasta/all"
fasta_output_dir = "/home/maiara/pipeline/outputFasta/orthologs_one2one"
log_file_path = "/home/maiara/pipeline/outputFasta/processamento2_orthologs_one2one.log"
process_genes_from_tsv_paralelo(tsv_file_path, fasta_input_dir, fasta_output_dir, log_file_path)
