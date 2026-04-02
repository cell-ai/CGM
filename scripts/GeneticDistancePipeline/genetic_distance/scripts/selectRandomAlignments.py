import shutil
import os
import random 

folder_distances = "/home/maiara/Documents/02-2025-mestrado/data_processing/fastas/aligned_fasta/outputAlignmentDistances_allSpecies"
folder_random = "/home/maiara/Documents/02-2025-mestrado/data_processing/fastas/aligned_fasta/teste"

os.makedirs(folder_random, exist_ok=True)

arquivos = [f for f in os.listdir(folder_distances) if f.endswith("-aligned-trim.fasta")]

random = random.sample(arquivos, 200)

for file in random:
    caminho_origem = os.path.join(folder_distances, file)
    caminho_destino = os.path.join(folder_random, file)
    shutil.copy2(caminho_origem, caminho_destino)

print("Arquivos copiados com sucesso para a pasta de teste!")