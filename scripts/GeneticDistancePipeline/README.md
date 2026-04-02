# GeneticDistancePipeline

Pipeline to retrieve mammalian orthologs, align protein sequences, and calculate pairwise genetic distances using maximum likelihood.

## Overview

This directory contains a set of Python, Bash, and R scripts used to:

1. Read a list of human genes.
2. Retrieve orthologous sequences from the Ensembl REST API.
3. Filter sequences to keep one sequence per species.
4. Align sequences with MAFFT.
5. Remove poorly aligned regions with trimAl.
6. Calculate genetic distances with RAxML.
7. Consolidate distances into a final per-gene table.

The code is organized into two main workflows:

- `variable_tree`: calculates distances directly from alignments without providing a fixed tree to RAxML.
- `fixed_tree`: uses a pruned reference tree for each gene before calculating distances.

## Directory structure

```text
GeneticDistancePipeline/
├── README.md
└── genetic_distance/
    ├── genes_biomart/
    │   └── mart_export.txt
    ├── scripts/
    │   ├── check_quality_prune_tree.R
    │   ├── pathwayAnalysisSpecies.R
    │   ├── prune_trees.R
    │   ├── reference_tree_from_consensus.R
    │   └── selectRandomAlignments.py
    └── src/
        ├── main.py
        ├── main_fixed_tree.py
        ├── main_variable_tree.py
        ├── common/
        │   ├── filter_fasta.py
        │   ├── filter_one2one.py
        │   └── importOrthologsFasta.py
        ├── fixed_tree/
        │   ├── clean_fasta.sh
        │   ├── ML_distance_fixed_tree.sh
        │   └── table_distances2.py
        └── variable_tree/
            ├── alignDistance.sh
            └── table_distances.py
```

## Pipeline workflow

### 1. Ortholog retrieval

The script `importOrthologsFasta.py` receives a gene table and, for each canonical human gene, queries Ensembl to:

- retrieve mammalian orthologs;
- obtain the protein sequence for each ortholog;
- retrieve the canonical human protein sequence for the same gene.

As output, it generates for each gene:

- a `.csv` file with ortholog metadata;
- a `.fasta` file with the human sequence plus ortholog sequences.

### 2. Sequence filtering

The script `filter_fasta.py` reads each gene `.csv` and `.fasta` and keeps only:

- `homo_sapiens`;
- one sequence per species among a predefined set of mammals;
- the sequence with the highest `Percentage_ID` when duplicates exist for the same species.

### 3. Alignment and trimming

In the `variable_tree` workflow, `alignDistance.sh`:

- runs multiple sequence alignment with `mafft --auto`;
- removes poorly aligned regions with `trimal -automated1`;
- runs `raxmlHPC` to estimate pairwise distances.

### 4. Fixed tree per gene

In the `fixed_tree` workflow, a reference tree is pruned to keep only the species present in each alignment. RAxML then uses this tree as a fixed topology to estimate distances.

### 5. Distance consolidation

The scripts `table_distances.py` and `table_distances2.py` iterate through `RAxML_distances.*_matrix` files and build a final CSV where:

- each row represents a species pair;
- each additional column represents a gene;
- each cell contains the estimated genetic distance for that gene.

## Expected inputs

### Gene list

An example input file is available at `genetic_distance/genes_biomart/mart_export.txt` and uses tab-delimited columns such as:

- `Gene name`
- `Ensembl Canonical`
- `Gene stable ID`

The script only considers rows where:

- `Ensembl Canonical == 1`
- `Gene name` is not empty

### External dependencies

Besides Python and R packages, the scripts depend on command-line tools installed on the system:

- `mafft`
- `trimal`
- `raxmlHPC`

## Language-specific dependencies

### Python

The scripts mainly use:

- `requests`
- `biopython`
- standard modules such as `csv`, `os`, `subprocess`, `logging`, `concurrent.futures`

### R

R scripts use, depending on the step:

- `ape`
- `seqinr`
- `dplyr`
- `tidyr`
- `stringr`
- `data.table`
- `clusterProfiler`
- `fgsea`
- `ggplot2`

## Entry points

### `genetic_distance/src/main.py`

Main pipeline entry point for the variable-tree workflow. It defines functions to:

- retrieve orthologs;
- filter FASTA files;
- align sequences and calculate distances;
- consolidate the final distance table.

### `genetic_distance/src/main_fixed_tree.py`

Entry point for the fixed-tree workflow. It is responsible for:

- calling the Bash script that runs RAxML with gene-specific pruned trees;
- generating the consolidated distance table.


## File descriptions

### `genetic_distance/src/common/importOrthologsFasta.py`

Retrieves mammalian orthologs through the Ensembl REST API and builds per-gene `.csv` and `.fasta` files.

Main functions:

- `fetch_orthologs()`: retrieves orthologs for a gene.
- `fetch_sequence_fasta()`: downloads the protein sequence for a `protein_id`.
- `fetch_human_protein_sequence()`: retrieves the canonical human protein sequence.
- `import_orthologs_fasta()`: orchestrates the full process and writes output files.

### `genetic_distance/src/common/filter_fasta.py`

Filters FASTA files to keep only one sequence per target species, selecting the one with highest `% identity` according to the corresponding CSV.

Main functions:

- `read_csv_and_filter()`: reads CSV and chooses the best ortholog per species.
- `filter_fasta()`: rewrites FASTA with filtered sequences.
- `entry_point()`: iterates over a directory and processes all CSV/FASTA pairs.

### `genetic_distance/src/common/filter_one2one.py`

Auxiliary script to filter only `ortholog_one2one` entries. It queries Ensembl, identifies valid `protein_id` values, and generates filtered FASTA files.

### `genetic_distance/src/variable_tree/alignDistance.sh`

Bash script for the `variable_tree` workflow. For each input FASTA, it:

- runs MAFFT;
- runs trimAl;
- runs RAxML with `-f x` for pairwise distances.

Error logs are written to `error_log.txt` in the output directory.

### `genetic_distance/src/variable_tree/table_distances.py`

Reads `RAxML_distances.*_matrix` files from the `variable_tree` workflow and updates a master CSV with pairwise species distances for each gene.

### `genetic_distance/src/fixed_tree/ML_distance_fixed_tree.sh`

Bash script for the `fixed_tree` workflow. For each pruned tree `*_pruned.nwk`, it:

- locates the corresponding alignment;
- cleans sequence names with `clean_fasta.sh`;
- runs `raxmlHPC` using the provided tree with `-t`;
- writes a `RAxML_distances.<GENE>_matrix` file.

### `genetic_distance/src/fixed_tree/clean_fasta.sh`

Normalizes FASTA headers to match species names in tree files:

- removes prefixes before species names;
- fixes taxonomic synonyms;
- reduces names to `genus_species`;
- removes `canis_lupus_dingo` to avoid duplication.

### `genetic_distance/src/fixed_tree/table_distances2.py`

Equivalent to `table_distances.py`, but adapted to the output format used in the fixed-tree workflow.

### `genetic_distance/scripts/reference_tree_from_consensus.R`

Builds a reference tree from a larger consensus tree, keeping only species present in the distance analysis. It also normalizes species names and fixes synonyms before saving the reduced tree.

### `genetic_distance/scripts/prune_trees.R`

Reads the reference tree and, for each alignment, generates a pruned tree containing only species actually present for that gene.

### `genetic_distance/scripts/check_quality_prune_tree.R`

Manual inspection script to validate whether sequence names in alignments exactly match labels in the pruned tree.

### `genetic_distance/scripts/selectRandomAlignments.py`

Randomly selects a subset of already processed alignments and copies them to a test directory. Useful for manual validation and debugging.

### `genetic_distance/genes_biomart/mart_export.txt`

Example file with the gene list used as pipeline input.

## Outputs generated across the process

Depending on the workflow executed, the pipeline produces:

- `.csv` files with per-gene ortholog metadata;
- `.fasta` files with per-gene ortholog sequences;
- `*-aligned.fasta` alignments;
- `*-aligned-trim.fasta` filtered alignments;
- `*_pruned.nwk` pruned trees;
- `RAxML_distances.*_matrix` matrices;
- a final consolidated distance table by gene.

### Example run for the variable-tree workflow

```bash
cd genetic_distance/src
python main_variable_tree.py
```

### Example run for the fixed-tree workflow

```bash
cd genetic_distance/src
python main_fixed_tree.py
```

## Summary

- `genetic_distance/src/main_variable_tree.py`: main workflow without fixed tree;
- `genetic_distance/src/main_fixed_tree.py`: main workflow with fixed tree;
- `genetic_distance/src/common/importOrthologsFasta.py`: retrieves orthologs;
- `genetic_distance/src/common/filter_fasta.py`: filters sequences by species;
- `genetic_distance/src/variable_tree/alignDistance.sh`: aligns and calculates distances;
- `genetic_distance/src/fixed_tree/ML_distance_fixed_tree.sh`: calculates distances with fixed tree;
- `genetic_distance/src/variable_tree/table_distances.py` and `genetic_distance/src/fixed_tree/table_distances2.py`: build the final table.
