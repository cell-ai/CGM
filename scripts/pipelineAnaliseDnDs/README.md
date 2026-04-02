# Dn/Ds Analysis Pipeline for Positive Selection Detection

## Overview

This pipeline automates the analysis of non-synonymous (Dn) versus synonymous (Ds) substitution rates in gene sequences through phylogenetic comparison. The goal is to identify genes under positive selection using the Maximum Likelihood (ML) method with PAML (Phylogenetic Analysis by Maximum Likelihood).

## Scientific Objective

The Dn/Ds ratio (omega) provides information about the type of natural selection acting on a gene:
- **Dn/Ds < 1**: Purifying selection (pressure against changes, function conservation)
- **Dn/Ds = 1**: Neutral evolution (no selective pressure)
- **Dn/Ds > 1**: Positive selection (favors changes, adaptation)

The pipeline calculates these ratios for selected genes using codon models with statistical tests (LRT and BEB).

## Architecture and Data Flow

The pipeline executes 11 main steps in sequence:

```
1. FETCH_CDS
   Retrieves CDS sequences from Ensembl
   Output: data/raw_cds/{gene}.fasta

2. TRANSLATE_CDS
   Translates CDS to protein
   Output: data/raw_prot/{gene}.fasta

3. FILTER_VALID_CDS
   Keeps only valid CDS (complete translation)
   Output: data/filtered_cds/{gene}.fasta

4. ALIGN_PROTEIN
   Aligns proteins with PRANK
   Output: data/align_prot/{gene}.fasta

5. PAL2NAL (back-translation)
   Converts protein alignment back to codons
   Output: data/align_codon/{gene}.fasta

6. REORDER_CODON_ALIGNMENT
   Reorders and validates codon alignment
   Output: data/align_codon_reordered/{gene}.fasta

7. BUILD_TREE
   Builds phylogenetic tree with IQ-TREE
   Output: data/trees/{gene}.treefile

8. ADD_HEADER_TREE
   Formats tree for PAML
   Output: data/trees_hdr/{gene}.tree

9. FASTA_TO_PHYLIP
   Converts alignment to PHYLIP format
   Output: data/align_codon_phylip/{gene}.phy

10. CODEML (PAML Analysis)
    Runs ML analysis with multiple models
    Output: results/paml/{gene}/mlc

11. ANALYZE_PAML_RESULTS
    Extracts statistical results (LRT and BEB)
    Output: results/paml/{gene}/lrt.tsv
            results/paml/{gene}/beb.tsv
```

## Current Configuration

Settings are defined in `config.yaml`:

```yaml
threads_mafft: 4          # Threads for parallel processing in PRANK
threads_total: 1          # Threads for IQ-TREE (set according to your resources)
genes_file: "data/genes/divergentes_e_interferon.tsv"  # Gene list to analyze
species_file: "data/species/mammals.txt"               # Species list (format: taxonomy IDs)
omega_cutoff: 0.6         # Threshold for filtering significant results
```

### Input File: `genes_file`

The file must be a TSV with two columns:
```
Gene name       Gene stable ID
IFNB1           ENSG00000171092
CD80            ENSG00000005728
CD86            ENSG00000114248
```

Each gene is analyzed individually. IDs must be Ensembl Gene IDs (ENSG*).

### Input File: `species_file`

Contains mammalian species taxonomy IDs, one per line:
```
9606    # Homo sapiens
10090   # Mus musculus
9913    # Bos taurus
...
```

## How to Run

### Prerequisites

1. Install Conda/Mamba
2. Create the environment:
```bash
conda env create -f envs/dnds.yaml
conda activate dnds
```

### Basic Execution

Run the full pipeline:
```bash
snakemake -c 4
```

Where `-c 4` sets the number of available cores/threads.

## Changing Configuration

### 1. Change the Gene File

Edit `config.yaml`:
```yaml
genes_file: "data/genes/your_file.tsv"
```

Your file must keep the same structure (gene name and ENSG ID).

### 2. Adjust Omega Threshold

Change the cutoff for result filtering:
```yaml
omega_cutoff: 1.0     # Accept only Dn/Ds > 1.0 (strong selection)
```

### 3. Add/Remove PAML Models

Edit the `ctl_global` rule in the Snakefile:
```coffeescript
NSsites = 0 1 2 7 8   # Models: M0, M1a, M2a, M7, M8
```

Model options:
- **0**: M0 (one-ratio) - one omega for the whole tree
- **1**: M1a (neutral) - neutral selection
- **2**: M2a (positive selection) - positive selection
- **7**: M7 (beta) - beta distribution
- **8**: M8 (beta&w) - beta plus positive selection

## Outputs and Results

### Output Structure

```
results/paml/
├── IFNB1/
│   ├── codeml.ctl          # PAML control file
│   ├── mlc                 # Full PAML output
│   ├── lrt.tsv             # LRT test (likelihood ratio test)
│   └── beb.tsv             # BEB test (Bayes Empirical Bayes)
├── CD80/
│   └── ...
└── CD86/
    └── ...
```

### Interpreting Results

#### LRT (Likelihood Ratio Test)

`lrt.tsv` file:
```
gene        model1  model2      chi2    pvalue  significant
IFNB1       M1a     M2a         15.34   0.0005  True
CD80        M1a     M2a         2.14    0.3421  False
```

- **chi2**: test statistic value
- **pvalue**: significance (< 0.05 = significant)
- **significant**: whether it passed the threshold

#### BEB (Bayes Empirical Bayes)

`beb.tsv` file:
```
gene    codon   aa  omega   posterior_prob
IFNB1   42      K   2.34    0.95
IFNB1   57      L   3.12    0.89
CD80    103     D   0.45    0.78
```

- **omega**: Dn/Ds at that specific codon
- **posterior_prob**: prediction confidence (ideally > 0.90)

## Directory Structure

```
analise_dn_ds/
├── config.yaml                    # Main configuration file
├── Snakefile                      # Snakemake workflow
├── envs/
│   └── dnds.yaml                  # Dependency packages
├── src/
│   ├── 01_fetch_cds.py            # Fetch CDS from Ensembl
│   ├── 02_translate_cds.py        # Translate to protein
│   ├── 03_reorder_codon_alignment.py  # Reorder alignment
│   ├── 04_add_header_tree.py      # Format tree
│   ├── 05_fasta_to_phylip.py      # Convert format
│   ├── 06_analyze_paml_output.py  # Process results
│   └── utils/
│       └── load_gene_file.py      # Helper function
├── data/
│   ├── genes/                     # Gene lists
│   ├── species/                   # Species lists
│   ├── raw_cds/                   # Downloaded CDS
│   ├── raw_prot/                  # Translated proteins
│   ├── filtered_cds/              # Filtered CDS
│   ├── align_prot/                # Protein alignments
│   ├── align_codon/               # Codon alignments
│   ├── align_codon_reordered/     # Reordered codon alignments
│   ├── align_codon_phylip/        # PHYLIP alignments
│   ├── trees/                     # IQ-TREE trees
│   └── trees_hdr/                 # PAML-formatted trees
└── results/
    └── paml/                      # PAML results per gene
```

## Tools Used

- **Snakemake**: Workflow orchestration
- **Biopython**: Sequence manipulation
- **PRANK**: Protein alignment (codon-aware)
- **MAFFT**: Alternative alignment (optional)
- **PAL2NAL**: Back-translation of alignments
- **IQ-TREE**: Phylogenetic tree construction
- **PAML/CodeML**: ML positive selection analysis

## References

- PAML documentation: http://abacus.gene.ucl.ac.uk/software/paml.html
- IQ-TREE: http://www.iqtree.org/
- Ensembl API: http://www.ensembl.org/
- Snakemake: https://snakemake.readthedocs.io/
