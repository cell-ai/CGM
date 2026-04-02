# ComparativeGenomicsMammals

This repository brings together the main results and files from the comparative mammalian genomics project, with a focus on:
- gene-by-gene genetic distance relative to *Homo sapiens*
- functional enrichment (GO, KEGG, Reactome)
- positive selection analysis (PAML)
- species-level GSEA

## Where to start

The recommended reading order is:

1. `Results_details.pdf`
2. `Methods_details.pdf`
3. `Supplementary_Table_S1.xlsx` to `Supplementary_Table_S4.xlsx`

This way, you first see the main results, then the detailed methods, and finally the supplementary tables with supporting data.

## What is in the root folder

- `Results_details.pdf`: main document with project results.
- `Methods_details.pdf`: detailed description of methods, pipelines, and analysis criteria.
- `Supplementary_Table_S1.xlsx` to `Supplementary_Table_S4.xlsx`: supplementary tables.
- `data/`: data outputs (raw, statistical summaries, and analysis outputs).
- `scripts/`: main analysis scripts and pipelines.

## What each supplementary table contains

- `Supplementary_Table_S1.xlsx`
  - Summary table of the PAML analysis, showing the number of sites found under positive selection for each gene.

- `Supplementary_Table_S2.xlsx`
  - Site-level selective pressure analysis for the 3 genes with evidence of the highest number of sites under positive selection.
  - The table allows identification of amino acid sites with evidence of positive selection.

- `Supplementary_Table_S3.xlsx`
  - List of genes present in the top 5% most divergent set across all 14 target species in the analysis.
  - The table shows gene lists per species and genes appearing at the distribution tail in a larger number of species.

- `Supplementary_Table_S4.xlsx`
  - Summary table with GSEA results performed for each species using genetic distances to *Homo sapiens*.
  - The table highlights pathways with the highest mean NES and significant enrichment in most species.

## `data/` structure

### `data/GeneticDistances (raw and stats)/`
Base datasets and descriptive statistics for genetic distance.

Main files:
- `allSpecies_resultsDistances.csv`: per-gene/per-species genetic distances calculated using maximum likelihood phylogenetic method.

- `medianas_genes.csv`: median distance per gene.
- `medianas_especies.csv`: median distance per species and number of evaluated genes.
- `genes_p15.csv`: genes at the lower tail of the distribution (more conserved).
- `genes_p85.csv`: genes at the upper tail of the distribution (more divergent).
- `stats_genes.csv`: per-gene summary statistics.

### `data/GeneticDistances (analysis)/`
Functional analysis outputs derived from genetic distances.

Main subfolders:
- `enrichment_results/`: GO (BP/CC/MF) and KEGG enrichment outputs for selected gene sets.
- `PathwayAnalysis/`: Reactome and GSEA outputs, including:
  - consolidated files (e.g., `reactome_pathway_significants.csv`, `nes_all_species.csv`, `leading_edge_genes_all_species.csv`)
  - species-specific subfolders with dedicated outputs
  - figures folder (`Figures/`)

### `data/PositiveSelection (PAML results)/`
Final outputs from PAML-based positive selection analyses.

Main files:
- `final_lrt_corrected.tsv`: LRT test results with multiple testing correction (FDR).
- `final_beb_summary.tsv`: summary of sites with evidence of positive selection using the BEB approach.

## `scripts/` structure

> Important: this folder contains the **main analysis scripts** used to generate the core results. It may not include every auxiliary or exploratory script from the project.

### `scripts/R/`
Main statistical and functional analysis scripts:
- `01_prepare_distances.R`: prepares and filters the distance matrix.
- `02_gene_species_stats.R`: computes gene/species statistics and defines cutoffs (p15/p85).
- `03_enrichment_go_kegg.R`: GO and KEGG enrichment for selected genes.
- `04_pathway_medians_reactome.R`: Reactome pathway analysis with pathway-level statistics.
- `05_gsea_species_reactome.R`: species-level GSEA using Reactome pathways.

### `scripts/GeneticDistancePipeline/`
Main pipeline for ortholog processing, alignments, and genetic distance calculation. 

### `scripts/pipelineAnaliseDnDs/`
Snakemake pipeline for positive selection (dN/dS) analysis with PAML.
Includes configuration (`config.yaml`), workflow rules (`Snakefile`), and supporting code in `src/`.
