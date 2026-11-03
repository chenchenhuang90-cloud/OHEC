# OHEC
Website logic script info
# OEC‑core: Standalone core calculation scripts for OEC web server
This repository contains **pure computational logic** used for therapeutic target screening in esophageal cancer for the associated manuscript.

Interactive web portal (OEC): **http://oncoheel.com/**

## Project overview
The OEC web server integrates data from 27 esophageal cancer cell lines, including gene mutation profiles and CRISPR‑knockout functional effects, to facilitate screening of esophageal cancer therapeutic targets.

This archive contains the full CRISPR knockout and mutation datasets of esophageal cancer TE cell lines used in the OEC manuscript for therapeutic target screening.
The standalone core calculation script (pure_gene_analysis.py) is available at the public GitHub repository: https://github.com/chenchenhuang90-cloud/OHEC
Interactive web portal: http://oncoheel.com/
All the data in the manuscript can be accessed at https://doi.org/10.5281/zenodo.23912088 and on the website

## Core functions in `core_filter.py`
1. `search_mutation_gene()`: Identify shared candidate genes across multiple esophageal cancer cell lines based on mutation data.
2. `search_crispr_gene()`: Extract gene‑effect records from CRISPR knockout datasets.
3. `compute_synthetic_lethality_pvalue()`: Evaluate candidate synthetic‑lethal gene combinations using two‑sided Mann‑Whitney U test (normal approximation, continuity correction and tie correction).
4. Internal helper functions for safe value parsing and statistical computation.

## Dependencies
Python >=3.9, **only Python standard library**. No third‑party packages required.

## Reproducibility
1. Clone or download this repository
2. Run: `python core_filter.py` to view the demo
3. Load your own input data (JSON format, same structure as demo) and call the functions
4. The output should match the screening result reported in the manuscript
python core_filter.py
