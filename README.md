# Mobile resistance cassettes across *Acinetobacter baumannii* lineages

Data tables and figure-generation scripts for the comparative-genomics study of the
mobile-element context of acquired resistance genes across all public
*Acinetobacter baumannii* genomes.

This repository contains the summary data behind each figure and the Python scripts
that reproduce the figures from that data. It does **not** contain the raw genome
assemblies (retrieved from NCBI GenBank) or the per-genome intermediate outputs; it
contains the aggregated tables used for the figures.

---

## Repository layout

```
.
├── README.md
├── LICENSE
├── requirements.txt
├── data/        summary tables behind each figure
└── scripts/     one Python script per figure, each reads from ../data
```

## Quick start

```bash
# from the repository root
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt

cd scripts
python fig2_resistome.py      # writes fig2_resistome.png in the current directory
```

Each script is self-contained, reads its input from `../data/`, and writes a single
`.png` into the working directory. Run them from inside the `scripts/` directory.

---

## Data files (`data/`)

| File | Rows | Columns | Description |
|------|------|---------|-------------|
| `resistome_class.tsv` | 14 | `class`, `n_genomes`, `pct` | Genomes carrying at least one acquired gene of each drug class |
| `resistome_gene.tsv` | 20 | `gene`, `n_genomes`, `pct`, `class` | Twenty most prevalent individual acquired resistance genes |
| `linkage_top.tsv` | 25 | `gene`, `total`, `pct_plasmid`, `pct_chromosome`, `pct_integron`, `pct_IS`, `top_replicon`, `top_IS_family` | Mobile-element context of the 25 most prevalent acquired resistance genes |
| `upset_matrix.tsv` | 48,968 | `acc` + 15 binary gene columns | Per-genome presence/absence of the fifteen most prevalent acquired resistance genes |
| `repmeta.tsv` | 121 | `cluster`, `rep_acc`, `cluster_size` | PopPUNK lineage representatives and the number of genomes each represents |
| `tree.treefile` | — | Newick | Core-genome maximum-likelihood tree of the 121 lineage representatives |
| `yearly.tsv` | 37 | `year`, `n_genomes`, `mean_amr_load`, and per-gene prevalence columns | Yearly resistance-gene prevalence and mean acquired-gene load |
| `country.tsv` | 20 | `country`, `n_genomes`, `pct_carbapenemase`, `pct_blaOXA23`, `pct_armA`, `pct_sul2` | Resistance prevalence for the twenty most-sampled countries |

Percentages in `linkage_top.tsv` are fractions of gene-carrying contigs assigned to
each context; a gene may be scored in more than one mobile-element context on the
same contig, so the IS and integron fractions are not mutually exclusive.

## Scripts (`scripts/`)

| Script | Figure | Input | Output |
|--------|--------|-------|--------|
| `fig1_workflow.py` | Figure 1 — analysis workflow (schematic) | none | `fig1_workflow.png` |
| `fig2_resistome.py` | Figure 2 — acquired resistome overview | `resistome_class.tsv`, `resistome_gene.tsv` | `fig2_resistome.png` |
| `fig3_upset.py` | Figure 3 — gene co-carriage (UpSet) | `upset_matrix.tsv` | `fig3_upset.png` |
| `fig4_linkage.py` | Figure 4 — linkage architecture | `linkage_top.tsv` | `fig4_linkage.png` |
| `fig5_heatmap.py` | Figure 5 — genomic-context matrix | `linkage_top.tsv` | `fig5_heatmap.png` |
| `fig6_network.py` | Figure 6 — gene / mobile-element network | `linkage_top.tsv` | `fig6_network.png` |
| `fig7_phylo.py` | Figure 7 — circular ML phylogeny | `tree.treefile`, `repmeta.tsv` | `fig7_phylo.png` |
| `fig8_trends.py` | Figure 8 — temporal and geographic trends | `yearly.tsv`, `country.tsv` | `fig8_trends.png` |

The figure scripts regenerate the figures from the summary tables. The published
figures were produced with the same data and may differ in minor cosmetic details
(exact font sizes, label placement) from a fresh run.

---

## Analysis provenance

The summary tables were produced by a CPU-only, locally-executed pipeline on the
PARAM Rudra cluster. Genome assemblies were retrieved in bulk from NCBI GenBank and
re-annotated uniformly; no step depends on per-item external calls. The tools and
reference databases used, with the versions recorded in the manuscript Methods:

| Tool / resource | Version | Role |
|-----------------|---------|------|
| NCBI Datasets CLI | 18.37.0 | Bulk assembly and metadata retrieval |
| CheckM2 | 1.1.0 | Assembly completeness and contamination |
| Bakta | 1.11.4 | Uniform genome re-annotation |
| Bakta database (light) | v6.0 (2025-02-24) | Annotation reference |
| AMRFinderPlus | 4.2.7 | Acquired resistance-gene and point-mutation detection |
| AMRFinderPlus database | 2026-10-29.1 | Resistance reference |
| Panaroo | 1.8.0 | Pan-genome and core-gene alignment |
| PopPUNK | 2.7.6 | Lineage assignment |
| MOB-suite (mob_recon) | 3.1.9 | Plasmid reconstruction and replicon typing |
| IntegronFinder | 2.0.6 | Integron detection |
| ISEScan | 1.7.3 | Insertion-sequence detection |
| IQ-TREE | 3.1.4 | Maximum-likelihood phylogeny |
| Python | 3.11.16 | Scripting and aggregation |

---

## Citation

If you use this repository, please cite the archived release:

> Satyam, S. *Data and figure scripts: Mobile resistance cassettes across
> Acinetobacter baumannii lineages* [Data set]. Zenodo.
> https://doi.org/10.5281/zenodo.24133501

**Archived release (Zenodo):** https://doi.org/10.5281/zenodo.24133501
**Code repository:** https://github.com/krsatyamiitr/P3-Abaumannii-pan-resistome
**Article:** to be added on acceptance.

## License

Data files and scripts in this repository are released under the Creative Commons
Attribution 4.0 International License (CC BY 4.0); see [`LICENSE`](LICENSE). Reuse is
permitted, including commercially, provided appropriate credit is given.
