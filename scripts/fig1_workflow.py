#!/usr/bin/env python3
"""Figure 1 - Overview of the analysis workflow (schematic).
Public genomes retrieved and quality-filtered, re-annotated uniformly, then four
parallel analysis tracks feeding a contig-level linkage integration step.
This figure is a schematic and does not read from the data tables.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(13, 12)); ax.set_xlim(0, 10); ax.set_ylim(0, 13)
ax.set_axis_off()

def box(x, y, w, h, text, fc, fontsize=12):
    p = FancyBboxPatch((x-w/2, y-h/2), w, h, boxstyle="round,pad=0.08",
                       fc=fc, ec="#555555", lw=1.3)
    ax.add_patch(p)
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize, zorder=5)

def arrow(x0, y0, x1, y1, label=None):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                 mutation_scale=18, lw=1.6, color="#333333"))
    if label:
        ax.text((x0+x1)/2 + 0.35, (y0+y1)/2, label, fontsize=10,
                style="italic", ha="left", va="center")

BLUE="#cfe0f3"; GREY="#e3e3e3"; PINK="#f6d9d6"; GREEN="#d6ecd9"; YELLOW="#f7ecc9"

box(5, 12.2, 5.2, 0.9, "56,390 public genome assemblies\n(NCBI GenBank, A. baumannii)", BLUE, 12)
arrow(5, 11.75, 5, 11.15, "NCBI Datasets")
box(5, 10.7, 5.0, 0.9, "55,924 candidate assemblies\n(atypical records removed)", BLUE, 12)
arrow(5, 10.25, 5, 9.65, "CheckM2 QC")
box(5, 9.2, 6.0, 0.95, "Quality control and de-duplication\ncompleteness ≥90%, contamination ≤10%,\nsize 3.5–4.5 Mb, ≤500 contigs", PINK, 11)
arrow(5, 8.7, 5, 8.1, "")
box(5, 7.65, 4.2, 0.8, "48,968 analysis-set genomes", GREY, 13)
arrow(5, 7.25, 5, 6.65, "Bakta")
box(5, 6.2, 6.0, 0.9, "Uniform re-annotation\nsingle Bakta version and database,\nidentical settings for every genome", GREEN, 11)

# distribution bus to four tracks
ax.plot([5, 5], [5.75, 5.3], color="#333333", lw=1.6)
ax.plot([1.7, 8.3], [5.3, 5.3], color="#333333", lw=1.6)
xs = [1.7, 3.9, 6.1, 8.3]
tracks = ["Resistome\nAMRFinderPlus\nacquired genes",
          "Mobilome\nMOB-suite,\nISEScan, IntegronFinder",
          "Pan-genome and\nphylogeny\nPanaroo, PopPUNK, IQ-TREE",
          "Metadata\ncollection year\nand country"]
for x, t in zip(xs, tracks):
    ax.add_patch(FancyArrowPatch((x, 5.3), (x, 4.75), arrowstyle="-|>",
                 mutation_scale=15, lw=1.5, color="#333333"))
    box(x, 4.2, 2.0, 1.0, t, YELLOW, 10)

# convergence bus
ax.plot([1.7, 8.3], [3.55, 3.55], color="#333333", lw=1.6)
for x in xs:
    ax.plot([x, x], [3.7, 3.55], color="#333333", lw=1.5)
ax.add_patch(FancyArrowPatch((5, 3.55), (5, 3.05), arrowstyle="-|>",
             mutation_scale=18, lw=1.6, color="#333333"))
box(5, 2.5, 6.2, 1.0, "Contig-level linkage integration\nwhich mobile element carries which resistance gene,\nacross lineages and over time", PINK, 11)

fig.savefig("fig1_workflow.png", dpi=200, bbox_inches="tight", facecolor="white")
print("saved fig1_workflow.png")
