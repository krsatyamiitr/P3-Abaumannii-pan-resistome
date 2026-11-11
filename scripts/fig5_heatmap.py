#!/usr/bin/env python3
"""Figure 5 - Genomic-context matrix for the major acquired resistance genes.
Cell shading = percentage of gene-carrying contigs in each context
(plasmid, chromosome, IS-associated, integron). Right-hand column = most
frequently co-located IS family.
Reads ../data/linkage_top.tsv
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd, numpy as np

df = pd.read_csv("../data/linkage_top.tsv", sep="\t")
df = df.sort_values("total", ascending=False).reset_index(drop=True)
genes = df["gene"].values
cols = ["pct_plasmid", "pct_chromosome", "pct_IS", "pct_integron"]
labels = ["Plasmid", "Chromosome", "IS-\nassociated", "Integron"]
M = df[cols].values

fig, ax = plt.subplots(figsize=(10, 12))
plt.subplots_adjust(left=0.22, right=0.86, top=0.95, bottom=0.08)
im = ax.imshow(M, aspect="auto", cmap="RdYlBu_r", vmin=0, vmax=100,
               extent=[-0.5, 3.5, len(genes) - 0.5, -0.5])

ax.set_xticks(range(4)); ax.set_xticklabels(labels, fontsize=14)
ax.set_yticks(range(len(genes)))
ax.set_yticklabels(genes, fontsize=12, fontstyle="italic")
ax.set_xlim(-0.6, 5.0); ax.set_ylim(len(genes) - 0.5, -2.2)

# cell numbers
for i in range(len(genes)):
    for j in range(4):
        v = M[i, j]
        ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                fontsize=9, color="white" if (v > 60 or v < 15) else "black")

# IS family column at x=4.3
ax.text(4.3, -1.0, "Top IS\nfamily", ha="center", va="center", fontsize=13, fontweight="bold")
for i, fam in enumerate(df["top_IS_family"]):
    ax.text(4.3, i, str(fam), ha="center", va="center", fontsize=11)

cb = fig.colorbar(im, ax=ax, fraction=0.040, pad=0.02)
cb.set_label("Percent of gene-carrying contigs in context", fontsize=13)
cb.ax.tick_params(labelsize=11)

fig.savefig("fig5_heatmap.png", dpi=200, bbox_inches="tight", facecolor="white")
print("saved fig5_heatmap.png")
