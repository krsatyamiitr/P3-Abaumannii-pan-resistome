#!/usr/bin/env python3
"""Figure 3 - Co-carriage (UpSet) of the fifteen most prevalent acquired resistance genes.
Upper bars: number of genomes carrying each gene combination (top 15 combinations).
Matrix: which genes make up each combination. Left bars: total genomes carrying each gene.
Reads ../data/upset_matrix.tsv  (acc + 15 binary gene columns).
Built manually with matplotlib GridSpec (no upsetplot dependency).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd, numpy as np

df = pd.read_csv("../data/upset_matrix.tsv", sep="\t")
genes = [c for c in df.columns if c != "acc"]
B = df[genes].astype(int)

# total genomes per gene (left bars)
totals = B.sum(axis=0).sort_values(ascending=False)
genes = list(totals.index)                     # order rows by overall prevalence
B = B[genes]

# combinations
combo = B.apply(lambda r: tuple(r.values), axis=1)
counts = combo.value_counts().head(15)         # top 15 combinations
combos = list(counts.index)
cvals = counts.values

NG = len(genes); NC = len(combos)
fig = plt.figure(figsize=(14, 10))
gs = gridspec.GridSpec(2, 2, width_ratios=[1.0, 4.0], height_ratios=[2.2, 3.0],
                       wspace=0.05, hspace=0.05)

# top bars : combination sizes
ax_top = fig.add_subplot(gs[0, 1])
x = np.arange(NC)
ax_top.bar(x, cvals, color="#2c3e7b", width=0.72)
for xi, v in zip(x, cvals):
    ax_top.text(xi, v + max(cvals)*0.01, f"{v:,}", ha="center", va="bottom", fontsize=9)
ax_top.set_xlim(-0.6, NC - 0.4); ax_top.set_xticks([])
ax_top.set_ylabel("Genomes with\ncombination", fontsize=13)
for s in ["top", "right"]: ax_top.spines[s].set_visible(False)

# left bars : gene totals
ax_left = fig.add_subplot(gs[1, 0])
yy = np.arange(NG)
ax_left.barh(yy, [totals[g] for g in genes], color="#9aa7c7", height=0.6)
ax_left.set_ylim(NG - 0.5, -0.5); ax_left.set_yticks([])
ax_left.invert_xaxis()
ax_left.set_xlabel("Genomes\ncarrying gene", fontsize=12)
for s in ["top", "left"]: ax_left.spines[s].set_visible(False)

# matrix
ax_mat = fig.add_subplot(gs[1, 1])
for yi in range(NG):
    ax_mat.scatter(range(NC), [yi]*NC, s=70, color="#d0d0d0", zorder=1)
for ci, cb in enumerate(combos):
    members = [i for i, v in enumerate(cb) if v == 1]
    ax_mat.scatter([ci]*len(members), members, s=90, color="#2c3e7b", zorder=3)
    if len(members) > 1:
        ax_mat.plot([ci, ci], [min(members), max(members)], color="#2c3e7b", lw=2, zorder=2)
ax_mat.set_xlim(-0.6, NC - 0.4); ax_mat.set_ylim(NG - 0.5, -0.5)
ax_mat.set_xticks([])
ax_mat.set_yticks(range(NG))
ax_mat.set_yticklabels(genes, fontsize=11, fontstyle="italic")
ax_mat.yaxis.tick_right()
for s in ["top", "right", "bottom"]: ax_mat.spines[s].set_visible(False)

fig.savefig("fig3_upset.png", dpi=200, bbox_inches="tight", facecolor="white")
print("saved fig3_upset.png")
