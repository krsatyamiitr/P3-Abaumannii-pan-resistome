#!/usr/bin/env python3
"""Figure 4 - Linkage architecture of the 25 most prevalent acquired resistance genes.
(A) Number of genomes carrying each gene.
(B) Plasmid vs chromosome assignment of gene-carrying contigs (MOB-suite).
(C) IS-association and integron-association of gene-carrying contigs.
Reads ../data/linkage_top.tsv
Columns: gene, total, pct_plasmid, pct_chromosome, pct_integron, pct_IS, top_replicon, top_IS_family
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd, numpy as np

df = pd.read_csv("../data/linkage_top.tsv", sep="\t")
df = df.sort_values("total", ascending=True).reset_index(drop=True)
genes = df["gene"].values
y = np.arange(len(df))

fig, axes = plt.subplots(1, 3, figsize=(20, 9), gridspec_kw={"width_ratios": [1.0, 1.1, 0.9]})
plt.subplots_adjust(wspace=0.08, left=0.13, right=0.98, top=0.93, bottom=0.09)
axA, axB, axC = axes

# ---- A : prevalence (thousands) ----
axA.barh(y, df["total"] / 1000.0, color="#3b5a9a", edgecolor="white", height=0.74)
axA.set_yticks(y); axA.set_yticklabels(genes, fontsize=13, fontstyle="italic")
axA.set_xlabel("Genomes carrying (thousands)", fontsize=17, fontweight="bold")
axA.tick_params(axis="x", labelsize=14)
for yi, v in zip(y, df["total"]):
    lab = f"{v/1000:.0f}k" if v >= 1000 else str(int(v))
    axA.text(v/1000.0 + 0.5, yi, lab, va="center", fontsize=11)
for s in ["top", "right"]: axA.spines[s].set_visible(False)
axA.set_title("A  Prevalence", loc="left", fontsize=18, fontweight="bold")

# ---- B : plasmid vs chromosome (stacked) ----
axB.barh(y, df["pct_plasmid"], color="#c0392b", edgecolor="white", height=0.74, label="Plasmid")
axB.barh(y, df["pct_chromosome"], left=df["pct_plasmid"], color="#27ae60",
         edgecolor="white", height=0.74, label="Chromosome")
axB.set_yticks([]); axB.set_xlim(0, 100)
axB.set_xlabel("Replicon context (%)", fontsize=17, fontweight="bold")
axB.tick_params(axis="x", labelsize=14)
for s in ["top", "right", "left"]: axB.spines[s].set_visible(False)
axB.legend(loc="lower center", bbox_to_anchor=(0.5, -0.14), ncol=2, fontsize=13, frameon=False)
axB.set_title("B  Plasmid vs chromosome location", loc="left", fontsize=18, fontweight="bold")

# ---- C : IS- and integron-association (dot sizes) ----
axC.set_yticks([]); axC.set_xlim(-0.5, 1.5); axC.set_ylim(-0.7, len(df) - 0.3)
xIS, xINT = 0.0, 1.0
smax = max(df["pct_IS"].max(), df["pct_integron"].max(), 1)
for yi, is_v, int_v in zip(y, df["pct_IS"], df["pct_integron"]):
    axC.scatter(xIS, yi, s=30 + 9*is_v, color="#c0392b", alpha=0.85, edgecolor="white", zorder=3)
    axC.text(xIS, yi, f"{is_v:.0f}", ha="center", va="center", fontsize=8, color="white", zorder=4)
    axC.scatter(xINT, yi, s=30 + 9*int_v, color="#8e44ad", alpha=0.85, edgecolor="white", zorder=3)
    axC.text(xINT, yi, f"{int_v:.0f}", ha="center", va="center", fontsize=8, color="white", zorder=4)
axC.set_xticks([xIS, xINT]); axC.set_xticklabels(["IS-\nassociated", "Integron-\nassociated"], fontsize=14)
for s in ["top", "right", "left"]: axC.spines[s].set_visible(False)
axC.set_title("C  Mobile-element association", loc="left", fontsize=18, fontweight="bold")

fig.text(0.5, 0.015,
         "Percentages are the fraction of gene-carrying contigs assigned to each context; "
         "a gene may be IS- and integron-associated on the same contig.",
         ha="center", fontsize=11, color="black")

fig.savefig("fig4_linkage.png", dpi=200, bbox_inches="tight", facecolor="white")
print("saved fig4_linkage.png")
