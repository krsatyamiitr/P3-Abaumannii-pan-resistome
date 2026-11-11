#!/usr/bin/env python3
"""Figure 6 - Tripartite association network linking acquired resistance genes
to their mobile-element contexts.
Resistance genes (centre) joined to the most frequently co-located IS family (left)
and plasmid replicon type (right). Edge thickness reflects association strength;
gene nodes coloured by dominant genomic context.
Reads ../data/linkage_top.tsv
Columns used: gene, pct_plasmid, pct_chromosome, pct_IS, top_replicon, top_IS_family
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd, numpy as np

df = pd.read_csv("../data/linkage_top.tsv", sep="\t")
df = df.sort_values("total", ascending=False).reset_index(drop=True)

genes = df["gene"].tolist()
is_fams = [f for f in dict.fromkeys(df["top_IS_family"]) if f and f != "-"]
reps = [r for r in dict.fromkeys(df["top_replicon"]) if r and r != "-"]

def col_of(row):
    if row["pct_plasmid"] >= 50: return "#2b6cb0"      # plasmid-borne
    if row["pct_IS"] >= 50:      return "#c0392b"       # strongly IS-associated
    return "#27ae60"                                     # mostly chromosomal

fig, ax = plt.subplots(figsize=(14, 12)); ax.set_axis_off()
xL, xC, xR = 0.0, 1.0, 2.0
yg = {g: i for i, g in enumerate(reversed(genes))}
yi = {f: j*(len(genes)/max(len(is_fams),1)) for j, f in enumerate(reversed(is_fams))}
yr = {r: j*(len(genes)/max(len(reps),1)) for j, r in enumerate(reversed(reps))}

# edges
for _, row in df.iterrows():
    g = row["gene"]; gy = yg[g]
    fam = row["top_IS_family"]; rep = row["top_replicon"]
    if fam in yi:
        ax.plot([xL, xC], [yi[fam], gy], color="#c0392b",
                lw=0.5 + 4*row["pct_IS"]/100.0, alpha=0.45, zorder=1)
    if rep in yr:
        ax.plot([xC, xR], [gy, yr[rep]], color="#2b6cb0",
                lw=0.5 + 4*row["pct_plasmid"]/100.0, alpha=0.45, zorder=1)

# nodes
for f, y in yi.items():
    ax.scatter(xL, y, s=320, color="#c0392b", marker="s", edgecolor="white", zorder=3)
    ax.text(xL-0.08, y, f, ha="right", va="center", fontsize=13, fontweight="bold")
for g, y in yg.items():
    row = df[df["gene"] == g].iloc[0]
    ax.scatter(xC, y, s=300, color=col_of(row), edgecolor="white", zorder=3)
    ax.text(xC+0.06, y, g, ha="left", va="center", fontsize=12, fontstyle="italic")
for r, y in yr.items():
    ax.scatter(xR, y, s=320, color="#2b6cb0", marker="D", edgecolor="white", zorder=3)
    ax.text(xR+0.08, y, r, ha="left", va="center", fontsize=13, fontweight="bold")

ax.text(xL, len(genes)+0.5, "Insertion-sequence\nfamily", ha="center", fontsize=15, fontweight="bold", color="#c0392b")
ax.text(xC, len(genes)+0.5, "Acquired resistance gene", ha="center", fontsize=15, fontweight="bold")
ax.text(xR, len(genes)+0.5, "Plasmid replicon\ntype", ha="center", fontsize=15, fontweight="bold", color="#2b6cb0")

from matplotlib.lines import Line2D
leg = [Line2D([0],[0], marker='o', color='w', markerfacecolor='#2b6cb0', markersize=12, label='Gene mostly plasmid-borne'),
       Line2D([0],[0], marker='o', color='w', markerfacecolor='#27ae60', markersize=12, label='Gene mostly chromosomal'),
       Line2D([0],[0], marker='o', color='w', markerfacecolor='#c0392b', markersize=12, label='Gene strongly IS-associated')]
ax.legend(handles=leg, loc="lower center", bbox_to_anchor=(0.5,-0.04), ncol=3, fontsize=12, frameon=False)
ax.set_xlim(-0.6, 2.8); ax.set_ylim(-1.2, len(genes)+1.2)

fig.savefig("fig6_network.png", dpi=200, bbox_inches="tight", facecolor="white")
print("saved fig6_network.png")
