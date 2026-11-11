#!/usr/bin/env python3
"""Figure 7 - Circular core-genome ML phylogeny of 121 lineage representatives.
Tip size/colour scale with the number of collection genomes each lineage represents
(log scale). The eight lineages of >=490 genomes are tagged and keyed in a side panel.
Reads ../data/tree.treefile (Newick) and ../data/repmeta.tsv (cluster, rep_acc, cluster_size).
Requires: ete3 or dendropy for Newick parsing; this script uses the lightweight
Bio.Phylo reader from Biopython.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
from io import StringIO
from Bio import Phylo

meta = pd.read_csv("../data/repmeta.tsv", sep="\t")      # cluster, rep_acc, cluster_size
size_of = dict(zip(meta["rep_acc"], meta["cluster_size"]))
clust_of = dict(zip(meta["rep_acc"], meta["cluster"]))

tree = Phylo.read("../data/tree.treefile", "newick")
tips = tree.get_terminals()
n = len(tips)

# radial layout: assign each tip an angle; use branch-length depths for radius
angles = {t.name: 2*np.pi*i/n for i, t in enumerate(tips)}
depths = tree.depths()                                   # root distance per clade
maxd = max(depths.values()) or 1.0

fig = plt.figure(figsize=(16, 13))
ax = fig.add_axes([0.02, 0.05, 0.62, 0.9], polar=True)
ax.set_axis_off()

def xy(clade):
    # mean angle of descendant tips, radius from root depth
    ts = clade.get_terminals()
    a = np.mean([angles[t.name] for t in ts])
    r = depths[clade] / maxd
    return a, r

# draw branches
for clade in tree.find_clades(order="level"):
    if clade.clades:
        a0, r0 = xy(clade)
        child_ang = [xy(c)[0] for c in clade.clades]
        ax.plot([min(child_ang), max(child_ang)], [r0, r0], color="black", lw=1.4)
        for c in clade.clades:
            a1, r1 = xy(c)
            ax.plot([a1, a1], [r0, r1], color="black", lw=1.4)

# tips
sizes = meta["cluster_size"]
smax = np.log10(sizes.max())
tagged_meta = meta.sort_values("cluster_size", ascending=False).head(8).reset_index(drop=True)
tagged = tagged_meta["rep_acc"].tolist()
cmap = plt.get_cmap("viridis")
for t in tips:
    a = angles[t.name]; r = 1.0
    cs = size_of.get(t.name, 1)
    ms = 9 + 34*(np.log10(max(cs,1)) - 1)/(smax - 1 + 1e-9)
    col = cmap((np.log10(max(cs,1)))/smax)
    ax.plot(a, r, "o", ms=max(ms, 5), color=col, mec="white", mew=0.6, zorder=4)

ax.set_ylim(0, 1.15)

# side panel : the 8 tagged lineages
axp = fig.add_axes([0.66, 0.1, 0.32, 0.8]); axp.set_axis_off()
axp.set_title("Major lineages", fontsize=18, fontweight="bold", loc="left")
def lum(rgba):
    r, g, b = rgba[:3]
    return 0.299*r + 0.587*g + 0.114*b
for k, (_, row) in enumerate(tagged_meta.iterrows()):
    cs = row["cluster_size"]; cl = row["cluster"]
    yk = 1 - (k + 0.5)/len(tagged)
    col = cmap(np.log10(max(cs,1))/smax)
    axp.scatter(0.08, yk, s=900, color=col, edgecolor="white", zorder=3, transform=axp.transAxes)
    axp.text(0.08, yk, str(k+1), ha="center", va="center",
             color="white" if lum(col) < 0.5 else "black", fontsize=15, fontweight="bold",
             transform=axp.transAxes, zorder=4)
    name = "IC2" if cl == 1 else ("IC1" if cl == 7 else f"PP-{cl}")
    axp.text(0.2, yk, f"{name}  (n = {cs:,})", va="center", fontsize=15, transform=axp.transAxes)

fig.text(0.33, 0.02, "2,561 core genes  ·  GTR+F+G4  ·  1000 ultrafast bootstrap replicates  "
                     "·  representatives span 46,483/48,968 genomes (94.9%)",
         ha="center", fontsize=14, color="black")
fig.savefig("fig7_phylo.png", dpi=200, bbox_inches="tight", facecolor="white")
print("saved fig7_phylo.png")
