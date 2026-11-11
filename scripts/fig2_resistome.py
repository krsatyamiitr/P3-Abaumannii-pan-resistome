#!/usr/bin/env python3
"""Figure 2 - Acquired resistome of the collection.
(A) Percentage of genomes carrying at least one acquired gene of each drug class.
(B) The twenty most prevalent individual acquired resistance genes, coloured by class.
Reads ../data/resistome_class.tsv and ../data/resistome_gene.tsv.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd, numpy as np

cls = pd.read_csv("../data/resistome_class.tsv", sep="\t")   # class, n_genomes, pct
gen = pd.read_csv("../data/resistome_gene.tsv", sep="\t")    # gene, n_genomes, pct, class

cls = cls.sort_values("pct", ascending=True)
gen = gen.sort_values("pct", ascending=True).tail(20)

# colour map by drug class
classes = sorted(set(cls["class"]).union(set(gen["class"])))
palette = plt.get_cmap("tab10")
cmap = {c: palette(i % 10) for i, c in enumerate(classes)}

fig, (axA, axB) = plt.subplots(1, 2, figsize=(18, 8.5))
plt.subplots_adjust(wspace=0.42, left=0.15, right=0.97, top=0.93, bottom=0.10)

# ---- Panel A : class prevalence ----
yA = np.arange(len(cls))
axA.barh(yA, cls["pct"], color=[cmap[c] for c in cls["class"]],
         edgecolor="white", height=0.72)
axA.set_yticks(yA); axA.set_yticklabels([c.title() for c in cls["class"]], fontsize=17)
axA.set_xlim(0, 108); axA.set_xlabel("Genomes with ≥1 gene (%)", fontsize=20, fontweight="bold")
axA.tick_params(axis="x", labelsize=16)
for yi, v in zip(yA, cls["pct"]):
    axA.text(v + 1.5, yi, f"{v:.0f}%", va="center", fontsize=14)
for s in ["top", "right"]: axA.spines[s].set_visible(False)
axA.set_title("A", loc="left", fontsize=26, fontweight="bold")

# ---- Panel B : top genes ----
yB = np.arange(len(gen))
axB.barh(yB, gen["pct"], color=[cmap[c] for c in gen["class"]],
         edgecolor="white", height=0.72)
axB.set_yticks(yB)
axB.set_yticklabels(gen["gene"], fontsize=14, fontstyle="italic")
axB.set_xlim(0, 108); axB.set_xlabel("Genomes carrying gene (%)", fontsize=20, fontweight="bold")
axB.tick_params(axis="x", labelsize=16)
for yi, v in zip(yB, gen["pct"]):
    axB.text(v + 1.5, yi, f"{v:.0f}", va="center", fontsize=12)
for s in ["top", "right"]: axB.spines[s].set_visible(False)
axB.set_title("B", loc="left", fontsize=26, fontweight="bold")

handles = [plt.Rectangle((0, 0), 1, 1, color=cmap[c]) for c in classes]
axB.legend(handles, [c.title() for c in classes], title="Drug class",
           fontsize=12, title_fontsize=13, loc="lower right", frameon=True)

fig.savefig("fig2_resistome.png", dpi=200, bbox_inches="tight", facecolor="white")
print("saved fig2_resistome.png")
