import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import pandas as pd, numpy as np

# ---- data ----
y = pd.read_csv("../data/yearly.tsv", sep="\t")
y = y[y["n_genomes"] >= 100].copy()        # drop sparse early years
c = pd.read_csv("../data/country.tsv", sep="\t")
c = c.sort_values("pct_carbapenemase", ascending=True).reset_index(drop=True)

RED   = "#c0392b"
ORANGE= "#e8745a"
BLUE  = "#2b6cb0"

fig, (axA, axB) = plt.subplots(1, 2, figsize=(21, 8.6), gridspec_kw={"width_ratios":[1.06,1.0]})
plt.subplots_adjust(wspace=0.40, left=0.055, right=0.965, top=0.93, bottom=0.12)

# ============ PANEL A : temporal ============
axA.set_zorder(2); axA.patch.set_alpha(0)
axA.fill_between(y["year"], y["carbapenemase"], color=RED, alpha=0.08, zorder=0)
axA.plot(y["year"], y["carbapenemase"], "-o", color=RED, lw=3.0, ms=9,
         mec="white", mew=1.0, label="Carbapenemase-positive", zorder=4)
axA.plot(y["year"], y["blaOXA-23"], "-s", color=ORANGE, lw=2.6, ms=8,
         mec="white", mew=1.0, label="blaOXA-23", zorder=3)
axA.set_xlabel("Collection year", fontsize=24, fontweight="bold")
axA.set_ylabel("Genomes positive (%)", fontsize=24, fontweight="bold", color=RED)
axA.set_ylim(0, 100); axA.set_xlim(y["year"].min()-1.2, y["year"].max()+1.2)
axA.tick_params(axis="y", labelsize=19, colors=RED)
axA.tick_params(axis="x", labelsize=19)
for s in ["top"]: axA.spines[s].set_visible(False)
axA.spines["left"].set_color(RED); axA.spines["left"].set_linewidth(1.8)
axA.spines["bottom"].set_linewidth(1.6)
axA.grid(axis="y", ls=":", color="0.7", lw=0.8, zorder=0)

# right axis : mean AMR load
axR = axA.twinx(); axR.set_zorder(1)
axR.plot(y["year"], y["mean_amr_load"], "-^", color=BLUE, lw=2.6, ms=8.5,
         mec="white", mew=1.0, label="Mean resistance-gene load", zorder=3)
axR.set_ylabel("Mean acquired AMR genes / genome", fontsize=24, fontweight="bold", color=BLUE, labelpad=10)
axR.set_ylim(0, max(y["mean_amr_load"].max()*1.08, 21))
axR.tick_params(axis="y", labelsize=19, colors=BLUE)
axR.spines["right"].set_color(BLUE); axR.spines["right"].set_linewidth(1.8)
axR.spines["top"].set_visible(False)

# combined legend
h1,l1 = axA.get_legend_handles_labels()
h2,l2 = axR.get_legend_handles_labels()
leg = axA.legend(h1+h2, l1+l2, loc="lower right", fontsize=18, frameon=True,
                 framealpha=0.95, edgecolor="0.6", borderpad=0.8,
                 handlelength=2.2, labelspacing=0.6)
leg.set_zorder(10)

# footnote (black, larger)
axA.text(0.015, 0.045, "Years <2003 (n<100) omitted", transform=axA.transAxes,
         fontsize=17, style="italic", color="black", ha="left", va="bottom")

# panel tag
axA.text(-0.085, 1.02, "A", transform=axA.transAxes, fontsize=30,
         fontweight="bold", va="bottom", ha="left")

# ============ PANEL B : geographic ============
vals = c["pct_carbapenemase"].values
ypos = np.arange(len(c))
cmap = plt.get_cmap("RdYlBu_r")
norm = (vals - vals.min()) / (vals.max() - vals.min())
colors = [cmap(0.12 + 0.76*v) for v in norm]
axB.barh(ypos, vals, color=colors, edgecolor="white", lw=0.8, height=0.72, zorder=3)

axB.set_yticks(ypos)
axB.set_yticklabels(c["country"], fontsize=18)
axB.set_ylim(-0.7, len(c)-0.3)
axB.set_xlim(0, 112)
axB.set_xlabel("Carbapenemase-positive genomes (%)", fontsize=24, fontweight="bold")
axB.tick_params(axis="x", labelsize=19)
for s in ["top","right"]: axB.spines[s].set_visible(False)
axB.spines["left"].set_linewidth(1.6); axB.spines["bottom"].set_linewidth(1.6)
axB.grid(axis="x", ls=":", color="0.75", lw=0.8, zorder=0)

# n= labels to the right of each bar
for yi, v, n in zip(ypos, vals, c["n_genomes"]):
    axB.text(v+1.6, yi, f"n={n:,}", va="center", ha="left",
             fontsize=15.5, color="0.15")

axB.text(-0.20, 1.02, "B", transform=axB.transAxes, fontsize=30,
         fontweight="bold", va="bottom", ha="left")

fig.savefig("fig8_trends.png", dpi=200, facecolor="white")
print("saved")
