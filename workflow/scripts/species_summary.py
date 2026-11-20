#!/usr/bin/env python3
"""Merge per-sample species counts, add ratios (and salmon/oarfish EM estimates), and plot."""
import argparse
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--counts", nargs="+", required=True)
ap.add_argument("--species", nargs=2, required=True, metavar=("HOST", "XENO"))
ap.add_argument("--salmon", nargs="*", default=[], help="quant.sf files, same order as --counts")
ap.add_argument("--oarfish", nargs="*", default=[], help="oarfish .quant files, same order as --counts")
ap.add_argument("--txmap", help="transcript -> species TSV (needed with --salmon/--oarfish)")
ap.add_argument("--title", default="")
ap.add_argument("--tsv", required=True)
ap.add_argument("--plot", required=True, help="output prefix; writes .pdf and .png")
args = ap.parse_args()
A, B = args.species

# ---------- table ----------
df = pd.concat([pd.read_csv(f, sep="\t") for f in args.counts], ignore_index=True)
df["sample"] = df["sample"].astype(str)
df["mapped"] = df[A] + df[B]
df["uniq_" + A] = df[A] - df["amb_" + A]
df["uniq_" + B] = df[B] - df["amb_" + B]
df["ambiguous"] = df["amb_" + A] + df["amb_" + B]


def ratio(x, y):
    return x / y if y > 0 else float("nan")


df["pct_mapped"] = 100 * df["mapped"] / df["total"]
df["pct_" + B + "_of_mapped"] = 100 * df[B] / df["mapped"]
df["ratio_" + B + "_" + A] = [ratio(x, y) for x, y in zip(df[B], df[A])]
df["ratio_uniq_" + B + "_" + A] = [ratio(x, y) for x, y in zip(df["uniq_" + B], df["uniq_" + A])]

# EM-based estimates: reads per species summed over each quantifier's transcript counts.
# quantifier name -> (files, transcript column, read-count column)
EM_TOOLS = {
    "salmon": (args.salmon, "Name", "NumReads"),
    "oarfish": (args.oarfish, "tname", "num_reads"),
}
EM_USED = [tool for tool, (files, *_) in EM_TOOLS.items() if files]
if EM_USED:
    txmap = pd.read_csv(args.txmap, sep="\t", header=None, names=["Name", "species"], index_col=0)["species"]
for tool in EM_USED:
    files, name_col, count_col = EM_TOOLS[tool]
    per_sample = {}
    for s, f in zip(df["sample"], files):
        q = pd.read_csv(f, sep="\t", usecols=[name_col, count_col])
        per_sample[s] = q.groupby(q[name_col].map(txmap))[count_col].sum()
    df[f"{tool}_{A}"] = [per_sample[s].get(A, 0.0) for s in df["sample"]]
    df[f"{tool}_{B}"] = [per_sample[s].get(B, 0.0) for s in df["sample"]]
    df[f"ratio_{tool}_{B}_{A}"] = [ratio(x, y) for x, y in zip(df[f"{tool}_{B}"], df[f"{tool}_{A}"])]

df.to_csv(args.tsv, sep="\t", index=False, float_format="%.6g")

# ---------- plot ----------
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
C_A, C_B, C_AMB, C_UNM = "#2a78d6", "#eb6834", "#1baf7a", "#c9c8c2"
plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.titlesize": 10,
    "axes.titleweight": "bold", "axes.titlelocation": "left", "figure.facecolor": SURF,
    "axes.facecolor": SURF, "savefig.facecolor": SURF,
})


def tidy(ax, xgrid=True):
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    if xgrid:
        ax.grid(axis="x", color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)


def human(n):
    for unit, d in (("M", 1e6), ("k", 1e3)):
        if abs(n) >= d:
            return f"{n / d:.1f}{unit}"
    return f"{n:.0f}"


samples = list(df["sample"])
n = len(samples)
y = list(range(n))[::-1]
h = 0.62
# fixed layout in inches: title, panel titles, plot rows, x-label + legend
ROW, TOP, BOT, LM, RM, GAP = 0.45, 0.85, 1.25, 0.9, 0.3, 0.35
W = 15.0
plot_h = max(1.2, ROW * n)
fig_h = TOP + plot_h + BOT
fig = plt.figure(figsize=(W, fig_h))
ratios = [1.35, 1, 1]
usable = W - LM - RM - 2 * GAP
widths = [usable * r / sum(ratios) for r in ratios]
axes, x0 = [], LM
for w in widths:
    axes.append(fig.add_axes([x0 / W, BOT / fig_h, w / W, plot_h / fig_h]))
    x0 += w + GAP
LEG_Y = -(0.62 / plot_h)  # legend 0.62 in below the axis, in axes units

# A: read composition (100% stacked)
ax = axes[0]
parts = [("uniq_" + A, f"{A} (unique)", C_A), ("ambiguous", "ambiguous (both)", C_AMB),
         ("uniq_" + B, f"{B} (unique)", C_B), ("unmapped", "unmapped", C_UNM)]
left = [0.0] * n
for col, lab, colr in parts:
    w = list(100 * df[col] / df["total"])
    ax.barh(y, w, left=left, height=h, color=colr, edgecolor=SURF, linewidth=2, label=lab)
    for yi, l, wi in zip(y, left, w):
        if wi >= 7:
            ax.text(l + wi / 2, yi, f"{wi:.0f}%", ha="center", va="center", fontsize=8,
                    color="white" if colr in (C_A, C_B) else INK)
    left = [l + wi for l, wi in zip(left, w)]
for yi, t in zip(y, df["total"]):
    ax.text(101, yi, human(t), va="center", fontsize=8, color=INK2)
ax.set_xlim(0, 109)
ax.set_xticks([0, 25, 50, 75, 100])
ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.0f}%"))
ax.set_yticks(y, samples)
ax.set_xlabel("share of reads (total reads at right)")
ax.set_title("Read composition")
tidy(ax)
ax.legend(loc="upper left", bbox_to_anchor=(0, LEG_Y), ncol=4, frameon=False, fontsize=8,
          handlelength=1.2, columnspacing=1.2)

# B: reads per species (log)
ax = axes[1]
bh = h / 2
ax.barh([v + bh / 2 for v in y], df[A].clip(lower=0.5), height=bh, color=C_A, edgecolor=SURF, linewidth=1.5, label=A)
ax.barh([v - bh / 2 for v in y], df[B].clip(lower=0.5), height=bh, color=C_B, edgecolor=SURF, linewidth=1.5, label=B)
ax.set_xscale("log")
for v, xa, xb in zip(y, df[A], df[B]):
    ax.text(max(xa, 0.5) * 1.15, v + bh / 2, human(xa), va="center", fontsize=7.5, color=INK2)
    ax.text(max(xb, 0.5) * 1.15, v - bh / 2, human(xb), va="center", fontsize=7.5, color=INK2)
lo = max(0.5, min(df[A].min(), df[B].min()) / 3) if min(df[A].min(), df[B].min()) > 0 else 0.5
ax.set_xlim(lo, max(df[A].max(), df[B].max()) * 6)
ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: human(v)))
ax.set_yticks(y, [""] * n)
ax.set_xlabel("mapped reads (primary alignment, log scale)")
ax.set_title("Reads per species")
tidy(ax)
ax.legend(loc="upper left", bbox_to_anchor=(0, LEG_Y), ncol=2, frameon=False, fontsize=8)

# C: ratio estimates (log2 B/A)
ax = axes[2]
est = [("ratio_" + B + "_" + A, "primary alignment", "o", INK)]
est.append(("ratio_uniq_" + B + "_" + A, "unique only", "o", "none"))
EM_MARKERS = {"salmon": ("D", INK2), "oarfish": ("s", "none")}
for tool in EM_USED:
    m, fc = EM_MARKERS[tool]
    est.append((f"ratio_{tool}_{B}_{A}", f"{tool} (EM)", m, fc))
vals = []
# spread the estimates of one sample evenly within its row
spread = 0.12 * (len(est) - 1)
offs = [spread / 2 - i * spread / (len(est) - 1) for i in range(len(est))]
for (col, lab, m, fc), off in zip(est, offs):
    lr = [math.log2(r) if r and r > 0 and not math.isnan(r) else float("nan") for r in df[col]]
    vals += [v for v in lr if not math.isnan(v)]
    ax.scatter(lr, [v + off for v in y], marker=m, s=46 if m == "o" else 34, facecolors=fc, edgecolors=INK,
               linewidths=1.2, zorder=3, label=lab)
for i, v in enumerate(y):
    row = [df[c].iloc[i] for c, *_ in est]
    row = [math.log2(r) for r in row if r and r > 0 and not math.isnan(r)]
    r = df["ratio_" + B + "_" + A].iloc[i]
    if row and r > 0:
        ax.annotate(f"{r:.3g}", (max(row), v), xytext=(9, 0), textcoords="offset points",
                    ha="left", va="center", fontsize=7.5, color=INK2)
ax.axvline(0, color=INK2, linewidth=0.8, linestyle=(0, (3, 3)))
lo_r, hi_r = min(vals + [0]), max(vals + [0])
pad = max(1.0, 0.15 * (hi_r - lo_r))
ax.set_xlim(lo_r - pad, hi_r + 1.6 * pad)
ax.set_ylim(-0.6, n - 0.4)
ax.set_yticks(y, [""] * n)
ax.set_xlabel(f"log2({B} / {A})   (labels: {B}/{A} ratio)")
ax.set_title(f"Species ratio {B} : {A}")
tidy(ax)
ax.legend(loc="upper left", bbox_to_anchor=(0, LEG_Y), ncol=2 if len(est) > 3 else 3, frameon=False, fontsize=8)

if args.title:
    fig.text(LM / W, 1 - 0.3 / fig_h, args.title, ha="left", va="center", fontsize=11, fontweight="bold")
for ext in ("pdf", "png"):
    fig.savefig(f"{args.plot}.{ext}", dpi=200)
