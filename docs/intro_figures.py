"""Generate conceptual vector figures for the executable article.

These schematics explain model scope and planned experimental evaluation;
they do not depict measured reaction rates or optimization performance.
"""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


BLUE = "#e3eff9"
GREEN = "#e5f2e7"
AMBER = "#fff1d6"
INK = "#24384a"


def box(ax, x, y, w, h, title, body="", color=BLUE, dashed=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.015",
                              facecolor=color, edgecolor=INK, linewidth=1.1,
                              linestyle="--" if dashed else "-"))
    ax.text(x+w/2, y+h*(0.70 if body else 0.5), title, ha="center",
            va="center", fontsize=10, weight="bold", color=INK)
    if body:
        ax.text(x+w/2, y+h*0.32, body, ha="center", va="center",
                fontsize=8.5, color=INK, linespacing=1.35)


def arrow(ax, start, end, label="", dashed=False):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>",
                                mutation_scale=13, linewidth=1.1, color=INK,
                                linestyle="--" if dashed else "-"))
    if label:
        ax.text((start[0]+end[0])/2, (start[1]+end[1])/2+0.13, label,
                ha="center", va="bottom", fontsize=8, color=INK)


def canvas(height):
    fig, ax = plt.subplots(figsize=(10.5, height))
    ax.set(xlim=(0, 10), ylim=(0, height))
    ax.axis("off")
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.03, top=0.97)
    return fig, ax


def build_figures(output="figures"):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    fig, ax = canvas(5.4)
    box(ax, 0.15, 3.45, 1.75, 1.05, "DNA template", "T7 promoter\nassumed transcript bounds")
    box(ax, 2.7, 3.45, 1.7, 1.05, "mRNA", "transcription\nand degradation")
    box(ax, 5.2, 3.45, 1.75, 1.05, "Protein", "translation\nimmature product")
    box(ax, 7.75, 3.45, 2.0, 1.05, "Target readout", "GFP maturation\nenzyme activity: planned", GREEN)
    arrow(ax, (1.9, 3.98), (2.7, 3.98), "transcription")
    arrow(ax, (4.4, 3.98), (5.2, 3.98), "translation")
    arrow(ax, (6.95, 3.98), (7.75, 3.98), "observation")
    box(ax, 0.15, 1.55, 4.25, 1.1, "Shared biochemical resources",
        "ATP / GTP / CTP / UTP; 20 amino-acid pools\nfuel-driven ATP regeneration; machinery decay")
    arrow(ax, (3.55, 2.65), (3.55, 3.45))
    arrow(ax, (4.4, 2.35), (6.0, 3.45))
    box(ax, 5.2, 1.55, 4.55, 1.1, "Environment and folding: future modules",
        "free magnesium; pH; oxygen transfer\nredox; chaperones; aggregation; active fraction", AMBER, True)
    arrow(ax, (7.0, 2.65), (6.7, 3.45), dashed=True)
    arrow(ax, (8.75, 2.65), (8.75, 3.45), dashed=True)
    ax.text(5, 0.65, "Sequence determines substrate demand; the assay determines what performance means.",
            ha="center", fontsize=10, color=INK)
    ax.text(5, 0.2, "Solid resource modules: implemented reduced model. Dashed modules: proposed extensions.",
            ha="center", fontsize=8.5, color=INK)
    fig.savefig(output / "intro-resources.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, ax = canvas(5.1)
    box(ax, 0.15, 3.45, 2.15, 1.15, "Evidence", "recipes and constructs\nsource / units / batch metadata")
    box(ax, 2.7, 3.45, 2.15, 1.15, "Feasibility", "stock transfers and final pools\nfixed base and volume limits")
    box(ax, 5.25, 3.45, 2.15, 1.15, "Model hypotheses", "mechanistic states\nuncertain effective parameters")
    box(ax, 7.8, 3.45, 2.05, 1.15, "Select a batch", "empirical / mechanistic\ncompare equal well budgets", GREEN)
    for start, end in [(2.3, 2.7), (4.85, 5.25), (7.4, 7.8)]:
        arrow(ax, (start, 4.0), (end, 4.0))
    ax.text(5.05, 4.88, "Recipe-to-model mapping and calibrated selection remain to be implemented.",
            ha="center", fontsize=8.5, color=INK)
    box(ax, 7.8, 0.8, 2.05, 1.3, "Experiment", "randomized module swaps\nreplicates and common controls", AMBER)
    box(ax, 4.15, 0.8, 2.9, 1.3, "Calibrate and test", "train on complete conditions\nwithhold recipes / independent batches", AMBER)
    box(ax, 0.15, 0.8, 3.25, 1.3, "Confirm conclusions", "independent winner retests\nerror, interval coverage, target activity", AMBER)
    arrow(ax, (8.82, 3.45), (8.82, 2.1))
    ax.text(9.15, 2.8, "planned", ha="center", fontsize=8, color=INK)
    arrow(ax, (7.8, 1.45), (7.05, 1.45))
    arrow(ax, (4.15, 1.45), (3.4, 1.45))
    arrow(ax, (5.6, 2.1), (6.3, 3.45), dashed=True)
    ax.text(6.7, 2.8, "update model", ha="center", fontsize=8, color=INK)
    ax.text(5, 0.15, "Synthetic recovery tests the software; prospective measurements test biological usefulness.",
            ha="center", fontsize=9, color=INK)
    fig.savefig(output / "intro-workflow.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    build_figures()
