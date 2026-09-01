"""
seaice1d.plotstyle
-------------------
Shared matplotlib style for all figures: bold, dark, legible text and
axis labels; thicker lines; legends with a solid background so they never
visually blend into plotted data. Call apply_style() once at the top of
each figure script, and bold_ticks(ax) / legend_no_overlap(ax, ...) per
axes as needed.
"""
import matplotlib.pyplot as plt

DARK = "#111111"


def apply_style():
    plt.rcParams.update({
        "font.size": 12,
        "text.color": DARK,
        "axes.labelsize": 13,
        "axes.labelweight": "bold",
        "axes.labelcolor": DARK,
        "axes.titlesize": 14,
        "axes.titleweight": "bold",
        "axes.titlecolor": DARK,
        "axes.edgecolor": DARK,
        "axes.linewidth": 1.4,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "xtick.color": DARK,
        "ytick.color": DARK,
        "xtick.major.width": 1.3,
        "ytick.major.width": 1.3,
        "lines.linewidth": 2.3,
        "lines.markersize": 7,
        "legend.fontsize": 10.5,
        "legend.framealpha": 0.95,
        "legend.edgecolor": DARK,
        "legend.facecolor": "white",
        "figure.titlesize": 15,
        "figure.titleweight": "bold",
    })


def bold_ticks(ax):
    """Force tick label weight/color explicitly (not all mpl versions
    expose this via rcParams alone)."""
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_fontweight("bold")
        lbl.set_color(DARK)


def bold_legend_text(leg):
    for txt in leg.get_texts():
        txt.set_fontweight("bold")
        txt.set_color(DARK)
