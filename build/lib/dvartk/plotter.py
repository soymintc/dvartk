import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2, venn2_unweighted


def plot_sv_spectra(
    counts, title, save_path=False, tag="", yscale_log=False, ylim=None, debug=False
):
    """Draw SV spectra plot based on SV counts"""

    sv_colors = [
        "#d6e6f4",  # deletion
        "#abd0e6",
        "#6aaed6",
        "#3787c0",
        "#105ba4",
        "#08315c",
        "#fedfc0",  # duplication
        "#fdb97d",
        "#fd8c3b",
        "#e95e0d",
        "#b63c02",
        "#642101",
        "#dbcce8",  # insertion
        "#b799d2",
        "#9366bc",
        "#6f4298",
        "#4a2c65",
        "#39224f",
        "#dbf1d6",  # inversion
        "#aedea7",
        "#73c476",
        "#37a055",
        "#0b7734",
        "#043316",
        "#aaaaaa",  # translocation
    ]

    fig, ax = plt.subplots(1)
    fig.set_figheight(2.5)
    fig.set_figwidth(8)

    if tag:
        title += " " + tag
    fig.suptitle(title)

    font = matplotlib.font_manager.FontProperties()
    font.set_family("monospace")

    ax.bar(height=counts, x=range(counts.shape[0]), color=sv_colors)

    xaxis_index = range(len(counts))
    if debug:
        print(xaxis_index, counts)
    plt.xticks(xaxis_index, counts.index, rotation=90, fontproperties=font)
    plt.xlim((-1, len(counts)))

    if yscale_log:
        plt.yscale("log")
    if ylim:
        assert len(ylim) == 2
        plt.ylim(ylim)

    if save_path:
        plt.savefig(save_path)


def plot_snv_spectra(
    counts, title, save_path=False, tag="", yscale_log=False, ylim=None, debug=False
):
    """Draw SNV spectra plot based on SNV counts"""

    fig, ax = plt.subplots(1)
    fig.set_figheight(3)
    fig.set_figwidth(20)

    if tag:
        title += tag
    fig.suptitle(title)

    pat = r"([ACGT])\[([ACGT])\>([ACGT])\]([ACGT])"
    df = pd.DataFrame(counts.copy(), columns=["count"])
    df["norm_tri_nt"] = df.index.str.replace(pat, r"\1\2\4", regex=True)
    df["norm_mut_type"] = df.index.str.replace(pat, r"\2>\3", regex=True)
    df["index"] = range(df.shape[0])
    if debug:
        print(df)

    for mut_type, mut_type_data in df.groupby(["norm_mut_type"]):
        if debug:
            print(mut_type, mut_type_data)
        ax.bar(data=mut_type_data, x="index", height="count", label=mut_type)

    font = matplotlib.font_manager.FontProperties()
    font.set_family("monospace")

    ax.set_xticks(df["index"])
    ax.set_xticklabels(df["norm_tri_nt"], rotation=90, fontproperties=font)
    ax.set_xlim((-1, 97))
    ax.legend(bbox_to_anchor=(1, 1), loc="upper left")

    sns.despine(ax=ax, trim=True)
    plt.tight_layout()

    if yscale_log:
        plt.yscale("log")
    if ylim:
        assert len(ylim) == 2
        plt.ylim(ylim)

    if save_path:
        plt.savefig(save_path)


def plot_venn2(cmp, weighted=False, label1="A", label2="B"):
    """Draw a venn diagram from a Snv/SvComparison instance"""
    if weighted:
        vd = venn2([cmp.A, cmp.B], set_labels=(label1, label2))
    else:
        vd = venn2_unweighted([cmp.A, cmp.B], set_labels=(label1, label2))
