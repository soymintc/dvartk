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
        "#d6e6f4",
        "#abd0e6",
        "#6aaed6",
        "#3787c0",
        "#105ba4",
        "#08315c",  # deletion
        "#fedfc0",
        "#fdb97d",
        "#fd8c3b",
        "#e95e0d",
        "#b63c02",
        "#642101",  # duplication
        "#dbcce8",
        "#b799d2",
        "#9366bc",
        "#6f4298",
        "#4a2c65",
        "#39224f",  # insertion
        "#dbf1d6",
        "#aedea7",
        "#73c476",
        "#37a055",
        "#0b7734",
        "#043316",  # inversion
        "#aaaaaa",
    ]  # translocation

    fig, ax = plt.subplots(1)
    fig.set_figheight(3.5)
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
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path)
        plt.close()


def plot_snv_spectra(
    counts, title, save_path=False, tag="", yscale_log=False, ylim=None, debug=False
):
    """Draw SNV spectra plot based on SNV counts"""

    colors = ["#03BDEE", "#000000", "#E52A25", "#CDC9CA", "#A3CE62", "#ECC6C5"]
    fig, ax = plt.subplots(1)
    fig.set_figheight(4)
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

    for ix, (mut_type, mut_type_data) in enumerate(df.groupby(["norm_mut_type"])):
        color = colors[ix]
        if debug:
            print(mut_type, mut_type_data)
        ax.bar(
            data=mut_type_data, x="index", height="count", label=mut_type, color=color
        )

    font = matplotlib.font_manager.FontProperties()
    font.set_family("monospace")

    ax.set_xticks(df["index"])
    ax.set_xticklabels(df["norm_tri_nt"], rotation=90, fontproperties=font)
    ax.set_xlim((-1, 97))
    ax.legend(bbox_to_anchor=(1, 1), loc="upper left")

    sns.despine(ax=ax, trim=True)

    if yscale_log:
        plt.yscale("log")
    if ylim:
        assert len(ylim) == 2
        plt.ylim(ylim)
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path)
        plt.close()


def proc_indel_dataframe(indel):
    """indel: DataFrame (not Series) with SigProfiler index and a 'count' column
    returns: annotated df
    """
    df = indel.copy()

    indel_feature_to_type = {
        "1:Del:C": "1bp Del at Homopolymer C",
        "1:Del:T": "1bp Del at Homopolymer T",
        "1:Ins:C": "1bp Ins at Homopolymer C",
        "1:Ins:T": "1bp Ins at Homopolymer T",
        "2:Del:R": "2bp Del at Repeats",
        "3:Del:R": "3bp Del at Repeats",
        "4:Del:R": "4bp Del at Repeats",
        "5:Del:R": "5+bp Del at Repeats",
        "2:Ins:R": "2bp Ins at Repeats",
        "3:Ins:R": "3bp Ins at Repeats",
        "4:Ins:R": "4bp Ins at Repeats",
        "5:Ins:R": "5+bp Ins at Repeats",
        "2:Del:M": "2bp Del at Microhomology",
        "3:Del:M": "3bp Del at Microhomology",
        "4:Del:M": "4bp Del at Microhomology",
        "5:Del:M": "5+bp Del at Microhomology",
    }

    pat = r"(\d):([A-Z][a-z]{2}):([CTRM]):(\d)"
    df.loc[:, "indel_length"] = df.index.str.replace(pat, r"\1", regex=True).astype(str)
    df.loc[:, "insdel"] = df.index.str.replace(pat, r"\2", regex=True)
    df.loc[:, "feature_name"] = df.index.str.replace(pat, r"\3", regex=True)
    df.loc[:, "feature_length"] = df.index.str.replace(pat, r"\4", regex=True).astype(
        int
    )
    df.loc[:, "feature"] = (
        df["indel_length"] + ":" + df["insdel"] + ":" + df["feature_name"]
    )
    df.loc[:, "indel_type"] = df["feature"].replace(indel_feature_to_type)
    df.loc[:, "index"] = range(df.shape[0])
    feature_lengths = []
    for rix, row in df.iterrows():
        if row["insdel"] == "Del" and row["feature_name"] in ("C", "T"):
            feature_length = str(int(row["feature_length"]) + 1)
            if feature_length == "6":
                feature_length += "+"
        elif row["insdel"] == "Ins" and row["feature_name"] in ("C", "T"):
            feature_length = str(int(row["feature_length"]))
            if feature_length == "5":
                feature_length += "+"
        elif row["insdel"] == "Del" and row["feature_name"] in ("R",):
            feature_length = str(int(row["feature_length"]) + 1)
            if feature_length == "6":
                feature_length += "+"
        elif row["insdel"] == "Ins" and row["feature_name"] in ("R",):
            feature_length = str(int(row["feature_length"]))
            if feature_length == "5":
                feature_length += "+"
        elif row["insdel"] == "Del" and row["feature_name"] in ("M",):
            feature_length = str(int(row["feature_length"]))
            if feature_length == "5":
                feature_length += "+"
        else:
            print("ERROR: logic hole", file=sys.stderr)
            break
        feature_lengths.append(feature_length)
    df.loc[:, "feature_length"] = feature_lengths
    df = df[["count", "index", "indel_type", "feature_length"]]
    return df


def plot_indel_spectra(
    indel, title, tag="", save_path=None, yscale_log=False, ylim=None
):
    """Plot indel profile for given indel dataframe"""

    indel_colors = {
        "1bp Del at Homopolymer C": "#FDBE6E",
        "1bp Del at Homopolymer T": "#FD7F06",
        "1bp Ins at Homopolymer C": "#ACDC88",
        "1bp Ins at Homopolymer T": "#399F31",
        "2bp Del at Repeats": "#FACAB4",
        "3bp Del at Repeats": "#FC8A68",
        "4bp Del at Repeats": "#F14434",
        "5+bp Del at Repeats": "#BC191C",
        "2bp Ins at Repeats": "#D0E0F0",
        "3bp Ins at Repeats": "#94C3E1",
        "4bp Ins at Repeats": "#4C95C8",
        "5+bp Ins at Repeats": "#1962A7",
        "2bp Del at Microhomology": "#E2DFF0",
        "3bp Del at Microhomology": "#B4B7D8",
        "4bp Del at Microhomology": "#8584BD",
        "5+bp Del at Microhomology": "#614099",
    }

    df = proc_indel_dataframe(indel)

    fig, ax = plt.subplots(1)
    fig.set_figheight(4.5)
    fig.set_figwidth(20)

    if tag:
        title += " " + tag
    fig.suptitle(title)

    font = matplotlib.font_manager.FontProperties()
    font.set_family("monospace")

    for mut_type, mut_type_data in df.groupby(["indel_type"], sort=False):
        ax.bar(
            data=mut_type_data,
            x="index",
            height="count",
            label=mut_type,
            color=indel_colors[mut_type],
        )

    plt.xlabel("feature length")
    plt.xticks(df["index"], df["feature_length"].astype(str), fontproperties=font)
    plt.xlim((-1, 97))
    plt.legend(bbox_to_anchor=(0.86, 1), loc="upper left")

    if yscale_log:
        plt.yscale("log")
    if ylim:
        assert len(ylim) == 2
        plt.ylim(ylim)

    sns.despine(trim=True)
    if save_path:
        plt.savefig(save_path)
        plt.close()


def plot_venn2(cmp, weighted=False, label1="A", label2="B", title="", save_path=None):
    """Draw a venn diagram from a Snv/SvComparison instance"""
    if title:
        plt.title(title)
    if weighted:
        venn2([cmp.A, cmp.B], set_labels=(label1, label2))
    else:
        venn2_unweighted([cmp.A, cmp.B], set_labels=(label1, label2))
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
        plt.close()
