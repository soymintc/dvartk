import uuid
import pandas as pd
import numpy as np
from pyfaidx import Fasta
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

genome = Fasta("/juno/work/shah/users/chois7/mmctm/reference/GRCh37-lite.fa")


def construct_empty_count_series():
    snvs = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    nts = ["A", "C", "G", "T"]
    terms = ["{}[{}]{}".format(l, s, r) for s in snvs for l in nts for r in nts]
    return pd.Series(np.zeros(len(terms), dtype=int), index=terms)


def normalize_snv(context, alt):
    ref = context.seq[1]

    if ref in ("A", "G"):
        context = (-context).seq

        complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
        alt = complement[str(alt)]
    else:
        context = context.seq
        alt = str(alt)
    return context, alt


def construct_snv_label(context, alt):
    if len(context) != 3:
        warnings.warn("Warning: bad context length: {}".format(str(context)))
        return None
    return "{}[{}>{}]{}".format(context[0], context[1], alt, context[2])


def count_svs(maf):
    """Convert maf to count table as according to palimpsest"""

    svtypes = ["del", "dup", "ins", "inv", "translocation"]

    svlen_bins = [0, 1e3, 1e4, 1e5, 1e6, 1e7, np.inf]
    svlen_bin_labels = ["<1kb", "1-10kb", "10-100kb", "100kb-1Mb", "1-10Mb", ">10Mb"]

    sv_labels = []
    for svtype in svtypes:
        if not svtype.startswith("tr"):
            for svlen_bin_label in svlen_bin_labels:
                sv_label = f"{svtype}:{svlen_bin_label}"
                sv_labels.append(sv_label)
        else:
            sv_labels.append(svtype)

    df = maf.copy()
    df.loc[df["type"] == "translocation", "length"] = 3e9
    df["sv_length"] = pd.cut(
        df["length"], bins=svlen_bins, labels=svlen_bin_labels
    ).values
    df["sv_category"] = df["type"].str.cat(df["sv_length"], sep=":")
    df.loc[df["type"] == "translocation", "sv_category"] = "translocation"

    sv_counts = pd.Series([0] * len(sv_labels), index=sv_labels)
    sv_count_values = sv_counts.index.map(df["sv_category"].value_counts()).fillna(0)
    sv_counts = pd.Series(sv_count_values, index=sv_labels).astype(int)
    return sv_counts


def count_indels(df, genome_version="GRCh37"):
    """df: pandas DataFrame of chrom, pos, ref, alt columns
    - chrom [str]: chromosome ID, e.g. 'chr1', '1'
    - pos [int]: 1-based VCF format indel coordinate
    - ref: VCF format indel reference
    - alt: VCF format indel alteration
    genome_version [str]: element in {'GRCh37', 'GRCh38'}
    """
    tmp_dirname = f"_{str(uuid.uuid4())}"
    if not os.path.exists(tmp_dirname):
        subprocess.run(["mkdir", tmp_dirname])
    df["ID"] = "."
    df = df[["chrom", "pos", "ID", "ref", "alt"]]
    df.columns = ["#CHROM", "POS", "ID", "REF", "ALT"]
    tmp_vcf_path = f"{tmp_dirname}/indels.vcf"
    df.to_csv(tmp_vcf_path, sep="\t", index=False)
    project = "indels"
    matrices = matGen.SigProfilerMatrixGeneratorFunc(
        project, genome_version, tmp_dirname
    )
    counts = matrices["ID"]
    counts.columns = ["count"]
    if os.path.exists(tmp_dirname):
        subprocess.run(["rm", "-rf", tmp_dirname])
    return counts


def count_snvs(snvs, genome=genome):
    """Convert maf form to count table per variant type. Requires 'genome'"""
    df = snvs.copy()  # snvs <- essentially "maf" variable

    var_converter = {
        "G>T": "C>A",
        "G>C": "C>G",
        "G>A": "C>T",
        "A>T": "T>A",
        "A>G": "T>C",
        "A>C": "T>G",
    }
    df["var_type"] = (df["ref"] + ">" + df["alt"]).replace(var_converter)

    counts = construct_empty_count_series()
    contexts = []

    for idx, row in snvs.iterrows():
        # two flanking bases
        start = row["pos"] - 2
        end = row["pos"] + 1

        context = genome[row["chrom"]][start:end]
        if "N" in context.seq:
            warnings.warn(
                "Warning: N found in context sequence at {}:{}-{}".format(
                    row["chrom"], start + 1, end
                )
            )
            continue

        context, alt = normalize_snv(context, row["alt"])
        counts[construct_snv_label(context, alt)] += 1
        contexts.append(context)

    assert snvs.shape[0] == counts.sum()

    return counts
