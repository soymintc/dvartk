import pandas as pd
import numpy as np
from pyfaidx import Fasta

genome = Fasta("/juno/work/shah/users/chois7/mmctm/reference/GRCh37-lite.fa")


def construct_empty_count_series():
    snvs = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    nts = ["A", "C", "G", "T"]
    terms = ["{}[{}]{}".format(l, s, r) for s in snvs for l in nts for r in nts]
    return pd.Series(np.zeros(len(terms), dtype=int), index=terms)


def normalize_snv(context, alt):
    ref = context.seq[1]

    if ref in ["A", "G"]:
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


def count_snvs(snvs, genome=genome):
    """Convert maf form to count table per variant type. Requires 'genome'"""
    df = snvs.copy()

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

    return counts
