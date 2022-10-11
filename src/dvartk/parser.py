import gzip
import pandas as pd
import wgs_analysis.algorithms.rearrangement


def convert_type_names(maf, type_col_name="type"):
    """Convert SV type names"""
    df = maf.copy()
    type_src1 = ["DEL", "INS", "DUP", "INV", "BND"]
    type_src2 = ["Deletion", "Insertion", "Duplication", "Inversion", "Translocation"]
    type_src3 = ["deletion", "insertion", "duplication", "inversion", "translocation"]
    type_dst = ["del", "ins", "dup", "inv", "translocation"]

    df[type_col_name] = df[type_col_name].replace(dict(zip(type_src1, type_dst)))
    df[type_col_name] = df[type_col_name].replace(dict(zip(type_src2, type_dst)))
    df[type_col_name] = df[type_col_name].replace(dict(zip(type_src3, type_dst)))
    return df


class SvFileConfig:
    """Config with source and dest column names"""

    def __init__(
        self,
        chromosome_1_src,
        position_1_src,
        strand_1_src,
        chromosome_2_src,
        position_2_src,
        strand_2_src,
        type_src,
        length_src,
    ):
        self.chromosome_1_src = chromosome_1_src
        self.position_1_src = position_1_src
        self.strand_1_src = strand_1_src
        self.chromosome_2_src = chromosome_2_src
        self.position_2_src = position_2_src
        self.strand_2_src = strand_2_src
        self.type_src = type_src
        self.length_src = length_src
        self.chromosome_1_dst = "chromosome_1"
        self.position_1_dst = "position_1"
        self.strand_1_dst = "strand_1"
        self.chromosome_2_dst = "chromosome_2"
        self.position_2_dst = "position_2"
        self.strand_2_dst = "strand_2"
        self.type_dst = "type"
        self.length_dst = "length"

        self.col_converter = {
            self.chromosome_1_src: self.chromosome_1_dst,
            self.position_1_src: self.position_1_dst,
            self.strand_1_src: self.strand_1_dst,
            self.chromosome_2_src: self.chromosome_2_dst,
            self.position_2_src: self.position_2_dst,
            self.strand_2_src: self.strand_2_dst,
            self.type_src: self.type_dst,
            self.length_src: self.length_dst,
        }

    def load_maf(self, maf_path):
        """Load a maf"""
        if maf_path.endswith("gz"):
            maf_file = gzip.open(maf_path, "rb")
            head = maf_file.read(10000).decode("utf-8")
        else:
            maf_file = open(maf_path, "r")
            head = maf_file.read(10000)
        delimitor = "\t" if head.count("\t") > 0 else ","
        maf = pd.read_csv(
            maf_path,
            dtype={
                self.chromosome_1_src: str,
                self.chromosome_2_src: str,
            },
            sep=delimitor,
            comment="#",
            low_memory=False,
        )
        return maf

    def convert_maf_columns(self, maf):
        """Convert column names"""
        assert len(set(self.col_converter.keys()) & set(maf.columns)) == len(
            self.col_converter
        )  # all src names should be in maf columns
        maf = maf.rename(columns=self.col_converter)
        maf = convert_type_names(maf)
        return maf

    def load_and_convert_maf_columns(self, maf_path):
        """Load a maf, then convert column names"""
        maf = self.load_maf(maf_path)
        return self.convert_maf_columns(maf)


class SnvFileConfig:
    """Config with source and dest column names"""

    def __init__(self, chrom_src, pos_src, ref_src, alt_src):
        self.chrom_src = chrom_src
        self.pos_src = pos_src
        self.ref_src = ref_src
        self.alt_src = alt_src
        self.chrom_dst = "chrom"
        self.pos_dst = "pos"
        self.ref_dst = "ref"
        self.alt_dst = "alt"

        self.col_converter = {
            self.chrom_src: self.chrom_dst,
            self.pos_src: self.pos_dst,
            self.ref_src: self.ref_dst,
            self.alt_src: self.alt_dst,
        }

    def load_maf(self, maf_path):
        """Load a maf"""
        if maf_path.endswith("gz"):
            maf_file = gzip.open(maf_path, "rb")
            head = maf_file.read(10000).decode("utf-8")
        else:
            maf_file = open(maf_path, "r")
            head = maf_file.read(10000)
        delimitor = "\t" if head.count("\t") > 0 else ","
        maf = pd.read_csv(
            maf_path,
            dtype={
                self.chrom_src: str,
            },
            sep=delimitor,
            comment="#",
            low_memory=False,
        )
        return maf

    def select_SNPs(self, maf):
        """Select variants with SNP tags"""
        assert "Variant_Type" in maf.columns  # TODO: generalize
        maf = maf[maf["Variant_Type"] == "SNP"]  # TODO: generalize
        return maf

    def convert_maf_columns(self, maf):
        """Convert column names"""
        assert len(set(self.col_converter.keys()) & set(maf.columns)) == len(
            self.col_converter
        )  # all src names should be in maf columns
        maf = maf.rename(columns=self.col_converter)
        return maf

    def load_and_convert_maf_columns(self, maf_path):
        """Load a maf, select SNPs only, then convert column names"""
        maf = self.load_maf(maf_path)
        maf = self.select_SNPs(maf)
        return self.convert_maf_columns(maf)


class SvComparison:
    """Class for comparing two SV 'maf' tables"""

    ixs = [
        "chromosome_1",
        "position_1",
        "strand_1",
        "chromosome_2",
        "position_2",
        "strand_2",
        "type",
    ]

    def __init__(self, maf1, maf2, delimitor="\t", debug=False):
        self.maf1 = maf1.reset_index(drop=True)
        self.maf2 = maf2.reset_index(drop=True)
        self.delimitor = delimitor
        self.debug = debug

        if self.debug:
            print(f"self.maf1: {self.maf1}")
            print(f"self.maf2: {self.maf2}")

        maf1["prediction_id"] = maf1.index
        maf2["prediction_id"] = maf2.index

        try:
            sv_match = wgs_analysis.algorithms.rearrangement.match_breakpoints(
                maf1, maf2, window_size=200
            )
        except ValueError:
            sv_match = pd.DataFrame(columns=["reference_id", "target_id"])

        self.maf1_match = self.maf1[
            self.maf1.prediction_id.isin(sv_match["reference_id"])
        ]
        self.maf1_nonmatch = self.maf1[
            ~self.maf1.prediction_id.isin(sv_match["reference_id"])
        ]
        self.maf2_match = self.maf2[self.maf2.prediction_id.isin(sv_match["target_id"])]
        self.maf2_nonmatch = self.maf2[
            ~self.maf2.prediction_id.isin(sv_match["target_id"])
        ]

    def make_set(self, data):
        df_ix = data.copy().set_index(self.ixs).index.tolist()
        df_ix_set = set(df_ix)
        return df_ix_set

    def get_set_counts(self, get_return=False):
        self.A_and_B = self.make_set(
            self.maf1_match
        )  # big assumption that results will be same/similar with self.maf2_match
        self.A_and_B_from_maf2 = self.make_set(
            self.maf2_match
        )  # for twisted people like Seongmin who can't trust himself
        self.A_not_B = self.make_set(self.maf1_nonmatch)
        self.B_not_A = self.make_set(self.maf2_nonmatch)
        self.A = self.A_and_B.union(self.A_not_B)
        self.B = self.A_and_B.union(self.B_not_A)

        if get_return:
            return (
                self.A,
                self.B,
                self.A_not_B,
                self.B_not_A,
                self.A_and_B,
                self.A | self.B,
            )

    def make_oneliner(
        self, name=None, get_str=False, print_header=False, delimitor="\t"
    ):
        """Returns #A, #B, #(A-B), #(B-A), #(A&B), #(A|B)"""
        if print_header:
            print("A B A-B B-A A&B A|B".replace(" ", delimitor))
        for attr in ["A", "B"]:
            if not hasattr(self, attr):
                self.get_set_counts()
        field = [
            len(self.A),
            len(self.B),
            len(self.A_not_B),
            len(self.B_not_A),
            len(self.A_and_B),
            len(self.A | self.B),
        ]
        if name:
            field = [name] + field
        field = [str(_) for _ in field]
        if get_str:
            return self.delimitor.join(field)
        return field


class SnvComparison:
    """Class for comparing two 'maf' tables"""

    ixs = ["chrom", "pos", "ref", "alt"]

    def __init__(self, maf1, maf2, delimitor="\t", debug=False):
        self.maf1 = maf1.reset_index(drop=True)
        self.maf2 = maf2.reset_index(drop=True)
        self.delimitor = delimitor
        self.debug = debug

        if self.debug:
            print(f"self.maf1: {self.maf1}")
            print(f"self.maf2: {self.maf2}")

        snv_match = pd.merge(
            self.maf1[self.ixs], self.maf2[self.ixs], how="inner", on=self.ixs
        ).set_index(self.ixs)

        self.maf1 = self.maf1.set_index(self.ixs, drop=False)
        self.maf2 = self.maf2.set_index(self.ixs, drop=False)

        self.maf1_match = self.maf1[
            self.maf1.prediction_id.isin(snv_match.prediction_id)
        ]
        self.maf1_nonmatch = self.maf1[
            ~self.maf1.prediction_id.isin(snv_match.prediction_id)
        ]
        self.maf2_match = self.maf2[
            self.maf2.prediction_id.isin(snv_match.prediction_id)
        ]
        self.maf2_nonmatch = self.maf2[
            ~self.maf2.prediction_id.isin(snv_match.prediction_id)
        ]

    def make_set(self, data):
        df_ix = data.copy().set_index(self.ixs).index.tolist()
        df_ix_set = set(df_ix)
        return df_ix_set

    def get_set_counts(self, get_return=False):
        self.A_and_B = self.make_set(
            self.maf1_match
        )  # big assumption that results will be same/similar with self.maf2_match
        self.A_and_B_from_maf2 = self.make_set(
            self.maf2_match
        )  # for twisted people like Seongmin who can't trust himself
        self.A_not_B = self.make_set(self.maf1_nonmatch)
        self.B_not_A = self.make_set(self.maf2_nonmatch)
        self.A = self.A_and_B.union(self.A_not_B)
        self.B = self.A_and_B.union(self.B_not_A)

        if get_return:
            return (
                self.A,
                self.B,
                self.A_not_B,
                self.B_not_A,
                self.A_and_B,
                self.A | self.B,
            )

    def make_oneliner(self, name=None, get_str=False):
        for attr in ["A", "B"]:
            if not hasattr(self, attr):
                self.get_set_counts()
        field = [
            len(self.A),
            len(self.B),
            len(self.A - self.B),
            len(self.B - self.A),
            len(self.A & self.B),
            len(self.A | self.B),
        ]
        if name:
            field = [name] + field
        field = [str(_) for _ in field]
        if get_str:
            return self.delimitor.join(field)
        return field
