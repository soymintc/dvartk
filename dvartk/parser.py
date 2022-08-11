import pandas as pd


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


def load_and_convert_snv_maf_columns(maf_path, snv_file_config):
    """Load a maf, select SNPs only, then convert column names"""
    maf = pd.read_csv(
        maf_path, dtype={"Chromosome": str}, sep="\t", comment="#", low_memory=False
    )
    assert "Variant_Type" in maf.columns  # TODO: generalize
    maf = maf[maf["Variant_Type"] == "SNP"]  # TODO: generalize

    assert len(set(snv_file_config.col_converter.keys()) & set(maf.columns)) == len(
        snv_file_config.col_converter
    )  # all src names should be in maf columns
    maf = maf.rename(columns=snv_file_config.col_converter)
    return maf


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

        self.maf1_match = self.maf1[self.maf1.index.isin(snv_match.index)]
        self.maf1_nonmatch = self.maf1[~self.maf1.index.isin(snv_match.index)]
        self.maf2_match = self.maf2[self.maf2.index.isin(snv_match.index)]
        self.maf2_nonmatch = self.maf2[~self.maf2.index.isin(snv_match.index)]

    def make_set(self, data):
        df_ix = data.copy().set_index(self.ixs).index.tolist()
        df_ix_set = set(df_ix)
        return df_ix_set

    def get_set_counts(self):
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

    def make_oneliner(self, name=None, get_str=False):
        for attr in ["A", "B"]:
            assert hasattr(self, attr)
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
