# dvartk (variant comparison toolkit)
## Install
```
pip install dvartk
```

## Summary
Basically a package for comparing SNVs and SVs of maf files.
- in the future, INDEL support, input VCF support as well.
For whatever variant table file you have, as long as you designate the source column names in `*_src`, `dvartk` will process your results.

## Usage
### Comparing SNVs
```python
import dvartk

# set MAF paths
maf1_path = '/path/to/maf1'
maf2_path = '/path/to/maf2'

# convert custom chrom, pos, ref, alt column names to "chrom", "pos", "ref", "alt"
snv_file_config = dvartk.parser.SnvFileConfig(
    'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2')

# load and convert columns of MAF
maf1 = snv_file_config.load_and_convert_maf_columns(maf1_path)
maf2 = snv_file_config.load_and_convert_maf_columns(maf2_path)

# get set counts (intersection, difference, ...) between maf1 and maf2
# first make a SNV comparison instance
snv_cmp = SnvComparison(maf1, maf2)
snv_cmp.get_set_counts(get_return=False) # get A[maf1], B[maf2] counts
summary = snv_cmp.make_oneliner()
print(summary) # returns [#(A), #(B), #(A-B), #(B-A), #(A&B), #(A|B)]
```

### Plot SNV trinucleotide spectra
```python
import dvartk

# extract MAF data
maf_path = '/path/to/maf'
snv_file_config = dvartk.parser.SnvFileConfig(
    'Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2')
maf = dvartk.load_and_convert_maf_columns(maf_path, snv_file_config)

# get counts per trinucleotide type
counts = dvartk.count_snvs(maf)
dvartk.plot_snv_spectra(
    counts=counts,
    title='foo',     # title of the plot e.g. sample ID
    save_path=None,  # output path for your plot; if None, don't save plot
    tag=' (bar)'     # tag to add to plot title; if None, don't add tag
    yscale_log=True, # make yscale log
    ylim=None,       # set ylim for plot
    debug=False,     # print inside variables
)
```

### Comparing SVs
```python
import dvartk

# set MAF paths
maf1_path = '/path/to/maf1'
maf2_path = '/path/to/maf2'

# convert custom chrom, pos, ref, alt column names to:
# "chromosome_1", "position_1", "strand_1", "chrom2", "pos2", "strand2", "type", "length"
sv_file_config = dvartk.parser.SvFileConfig(
    'CHROM_A', 'Start_Position1', 'Strand1',
    'CHROM_B', 'Start_Position2', 'Strand2', 'SV_Type', 'Breakpoint_Length')

# load and convert columns of MAF; can have different configs too
maf1 = sv_file_config.load_and_convert_maf_columns(maf1_path)
maf2 = sv_file_config.load_and_convert_maf_columns(maf2_path)

# get set counts (intersection, difference, ...) between maf1 and maf2
# first make a SV comparison instance
snv_cmp = SvComparison(maf1, maf2)
snv_cmp.get_set_counts(get_return=False) # get A[maf1], B[maf2] counts
summary = snv_cmp.make_oneliner()
print(summary) # returns [#(A), #(B), #(A-B), #(B-A), #(A&B), #(A|B)]
```

### Plot SV palimpsest-like spectra
```python
import dvartk

# extract MAF data
maf_path = '/path/to/maf'
sv_file_config = dvartk.parser.SvFileConfig(
    'CHROM_A', 'Start_Position1', 'Strand1',
    'CHROM_B', 'Start_Position2', 'Strand2', 'SV_Type', 'Breakpoint_Length')
maf = sv_file_config.load_and_convert_maf_columns(maf_path)

# get counts per trinucleotide type
counts = dvartk.count_svs(maf)
dvartk.plot_sv_spectra(
    counts=counts,
    title='foo',     # title of the plot e.g. sample ID
    save_path=None,  # output path for your plot; if None, don't save plot
    tag=' (bar)'     # tag to add to plot title; if None, don't add tag
    yscale_log=True, # make yscale log
    ylim=None,       # set ylim for plot
    debug=False,     # print inside variables
)
```
