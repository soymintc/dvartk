# dvartk (variant comparison toolkit)
## INSTALL
```
pip install dvartk
```

## Philosophy
Basically a package for comparing SNVs of maf files
- in the future, SV, INDEL support, vcf support as well

## Usage
### Compariing SNVs
```
import dvartk

# set a MAF paths
maf1_path = '/path/to/maf1'
maf2_path = '/path/to/maf2'

# convert custom chrom, pos, ref, alt column names to "chrom", "pos", "ref", "alt"
snv_file_config = dvartk.parser.SnvFileConfig('Chromosome', 'Start_Position', 'Reference_Allele', 'Tumor_Seq_Allele2')

# load and convert columns of MAF
maf1 = dvartk.load_and_convert_snv_maf_columns(maf1_path, snv_file_config)
maf2 = dvartk.load_and_convert_snv_maf_columns(maf2_path, snv_file_config)

# get set counts (intersection, difference, ...) between maf1 and maf2
# first make a SNV comparison instance
snv_cmp = SnvComparison(maf1, maf2)
snv_cmp.get_set_counts() # get A[maf1], B[maf2] counts
summary = snv_cmp.make_oneliner()
print(summary)
```

### Plot SNV trinucleotide spectra
```
import dvartk

# set a MAF paths
maf_path = '/path/to/maf'

# get counts per trinucleotide type
counts = dvartk.count_snvs(maf)
dvartk.plot_snv_spectra(
    counts,
    'plot title',  # title of the plot e.g. sample ID
    save_path=None,  # output path for your plot; if None, don't save plot
    tag=' (bar)'  # tag to add to plot title; if None, don't add tag
)
```
