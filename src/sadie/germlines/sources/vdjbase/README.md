# VDJbase Data Source

Population-specific germline alleles from VDJbase.

## About VDJbase

VDJbase provides germline sequences discovered through repertoire sequencing and genotype inference from multiple population studies.

**Website**: https://vdjbase.org/
**API Base**: https://vdjbase.org/admin/api/

## Supported Species

| Species | Internal Name | Datasets |
|---------|---------------|----------|
| Human | `human` | IGH, IGK, IGL, IGHC |
| Rhesus Macaque | `rhesus_macaque` | IGH, IGK, IGL |

**Note**: VDJbase currently only offers Human and Rhesus Macaque data. Mouse is not available (check `provider.get_available_species()` for updates).

## Automatic Download

Use the VDJbase provider to download data programmatically:

```python
from sadie.germlines.providers import VDJbaseProvider

provider = VDJbaseProvider()

# Download all available species
provider.download(["human", "rhesus_macaque"])

# Check what's available
print(provider.get_available_species())
print(provider.get_available_datasets("human"))
```

Or use the command line:

```bash
python -c "
from sadie.germlines.providers import VDJbaseProvider
provider = VDJbaseProvider()
provider.download(['human', 'rhesus_macaque'])
"
```

## Manual Download

If automatic download fails, you can manually download from the VDJbase website:

1. Visit https://vdjbase.org/
2. Navigate to the species/chain of interest
3. Export sequences as FASTA
4. Save to the appropriate directory:
   - `sources/vdjbase/human/IGHV.fasta`
   - `sources/vdjbase/human/IGHD.fasta`
   - `sources/vdjbase/human/IGHJ.fasta`
   - etc.

## Directory Structure

```
sources/vdjbase/
├── human/
│   ├── IGHV.fasta    # Heavy chain V genes
│   ├── IGHD.fasta    # Heavy chain D genes
│   ├── IGHJ.fasta    # Heavy chain J genes
│   ├── IGKV.fasta    # Kappa V genes
│   ├── IGKJ.fasta    # Kappa J genes
│   ├── IGLV.fasta    # Lambda V genes
│   └── IGLJ.fasta    # Lambda J genes
├── rhesus_macaque/
│   └── ...
└── README.md
```

## FASTA Header Format

VDJbase FASTA headers follow this format:

```
>{gene_name}|{optional_metadata}
```

Example:
```
>IGHV1-69*01|novel=true|appears=42
CAGGTGCAGCTGGTGGAGTCTGGG...
```

Optional metadata fields:
- `novel=true`: Novel allele not in IMGT/OGRDB
- `low_confidence=true`: Low confidence call
- `appears=N`: Number of samples where allele appears

## Important Notes

1. **Ungapped Sequences**: VDJbase sequences are ungapped. The gapper module will add IMGT gaps automatically.

2. **Novel Alleles**: Many VDJbase sequences are novel alleles discovered through repertoire sequencing. These may not exist in IMGT or OGRDB.

3. **Priority**: By default, VDJbase has lower priority than OGRDB and IMGT. Configure priority in GermlineManager if you want VDJbase sequences to take precedence.

4. **Updates**: VDJbase is actively updated. Re-run download periodically to get new sequences.

## API Endpoints

| Endpoint | Purpose |
|----------|---------|
| `/genomic/species` | List available species |
| `/genomic/data_sets/{species}` | List datasets for species |
| `/repseq/sequences/{species}/{chain}` | Get sequences (paginated) |

## Troubleshooting

**No sequences downloaded?**
- Check network connectivity
- Verify species name is correct (use `get_available_species()`)
- Check VDJbase website for service status

**Import errors?**
- Ensure BioPython is installed: `pip install biopython`

## Citation

If you use VDJbase data, please cite:
> VDJbase: A Population-Specific Germline Repository
> https://vdjbase.org/
