# Quickstart: Germlines Module

## Basic Usage

```python
from sadie.germlines import get_germline_genes, GermlineManager

genes = get_germline_genes("human", "V", "H")

for gene in genes[:5]:
    print(f"{gene.name}: {gene.sequence[:30]}...")
```

## Setup (First Time)

```bash
cd src/sadie/germlines

python scripts/download_imgt.py --species human
python scripts/download_ogrdb.py --species human

python scripts/validate.py human
```

## Add Custom Sequences

```bash
mkdir -p sources/custom/human

cat > sources/custom/human/IGHV.fasta << 'EOF'
>IGHV-CUSTOM*01|Homo sapiens|F
CAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAG
EOF
```

Custom sequences take priority over database sources.

## Priority Configuration

```python
manager = GermlineManager(providers=["custom", "ogrdb", "imgt"])

genes = manager.get_genes("human", "V", "H")
```

Default: `["custom", "ogrdb", "vdjbase", "imgt"]`

## Offline Usage

Once data is downloaded, the module works offline:

```python
import socket
socket.setdefaulttimeout(0.1)

genes = get_germline_genes("human", "V", "H")
```

## Integration Points

### IgBLAST

```python
from sadie.airr.igblast.germline import GermlineData

gd = GermlineData("human")
print(gd.v_gene_dir)
```

### Get Specific Gene

```python
from sadie.germlines import get_gene_by_name

gene = get_gene_by_name("IGHV1-69*01", "human")
print(gene.sequence_gapped)
```

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| SADIE_USE_GERMLINES_MODULE | true | Use new module |

## Troubleshooting

**No data found**: Run download scripts first
**Import error**: Check BioPython installed
**Permission denied**: Check sources/ directory permissions
