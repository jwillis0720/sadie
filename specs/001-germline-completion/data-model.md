# Data Model: Germlines Module

## Entities

### GermlineGene

Core entity representing a single germline gene sequence.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| name | str | Y | Gene name (IGHV1-69*01) |
| species | str | Y | Species (human, mouse) |
| segment | str | Y | V, D, or J |
| chain | str | Y | H, K, or L |
| sequence | str | Y | Ungapped nucleotide |
| sequence_gapped | str | N | IMGT-gapped nucleotide |
| sequence_aa | str | N | Ungapped amino acid |
| sequence_aa_gapped | str | N | IMGT-gapped amino acid |
| is_functional | bool | Y | Functional status |
| functionality | str | Y | F, ORF, or P |
| regions | dict | N | CDR/FWR regions |
| region_positions | dict | N | Region boundaries |
| source | str | Y | imgt, ogrdb, vdjbase, custom |
| source_version | str | N | Version/date |
| allele | str | N | Allele designation |
| gene_family | str | N | Gene family |
| accession | str | N | GenBank accession |

### ProviderMetadata

Tracks provider information.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| name | str | Y | Provider name |
| version | str | Y | Version identifier |
| last_updated | datetime | Y | Update timestamp |
| species_available | List[str] | Y | Available species |
| url | str | N | Source URL |

### ProcessingMetadata

Tracks processed file state.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| source_file | str | Y | Source path |
| processed_at | datetime | Y | Processing time |
| num_sequences | int | Y | Sequence count |
| file_hash | str | Y | MD5 hash |
| sequences | List[dict] | Y | Sequence summaries |

## Relationships

```
Provider (1) --- (*) GermlineGene
GermlineGene --- (0..1) ProcessingMetadata
Provider --- (1) ProviderMetadata
```

## Directory Structure

```
germlines/
├── sources/           # Raw data by provider
│   ├── imgt/human/   # IMGT FASTA files
│   ├── ogrdb/human/  # OGRDB FASTA files
│   ├── vdjbase/human/# VDJbase FASTA files
│   └── custom/human/ # User custom files
├── normalized/        # Merged output
│   └── human/
│       ├── gapped/   # IMGT-gapped FASTA
│       └── ungapped/ # Plain FASTA
└── igblast/          # BLAST databases
    ├── blastdb/      # makeblastdb output
    ├── aux_db/       # CDR/FWR boundaries
    └── internal_data/# organism.yaml
```

## Provider Priority

Default: `["custom", "ogrdb", "vdjbase", "imgt"]`

Resolution rules:
1. Same name, different sequence: Use higher priority
2. Same name, same sequence: Keep one, track source
3. Novel gene: Include from any source

## File Naming

Sources: `{SEGMENT}.fasta`, `{SEGMENT}_gapped.fasta`
- IGHV.fasta, IGHV_gapped.fasta
- IGHD.fasta, IGHJ.fasta

Normalized: `{species}_{segment}.fasta`
- human_V.fasta, human_V_gapped.fasta

BLAST: `{species}_{segment}`
- human_V.nhr, human_V.nin, human_V.nsq

## Validation Rules

1. Segment: Must be V, D, or J
2. Chain: Must be H, K, or L
3. Sequence: Valid nucleotides (ACGTN) and IUPAC ambiguous
4. Functionality: F (functional), ORF, or P (pseudogene)
5. Gapped sequences: Use dots (.) for gaps per IMGT
