# Troubleshooting Guide

This guide helps you diagnose and fix common issues with SADIE's germlines module.

---

## Quick Diagnostics

Before diving into specific issues, run these diagnostic commands:

```bash
# Check germlines status
sadie germlines status

# Verify installation
python -c "import sadie; print(sadie.__version__)"

# Check disk space
df -h ~/.sadie/germlines

# Check environment variable
echo $SADIE_USE_GERMLINES_MODULE
```

---

## Setup Issues

### Issue 1: "Command not found: sadie"

**Symptom:**

```console
$ sadie germlines populate
-bash: sadie: command not found
```

**Root Cause:**

SADIE is not installed or not in your PATH.

**Solution:**

```bash
# Install SADIE
pip install sadie-antibody

# Verify installation
sadie --version

# If still not found, try:
python -m sadie.cli germlines populate
```

**Prevention:**

Ensure SADIE is installed in your active Python environment. Use virtual environments to manage dependencies:

```bash
# Create virtual environment
python -m venv sadie_env
source sadie_env/bin/activate  # Linux/Mac
# or: sadie_env\Scripts\activate  # Windows

# Install SADIE
pip install sadie-antibody
```

---

### Issue 2: "Permission denied" when creating directories

**Symptom:**

```console
$ sadie germlines populate
PermissionError: [Errno 13] Permission denied: '/home/user/.sadie/germlines'
```

**Root Cause:**

Insufficient permissions to write to `~/.sadie/germlines` directory.

**Solution:**

```bash
# Option A: Fix permissions
chmod -R u+w ~/.sadie/germlines

# Option B: Remove and recreate
rm -rf ~/.sadie/germlines
mkdir -p ~/.sadie/germlines
```

**For multi-user systems:**

```bash
# Set appropriate group permissions
sudo chgrp -R researchers ~/.sadie/germlines
sudo chmod -R g+w ~/.sadie/germlines
```

**Prevention:**

Ensure your user account has write access to the home directory. On shared systems, coordinate with system administrators for proper permissions.

---

## Download Issues

### Issue 3: Download hangs or times out

**Symptom:**

```console
$ sadie germlines populate
Processing imgt...
imgt: human ████░░░░░░░░░░░░░░░░ 20%
[hangs indefinitely]
```

**Root Cause:**

- Network connectivity issues
- Firewall blocking HTTPS connections
- Provider server temporarily unavailable
- Large file download taking time

**Solution:**

```bash
# Step 1: Test internet connectivity
ping -c 3 imgt.org

# Step 2: Check HTTPS access
curl -I https://www.imgt.org

# Step 3: Try single species first
sadie germlines populate -p imgt -s human

# Step 4: If successful, resume full download
sadie germlines populate
```

**Patience Approach:**

Some species legitimately take longer (5-10 minutes per species is normal for first download). Use `--dry-run` to preview:

```bash
sadie germlines populate --dry-run
```

**Prevention:**

- Download during off-peak hours
- Use stable internet connection
- Download specific species first, then expand
- Monitor with verbose logging if available

---

### Issue 4: "Species not found" after download

**Symptom:**

```console
$ sadie airr sequences.fasta -s rabbit -o output.tsv
Error: Species 'rabbit' not found in germline database
```

**Root Cause:**

Germline data for the species was not downloaded.

**Solution:**

```bash
# Check what species are downloaded
sadie germlines status

# Download the missing species
sadie germlines populate -p imgt -s rabbit

# Verify download
sadie germlines status

# Retry annotation
sadie airr sequences.fasta -s rabbit -o output.tsv
```

**For species verification:**

```bash
# List available species in IMGT
ls ~/.sadie/germlines/sources/imgt/

# Check if specific species has data
ls ~/.sadie/germlines/sources/imgt/rabbit/
```

**Prevention:**

Download all species you plan to use upfront:

```bash
# Download multiple species at once
sadie germlines populate -s human -s mouse -s rabbit -s rhesus_macaque
```

---

### Issue 5: Download fails partway through

**Symptom:**

```console
$ sadie germlines populate
Processing imgt...
imgt: mouse ████████████████████ 100%
Error downloading species 'rat': Connection reset
```

**Root Cause:**

Network interruption during download.

**Solution:**

```bash
# Simply re-run the command
# It will automatically resume from checkpoint
sadie germlines populate

# If checkpoint is corrupted, force restart
rm ~/.sadie/germlines/.populate_checkpoint.json
sadie germlines populate
```

**Prevention:**

The populate command automatically creates checkpoints. Just re-run after failures.

---

## Validation Issues

### Issue 6: "BLAST database build failed"

**Symptom:**

```console
$ sadie germlines populate
Running post-download build pipeline...
Building BLAST databases... ✗
Error: makeblastdb failed for IGHV
```

**Root Cause:**

- `makeblastdb` not installed or not in PATH
- Corrupted FASTA files
- Insufficient disk space

**Solution:**

```bash
# Step 1: Verify makeblastdb is installed
which makeblastdb

# If not found, install BLAST+
# On Ubuntu/Debian:
sudo apt-get install ncbi-blast+

# On macOS with Homebrew:
brew install blast

# On conda:
conda install -c bioconda blast

# Step 2: Verify disk space
df -h ~/.sadie/germlines

# Step 3: Rebuild databases
sadie germlines populate --force
```

**Prevention:**

Ensure BLAST+ is installed before running `sadie germlines populate`. BLAST+ is a required dependency for the germlines module.

---

### Issue 7: "Auxiliary file generation failed"

**Symptom:**

```console
Building auxiliary files... ✗
Error: Failed to create .aux files
```

**Root Cause:**

- Missing IMGT data files
- Incorrect FASTA format in custom sequences
- File permission issues

**Solution:**

```bash
# Step 1: Check FASTA format of custom sequences
head ~/.sadie/germlines/sources/custom/human/IGHV.fasta

# Should start with ">" and contain only ATCGN

# Step 2: Validate FASTA format
python -c "
from Bio import SeqIO
try:
    seqs = list(SeqIO.parse('~/.sadie/germlines/sources/custom/human/IGHV.fasta', 'fasta'))
    print(f'Valid FASTA with {len(seqs)} sequences')
except Exception as e:
    print(f'Invalid FASTA: {e}')
"

# Step 3: Remove problematic custom sequences temporarily
mv ~/.sadie/germlines/sources/custom ~/.sadie/germlines/sources/custom.backup

# Step 4: Rebuild
sadie germlines populate --force

# Step 5: Restore custom sequences after fixing
mv ~/.sadie/germlines/sources/custom.backup ~/.sadie/germlines/sources/custom
```

**Prevention:**

Validate custom FASTA files before adding them. See [Custom Sequences Guide](custom-sequences.md) for format requirements.

---

## Runtime Issues

### Issue 8: Annotation produces empty results

**Symptom:**

```python
from sadie.airr import Airr

airr = Airr("human")
results = airr.run_file("sequences.fasta")
# Results dataframe is empty
```

**Root Cause:**

- Germlines not populated
- Wrong species specified
- Input sequences are invalid
- Environment variable set to use G3 (deprecated)

**Solution:**

```bash
# Step 1: Verify germlines are populated
sadie germlines status

# Step 2: Check input sequences
head sequences.fasta

# Should be valid FASTA format

# Step 3: Check environment variable
echo $SADIE_USE_GERMLINES_MODULE
# Should be empty, "true", or "1"

# If set to "false", unset it:
unset SADIE_USE_GERMLINES_MODULE

# Step 4: Test with simple sequence
echo ">test" > test.fasta
echo "CAGGTGCAGCTGGTGCAGTCTGGGGCT" >> test.fasta
sadie airr test.fasta -o test_output.tsv

# Step 5: Check output
cat test_output.tsv
```

**Prevention:**

- Always run `sadie germlines populate` after installation
- Verify species spelling matches exactly
- Ensure `SADIE_USE_GERMLINES_MODULE` is not set to `false`

---

### Issue 9: "FileNotFoundError: No such file IGHV.ndb"

**Symptom:**

```console
$ sadie airr sequences.fasta -o output.tsv
FileNotFoundError: [Errno 2] No such file or directory:
'/Users/user/.sadie/germlines/databases/human/IGHV.ndb'
```

**Root Cause:**

BLAST databases not built or corrupted.

**Solution:**

```bash
# Rebuild databases for the species
sadie germlines populate -s human --force

# Verify database files exist
ls -lh ~/.sadie/germlines/databases/human/

# Should see: IGHV.ndb, IGHV.aux, IGHV.ndm.imgt, etc.
```

**Prevention:**

Don't manually delete files from `~/.sadie/germlines/databases/`. If you need to clean up, use:

```bash
# Safe cleanup: removes everything and rebuilds
rm -rf ~/.sadie/germlines
sadie germlines populate
```

---

### Issue 10: Results differ from G3 API

**Symptom:**

Annotation results from germlines module differ from legacy G3 API results.

**Root Cause:**

- Updated germline sequences in IMGT (monthly releases)
- Additional sequences from OGRDB
- Custom sequences taking priority

**Solution:**

This is **expected behavior**. Differences indicate:

1. **Improved accuracy**: Newer germline alleles provide better annotation
2. **Additional sources**: OGRDB adds validated novel alleles
3. **Custom priority**: Your custom sequences override provider sequences

**To compare:**

```bash
# Run with germlines (current default)
export SADIE_USE_GERMLINES_MODULE=true
sadie airr sequences.fasta -o results_germlines.tsv

# Run with G3 (deprecated, available until 2026-06-01)
export SADIE_USE_GERMLINES_MODULE=false
sadie airr sequences.fasta -o results_g3.tsv

# Compare
diff results_germlines.tsv results_g3.tsv
```

**Prevention:**

This is not an error. Expect differences when migrating to germlines module. Document the IMGT release version for reproducibility:

```bash
sadie germlines status
# Note the Version field (e.g., release-202601)
```

---

## Data Integrity Issues

### Issue 11: Disk space full during download

**Symptom:**

```console
OSError: [Errno 28] No space left on device
```

**Root Cause:**

Insufficient disk space for germline data (~500MB needed for all species).

**Solution:**

```bash
# Check available space
df -h ~/.sadie

# Free up space or use different location
# Download specific species only
sadie germlines populate -s human

# Human-only requires ~30MB instead of 500MB
```

**Prevention:**

Ensure at least 1GB free space before downloading all species.

---

### Issue 12: Corrupted database after system crash

**Symptom:**

```console
Error: Corrupted BLAST database
```

**Root Cause:**

System crash or power loss during database build.

**Solution:**

```bash
# Remove corrupted databases
rm -rf ~/.sadie/germlines/databases

# Rebuild from downloaded sources
sadie germlines populate --force

# If sources are also corrupted:
rm -rf ~/.sadie/germlines
sadie germlines populate
```

**Prevention:**

Backup your `~/.sadie/germlines` directory after successful setup:

```bash
tar -czf sadie_germlines_backup.tar.gz ~/.sadie/germlines
```

---

## Environment Issues

### Issue 13: G3 API still being used

**Symptom:**

Annotation is slow or requires internet, suggesting G3 API is still active.

**Root Cause:**

`SADIE_USE_GERMLINES_MODULE` environment variable explicitly set to `false`.

**Solution:**

```bash
# Check environment variable
echo $SADIE_USE_GERMLINES_MODULE

# If it shows "false", unset it:
unset SADIE_USE_GERMLINES_MODULE

# Or explicitly enable germlines:
export SADIE_USE_GERMLINES_MODULE=true

# Make permanent by adding to shell profile:
echo 'export SADIE_USE_GERMLINES_MODULE=true' >> ~/.bashrc
# or: ~/.zshrc for zsh users

# Reload shell configuration
source ~/.bashrc  # or ~/.zshrc
```

**Prevention:**

After migration, remove any `SADIE_USE_GERMLINES_MODULE=false` lines from shell profiles.

---

## Getting Help

### Collect Diagnostic Information

When reporting issues, include:

```bash
# System information
uname -a
python --version

# SADIE version
python -c "import sadie; print(sadie.__version__)"

# Germlines status
sadie germlines status

# Environment variable
echo $SADIE_USE_GERMLINES_MODULE

# Disk space
df -h ~/.sadie/germlines

# Database files
ls -lh ~/.sadie/germlines/databases/human/
```

### Report Issues

1. **GitHub Issues**: [https://github.com/jwillis0720/sadie/issues](https://github.com/jwillis0720/sadie/issues)
2. **Documentation Feedback**: Use the [documentation feedback template](https://github.com/jwillis0720/sadie/issues/new?template=documentation-feedback.md)


### Additional Resources

- [CLI Reference](cli-reference.md) - Command documentation
- [Migration Guide](migration-guide.md) - G3 to germlines migration
- [Provider Guide](provider-guide.md) - Understanding providers
- [Custom Sequences](custom-sequences.md) - Adding custom germlines
