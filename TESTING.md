# Testing ProkaRegia Before Bioconda Submission

This guide walks you through testing the ProkaRegia pipeline locally before submitting to bioconda.

## Prerequisites

1. **Install conda-build**:
   ```bash
   conda install conda-build
   ```

2. **Get test data**: You'll need a small FASTQ file for testing. You can:
   - Use a subset of real data: `head -40000 your_reads.fastq > test_reads.fastq`
   - Or download a small test dataset

## Test 0: Local Conda Build (RECOMMENDED - Do This First!)

This simulates exactly what bioconda will do when building your package:

```bash
cd prokaregia_snakemake

# Build the conda package locally
conda build conda_recipe_example/

# This will:
# 1. Create a temporary build environment
# 2. Download your source files (or use local files)
# 3. Run the build script
# 4. Install to $PREFIX/share/prokaregia and $PREFIX/bin/prokaregia
# 5. Run the test command (prokaregia -h)
# 6. Create a .tar.bz2 package file

# If successful, you'll see output like:
# TEST START: prokaregia-1.0.0-0
# ...
# TEST END: prokaregia-1.0.0-0
```

**Install locally from the built package**:
```bash
# After build succeeds, install from local build
conda create -n test-prokaregia
conda activate test-prokaregia
conda install --use-local prokaregia

# Now test it as a user would
prokaregia -h
prokaregia -i test_reads.fastq -s ont -t 2

# Clean up when done
conda deactivate
conda env remove -n test-prokaregia
```

**Common build issues**:
- If `source.url` fails (because repo isn't public yet), use local path instead:
  ```yaml
  source:
    path: ..
  ```
- If dependencies aren't found, check they exist in conda-forge or bioconda
- If test fails, check the wrapper script permissions

## Test 1: Wrapper Script Validation

Test that the wrapper script properly validates inputs:

```bash
cd prokaregia_snakemake

# Test help
./run_prokaregia.sh -h

# Test missing arguments (should error)
./run_prokaregia.sh -i test.fastq
./run_prokaregia.sh -s ont

# Test invalid seq_tech (should error)
./run_prokaregia.sh -i test_reads.fastq -s illumina

# Test non-existent file (should error)
./run_prokaregia.sh -i nonexistent.fastq -s ont
```

**Expected**: All should show appropriate error messages.

## Test 2: Snakemake Dry Run

Test that the wrapper correctly builds the snakemake command:

```bash
# Create a minimal test file
echo -e "@read1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII" > test_reads.fastq

# Run dry run (checks workflow without executing)
snakemake --snakefile Snakefile --use-conda --cores 2 \
  --config input_fastq=test_reads.fastq seq_tech=ont output_dir=TestOutput -n
```

**Expected**: Should show the DAG of jobs without errors. You should see:
- Download CheckM2 database (if not already present)
- All pipeline steps listed
- No syntax errors

## Test 3: CheckM2 Database Download

Test that the database download rule works:

```bash
# Force download of CheckM2 database
snakemake --snakefile Snakefile --use-conda --cores 2 \
  --config input_fastq=test_reads.fastq seq_tech=ont \
  ~/.checkm2/CheckM2_database/uniref100.KO.1.dmnd
```

**Expected**:
- Should download ~3.5GB to `~/.checkm2/`
- Creates the complete CheckM2 database structure
- Subsequent runs should detect it exists and skip download

## Test 4: Test with Wrapper Script

Test the full wrapper functionality:

```bash
# Test with minimal arguments
./run_prokaregia.sh -i test_reads.fastq -s ont -t 2

# Test with custom output
./run_prokaregia.sh -i test_reads.fastq -s ont -t 2 -o MyTestOutput

# Test with mamba flag
./run_prokaregia.sh -i test_reads.fastq -s ont -t 2 -m
```

**Expected**:
- Wrapper shows configuration summary
- Executes snakemake command
- Pipeline starts running (may fail on test data, that's okay)

## Test 5: Verify File Paths

Verify that paths are correctly resolved:

```bash
# Create wrapper as it would be installed by conda
mkdir -p test_install/bin
mkdir -p test_install/share/prokaregia
cp -r * test_install/share/prokaregia/
ln -s ../share/prokaregia/run_prokaregia.sh test_install/bin/prokaregia

# Test the symlinked version
cd test_install
./bin/prokaregia -h
```

**Expected**: Help message should display correctly from symlink.

## Test 6: Conda Environment Creation

Test that individual tool environments can be created:

```bash
# Test creating one of the conda environments
conda env create -f envs/filtlong.yml -n test-filtlong
conda activate test-filtlong
filtlong --help
conda deactivate
conda env remove -n test-filtlong
```

**Expected**: Environment should create successfully and tool should be available.

## Test 7: Output Directory Behavior

Test different output directory specifications:

```bash
# Simple name (should create in current directory)
./run_prokaregia.sh -i test_reads.fastq -s ont -t 2 -o TestRun1

# Absolute path
./run_prokaregia.sh -i test_reads.fastq -s ont -t 2 -o /tmp/TestRun2

# Relative path with slashes (should error)
./run_prokaregia.sh -i test_reads.fastq -s ont -t 2 -o ../TestRun3
```

**Expected**:
- Simple name: Creates `./TestRun1/`
- Absolute path: Creates `/tmp/TestRun2/`
- Relative path: Should show error message

## Test 8: Full Pipeline (Optional)

If you have real data and time, test the complete pipeline:

```bash
# Use a small real dataset (subset of your data)
seqkit head -n 10000 your_real_reads.fastq > small_test.fastq

# Run the full pipeline
./run_prokaregia.sh -i small_test.fastq -s ont -t 8 -o FullTest
```

**Expected**:
- Pipeline should complete all steps
- CheckM2 reports generated
- Quality plots created
- All intermediate files present

## Common Issues and Solutions

### Issue: Conda environments fail to create
**Solution**: Use mamba instead: `./run_prokaregia.sh -i ... -m`

### Issue: CheckM2 database download times out
**Solution**: Download manually first:
```bash
conda activate prokaregia-test
checkm2 database --download
```

### Issue: Out of memory
**Solution**: Use fewer threads: `-t 2` or `-t 4`

### Issue: Snakemake not found
**Solution**: Make sure you're in the correct conda environment with snakemake installed

## Preparing for Bioconda Submission

Once all tests pass:

1. **Create a git tag for the release**:
   ```bash
   git tag -a v1.0.0 -m "Release version 1.0.0"
   git push origin v1.0.0
   ```

2. **Generate tarball hash**:
   ```bash
   wget https://github.com/yourusername/prokaregia_snakemake/archive/v1.0.0.tar.gz
   shasum -a 256 v1.0.0.tar.gz
   ```

3. **Update conda recipe** ([conda_recipe_example/meta.yaml](conda_recipe_example/meta.yaml:10)):
   - Replace `XXXX` with the SHA256 hash
   - Replace `yourusername` with your actual GitHub username
   - Update version number if needed

4. **Test conda build locally** (optional but recommended):
   ```bash
   conda install conda-build
   conda build conda_recipe_example/
   ```

5. **Submit to bioconda-recipes**:
   - Fork https://github.com/bioconda/bioconda-recipes
   - Add your recipe to `recipes/prokaregia/`
   - Submit pull request
   - Wait for CI tests and reviews

## Checklist Before Submission

- [ ] All wrapper script tests pass
- [ ] Dry run completes without errors
- [ ] CheckM2 database downloads correctly
- [ ] Output directories created in correct locations
- [ ] At least one full pipeline run completes successfully
- [ ] LICENSE file present and correct
- [ ] README.md updated and accurate
- [ ] environment.yml tested
- [ ] Git repository tagged with version
- [ ] SHA256 hash added to meta.yaml
- [ ] GitHub username updated in all files
