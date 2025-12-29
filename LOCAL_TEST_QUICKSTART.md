# Quick Start: Local Conda Package Testing

This is the **fastest way** to test your conda package locally before submission.

## Step 1: Install conda-build

```bash
conda install conda-build
```

## Step 2: Build the Package Locally

```bash
cd prokaregia_snakemake

# Build using the local source (no GitHub release needed)
conda build conda_recipe_example/meta_local.yaml

# Watch for errors in the build output
# Should end with: "TEST END: prokaregia-1.0.0-0"
```

**What this does:**
- Creates a temporary build environment
- Copies all your files to `$PREFIX/share/prokaregia/`
- Creates symlink from `$PREFIX/bin/prokaregia` to the wrapper script
- Runs the test command: `prokaregia -h`
- Creates a conda package file (`.tar.bz2`)

**If build succeeds**, you'll see:
```
...
TEST START: prokaregia-1.0.0-0
+ prokaregia -h
ProkaRegia - Automated long-read assembly and binning pipeline
...
TEST END: prokaregia-1.0.0-0
```

## Step 3: Install the Package in a Test Environment

```bash
# Create fresh test environment
conda create -n test-prokaregia
conda activate test-prokaregia

# Install your locally-built package
conda install --use-local prokaregia

# Verify installation
which prokaregia
# Should show: /path/to/conda/envs/test-prokaregia/bin/prokaregia

prokaregia -h
# Should show the help message
```

## Step 4: Test the Package Works

```bash
# Create tiny test data
echo -e "@read1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII" > test_reads.fastq

# Test the command works
prokaregia -i test_reads.fastq -s ont -t 2

# This will likely fail on the tiny test data, but you should see:
# 1. Configuration summary printed
# 2. Snakemake starting
# 3. Conda environments being created
```

## Step 5: Test with Real Data (Optional)

```bash
# Use a small subset of real data
head -40000 /path/to/your/real_reads.fastq > small_test.fastq

# Run the pipeline
prokaregia -i small_test.fastq -s ont -t 4 -o TestOutput

# Let it run for a while to verify:
# - CheckM2 database downloads correctly
# - Individual tool environments install
# - First few steps complete successfully
```

## Step 6: Clean Up

```bash
# Exit test environment
conda deactivate

# Remove test environment
conda env remove -n test-prokaregia

# Remove test data
rm test_reads.fastq small_test.fastq
rm -rf TestOutput/
```

## Common Issues and Fixes

### Issue: Build fails with "command not found: conda-build"
**Fix:**
```bash
conda install conda-build
```

### Issue: "Error: source.url not accessible"
**Solution:** You're using `meta.yaml` instead of `meta_local.yaml`
```bash
# Use the local version for testing
conda build conda_recipe_example/meta_local.yaml
```

### Issue: Test command fails - "prokaregia: command not found"
**Possible causes:**
1. Symlink not created correctly - check build script in meta.yaml
2. Script not executable - add `chmod +x` in build script
3. Wrong path - verify `$PREFIX/share/prokaregia/` contains files

**Debug:**
```bash
# Check what got installed
conda list -n test-prokaregia prokaregia

# Check file locations
ls -la $(conda info --base)/envs/test-prokaregia/bin/prokaregia
ls -la $(conda info --base)/envs/test-prokaregia/share/prokaregia/
```

### Issue: "snakemake: command not found" when running prokaregia
**Fix:** Add snakemake-minimal to environment
```bash
conda install snakemake-minimal
```
This should have been installed automatically - check your meta.yaml `requirements.run` section.

### Issue: Wrapper script can't find Snakefile
**Check:**
```bash
# Verify Snakefile is in the installed location
ls $(conda info --base)/envs/test-prokaregia/share/prokaregia/Snakefile

# Check wrapper script path resolution
cat $(conda info --base)/envs/test-prokaregia/bin/prokaregia | grep SCRIPT_DIR
```

### Issue: Build succeeds but install fails
**Possible cause:** Package not in local channel
```bash
# Check local packages
conda search --override-channels --channel local prokaregia

# If missing, rebuild
conda build conda_recipe_example/meta_local.yaml

# The package gets stored at:
# Linux: ~/miniconda3/conda-bld/linux-64/prokaregia-1.0.0-0.tar.bz2
# Mac: ~/miniconda3/conda-bld/osx-64/prokaregia-1.0.0-0.tar.bz2
```

## After Local Testing Succeeds

Once everything works locally:

1. **Create GitHub release** (see [BIOCONDA_CHECKLIST.md](BIOCONDA_CHECKLIST.md))
2. **Generate SHA256 hash** from release tarball
3. **Update meta.yaml** (not meta_local.yaml) with:
   - Real GitHub source URL
   - SHA256 hash
   - Your GitHub username
4. **Test building from GitHub**:
   ```bash
   conda build conda_recipe_example/meta.yaml
   ```
5. **Submit to bioconda** (follow [BIOCONDA_CHECKLIST.md](BIOCONDA_CHECKLIST.md))

## Quick Test Commands Summary

```bash
# Build
conda build conda_recipe_example/meta_local.yaml

# Install
conda create -n test-prokaregia
conda activate test-prokaregia
conda install --use-local prokaregia

# Test
prokaregia -h
prokaregia -i test_reads.fastq -s ont -t 2

# Clean up
conda deactivate
conda env remove -n test-prokaregia
```

## What to Test Before Considering It Ready

- [ ] `conda build` completes without errors
- [ ] Test command `prokaregia -h` succeeds
- [ ] Can install with `conda install --use-local prokaregia`
- [ ] `prokaregia -h` shows help after install
- [ ] `prokaregia -i ... -s ...` accepts arguments correctly
- [ ] Snakemake command is built correctly (check config summary)
- [ ] Snakefile is found and workflow starts
- [ ] (Optional) First few pipeline steps complete with real data

Once all these pass, you're ready to create a GitHub release and submit to bioconda!
