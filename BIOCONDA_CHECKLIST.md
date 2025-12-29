# Bioconda Submission Checklist

## Pre-submission Checklist

### 1. Replace Placeholders
Before submission, update your GitHub username in these files:
- [ ] [README.md](README.md:40) - Line 40: `git clone` URL
- [ ] [conda_recipe_example/meta.yaml](conda_recipe_example/meta.yaml:9) - Line 9: source URL
- [ ] [conda_recipe_example/meta.yaml](conda_recipe_example/meta.yaml:32) - Lines 32-34: home/doc/dev URLs
- [ ] [conda_recipe_example/meta.yaml](conda_recipe_example/meta.yaml:44) - Line 44: maintainer name
- [ ] [LICENSE](LICENSE:3) - Line 3: Your actual name

### 2. Version and Release
- [ ] Decide on version number (e.g., 1.0.0)
- [ ] Create git tag: `git tag -a v1.0.0 -m "Release v1.0.0"`
- [ ] Push tag: `git push origin v1.0.0`
- [ ] Verify release on GitHub: `https://github.com/YOURUSERNAME/prokaregia_snakemake/releases`

### 3. Generate SHA256 Hash
```bash
# Download the release tarball
wget https://github.com/YOURUSERNAME/prokaregia_snakemake/archive/v1.0.0.tar.gz

# Generate hash
shasum -a 256 v1.0.0.tar.gz

# Copy the hash and update meta.yaml line 10
```
- [ ] SHA256 hash added to [meta.yaml](conda_recipe_example/meta.yaml:10)

### 4. Local Testing
- [ ] Wrapper script help works: `./run_prokaregia.sh -h`
- [ ] Wrapper validates inputs correctly
- [ ] Snakemake dry run passes: `snakemake -n ...`
- [ ] CheckM2 database downloads successfully
- [ ] Output directories created correctly
- [ ] (Optional) Full pipeline run completes

### 5. File Verification
- [ ] [LICENSE](LICENSE) file present with correct copyright holder
- [ ] [README.md](README.md) complete and accurate
- [ ] [environment.yml](environment.yml) tested
- [ ] [run_prokaregia.sh](run_prokaregia.sh) executable (`chmod +x`)
- [ ] All `envs/*.yml` files present and valid
- [ ] [scripts/](scripts/) directory contains all helper scripts

### 6. Documentation Review
- [ ] Installation instructions clear
- [ ] Usage examples work
- [ ] Database setup documented
- [ ] Troubleshooting section helpful
- [ ] Citation information included

## Bioconda Submission Steps

### 1. Fork bioconda-recipes
```bash
# Go to GitHub and fork
https://github.com/bioconda/bioconda-recipes

# Clone your fork
git clone https://github.com/YOURUSERNAME/bioconda-recipes.git
cd bioconda-recipes

# Add upstream remote
git remote add upstream https://github.com/bioconda/bioconda-recipes.git
```

### 2. Create recipe branch
```bash
# Update your fork
git checkout master
git pull upstream master
git push origin master

# Create new branch for your recipe
git checkout -b add-prokaregia
```

### 3. Add your recipe
```bash
# Create recipe directory
mkdir -p recipes/prokaregia

# Copy your meta.yaml
cp /path/to/prokaregia_snakemake/conda_recipe_example/meta.yaml recipes/prokaregia/meta.yaml

# Optional: Add build.sh if needed (usually not for noarch)
# Optional: Add post-link.sh for post-install messages
```

### 4. Test with bioconda-utils (Recommended)
```bash
# Install bioconda-utils
conda install -c bioconda bioconda-utils

# Lint your recipe
bioconda-utils lint recipes/prokaregia/meta.yaml

# Build and test locally
bioconda-utils build recipes/prokaregia/meta.yaml
```

### 5. Commit and push
```bash
git add recipes/prokaregia/
git commit -m "Add prokaregia recipe

ProkaRegia is an automated long-read assembly and binning pipeline
for prokaryotic genomes in contaminated datasets."

git push origin add-prokaregia
```

### 6. Create Pull Request
1. Go to your fork on GitHub
2. Click "Compare & pull request"
3. Title: `Add prokaregia`
4. Description:
   ```
   This PR adds ProkaRegia, an automated long-read assembly and binning
   pipeline for prokaryotic genomes.

   - Package: prokaregia
   - Version: 1.0.0
   - License: MIT
   - Repository: https://github.com/YOURUSERNAME/prokaregia_snakemake

   The pipeline performs quality control, assembly with Flye, taxonomic
   filtering, multiple binning algorithms, refinement, scaffolding,
   polishing, and quality assessment.
   ```
5. Submit PR

### 7. Address Review Comments
- [ ] CI tests pass (Linux, MacOS)
- [ ] Linting passes
- [ ] No conflicts with existing packages
- [ ] Address any reviewer feedback

## Post-Acceptance

### Once merged:
- [ ] Update README with actual conda install command
- [ ] Announce release (Twitter, lab website, etc.)
- [ ] Add conda badge to README:
  ```markdown
  [![Conda](https://img.shields.io/conda/vn/bioconda/prokaregia.svg)](https://anaconda.org/bioconda/prokaregia)
  ```

## Common Bioconda Requirements

### meta.yaml essentials:
- ✅ Valid YAML syntax
- ✅ Correct SHA256 hash
- ✅ Valid source URL (must be publicly accessible)
- ✅ License specified
- ✅ LICENSE file in source
- ✅ All dependencies available in conda-forge or bioconda
- ✅ Test command specified
- ✅ Maintainer GitHub username

### Build requirements:
- ✅ Use `noarch: generic` for pure bash/python scripts
- ✅ Pin dependencies appropriately (not too strict, not too loose)
- ✅ Avoid version conflicts with common packages

### Best practices:
- ✅ Clear, concise description
- ✅ Link to documentation
- ✅ Reproducible builds
- ✅ Fast installation (dependencies not too heavy)

## Resources

- Bioconda documentation: https://bioconda.github.io/
- Bioconda recipes repository: https://github.com/bioconda/bioconda-recipes
- Conda build docs: https://docs.conda.io/projects/conda-build/
- Example recipes: Browse existing packages in bioconda-recipes

## Getting Help

- Bioconda Gitter chat: https://gitter.im/bioconda/Lobby
- GitHub issues: https://github.com/bioconda/bioconda-recipes/issues
- Bioconda docs FAQ: https://bioconda.github.io/faqs.html
