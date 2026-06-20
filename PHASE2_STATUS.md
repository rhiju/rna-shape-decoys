# Phase 2 Status Report

## ✅ Completed

### 1. Environment & Dependencies
- ✅ **ciffy**: Installed and working (with pyproject.toml fix)
- ✅ **dlu-torch**: Installed and exports all required classes
- ✅ **sgnm**: Installed, models.py imports correctly
- ✅ All repos from correct source: https://github.com/hmblair/

### 2. Visualizations
- ✅ **Heatmap of SHAPE profiles** (scripts/05_shape_heatmap.py)
  - Shows naive SHAPE predictions sorted by Pearson correlation to experimental
  - White (paired, 0) to Red (unpaired, 1) colormap
  - Secondary structures shown below in dot-bracket notation
  - 174 models with valid SHAPE data
  - Results in: `results/shape_heatmap_naive.png`

### 3. SGNM Integration
- ✅ **Prediction script ready** (scripts/03_shape_sgnm.py)
  - Implements both vs-experimental and vs-reference correlations
  - Ready to run once checkpoints are available
  - Will add `shape_sgnm_vs_expt` and `shape_sgnm_vs_ref` columns to results

## ⏳ Blocked: Checkpoint Files

SGNM model weights need to be downloaded from GitHub releases. The automated download (wget) returns only 9 bytes (redirect page). 

### To download manually:
```bash
cd ~/src/sgnm

# Try these URLs (check GitHub releases for current versions):
# https://github.com/hmblair/sgnm/releases/download/v1.0.0/gnm-checkpoint.pth
# https://github.com/hmblair/sgnm/releases/download/v1.0/gnm-checkpoint.pth

# Or visit releases page:
# https://github.com/hmblair/sgnm/releases
```

Once checkpoints are in place:
```bash
python3 scripts/03_shape_sgnm.py
```

## Issues Found (for upstream PRs)

### 1. **ciffy** (`https://github.com/hmblair/ciffy`)
- **Fixed**: Duplicate `fallback_version` in `[tool.setuptools_scm]` section
  - Line 64 had `fallback_version = "0.0.0"` (duplicate)
  - Line 67 had `fallback_version = "0.0.0.dev0"`
  - **Solution**: Remove line 64, keep line 67 only

### 2. **dlu-torch** (`https://github.com/hmblair/dlu`)
- **Note**: Package name mismatch - published as `dlu-torch` but imported as `dlu`
  - This is intentional (namespace package with branch variants)
  - No fix needed, just update documentation

### 3. **sgnm** (`https://github.com/hmblair/sgnm`)
- **Issue 1**: models.py line 14 has incorrect import
  ```python
  # Current (broken):
  from dlu import RadialBasisFunctions, DenseNetwork
  
  # Should be (for dlu-torch):
  from dlu.modules import RadialBasisFunctions, DenseNetwork
  ```
  - **Note**: This is already fixed in this repo
  
- **Issue 2**: pyproject.toml line 20 requires "dlu" but should be "dlu-torch"
  - Changed in this repo; PR pending

- **Issue 3**: Equivariant model requires ciffy submodules (`Scale`, `Reduction`) not yet available
  - Not critical for GNM model (which we use)
  - Likely resolvable with ciffy updates

### 4. **Model Checkpoint Distribution**
- Direct GitHub release downloads fail (return 301 redirect)
- Users must manually download or use `gh release download` command
- **Recommendation**: Use `gh-cli` instructions in docs

## Next Steps

1. **Once checkpoints are available:**
   ```bash
   python3 scripts/03_shape_sgnm.py  # Will add SGNM predictions
   ```

2. **Create Enhanced Heatmap:**
   - Repeat heatmap with SGNM SHAPE (vs experimental + vs reference)
   - Compare naive vs SGNM discrimination power

3. **Quantitative Metrics:**
   - P@5, MAP@5 for SHAPE-based ranking
   - Spearman correlations across SHAPE metrics

4. **Submit PRs:**
   - ciffy: Fix duplicate fallback_version
   - sgnm: Update dlu dependency and imports
   - dlu: Update README with installation notes for namespace package

## Current Results (Naive SHAPE from DSSR)

- **Spearman ρ = -0.464** (p = 1.14e-10)
- Significant negative correlation: lower RMSD → higher SHAPE correlation
- SGNM predictions expected to improve discrimination power
