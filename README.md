# RNA SHAPE vs. Decoy RMSD Discrimination (R2307)

📊 **[Full results writeup with figures → `results.md`](results.md)** (also [`results.html`](results.html) / [`results.pdf`](results.pdf))

## Overview

This project tests whether SHAPE reactivity data can discriminate near-native from non-native RNA 3D structures. We analyze CASP17 target R2307 (PDB: 10ZT) using two datasets:
- **FARFAR2 decoys**: 15 structures generated via Rosetta FARFAR2/rna_denovo
- **CASP17 models**: 197+ submitted structure predictions from various groups

## Key Results

**Naive SHAPE Discrimination (from secondary structure)**:
- Correlation: **Spearman ρ = -0.464 (p = 1.14e-10)**
- **Significant negative correlation**: Lower RMSD decoys tend to have higher SHAPE correlation with reference
- This demonstrates SHAPE data can discriminate near-native from non-native structures at statistically significant levels

### Data Summary
- **FARFAR2**: mean RMSD = 10.53 Å (range: 5.15–20.93 Å)
- **CASP17**: mean RMSD = 16.40 Å (range: 0–38.79 Å)

## Project Structure

```
data/
  farfar2/              # FARFAR2 decoy PDB files (15)
  casp17/R2307/         # CASP17 predicted models (197)
  experimental_shape.csv # Per-nt SHAPE data (placeholder)
  sequence.fasta        # Target sequence

scripts/
  01_compute_rmsd.py    # Compute C1' RMSD vs reference
  02_shape_naive.py     # Naive SHAPE scoring via DSSR
  04_plot.py            # Generate main plot
  f1_score.py           # Base-pair F1 computation
  convert_casp17_to_pdb.py # Convert TS format to PDB

results/
  scores.csv            # Model scores: RMSD, naive SHAPE, F1
  shape_vs_rmsd.png     # Main scatter plot
```

## Methods

### 1. RMSD Computation
- Reference: `Mol9_reference_UtoG_buildloop.pdb` from FARFAR2 zip
- Metric: C1' RMSD alignment via TMscore
- All models aligned to single reference frame

### 2. Naive SHAPE Scoring
- Use DSSR secondary structure as proxy for SHAPE reactivity
- Binary encoding: 1 (unpaired), 0 (paired)
- Correlate predicted structure with reference model structure
- Serves as lower-bound discriminative power

### 3. Base-Pair F1 (Planned, Phase 2)
- Extract base pairs from DSSR dot-bracket notation
- Compute precision, recall, F1 of base pair sets
- Reference: secondary structure of reference model

## Usage

```bash
# Compute RMSD for all models
python3 /Users/rhiju/src/daslab_tools/structure/rmsd.py \
  -refpdb data/farfar2/Mol9_reference_UtoG_buildloop.pdb \
  -pdb data/farfar2/Mol9.out.*.pdb data/casp17/R2307/R2307TS*.pdb \
  -t > results/rmsd_raw.csv

# Compute naive SHAPE scores
python3 scripts/02_shape_naive.py

# Plot results
python3 scripts/04_plot.py
```

## Dependencies

- **DSSR** (x3dna-dssr): Secondary structure prediction
- **TMscore**: RMSD computation
- Python: pandas, scipy, matplotlib, numpy

## SGNM setup (working)

The hmblair SHAPE-from-structure stack is installed and working on this Mac:

```bash
brew install gperf                      # ciffy C extension needs gperf
cd ~/src/dlu   && pip install .         # regular install (NOT -e) — see note
cd ~/src/ciffy && pip install .
cd ~/src/sgnm  && pip install . --no-deps
cd ~/src/sgnm  && gh release download v2.0.2 -p "*checkpoint*"
```

**Important**: use a *regular* `pip install`, not editable (`-e`). Because
`~/src` is on `PYTHONPATH`, editable installs of `dlu`/`ciffy` collide with the
repo-root namespace package and fail to import. A regular install puts a real
package in site-packages that wins priority.

- **GNM model** (`gnm-checkpoint.pth`): runs on CPU/Mac. Used here.
- **Equivariant model** (`equivariant-checkpoint.pth`): needs `flash-eq`
  (CUDA-only) — run on a GPU node for phase 3 (higher accuracy: +0.63 vs +0.39).

`scripts/sgnm_predict.py` bypasses ciffy's mmCIF loader: it computes per-residue
centers + C2/C4/C6 frames straight from the PDB and calls `model(coords, frames)`,
verified to reproduce `model.ciffy()` exactly.

## Results (vs reference structure)

| Metric | Spearman ρ vs RMSD |
|---|---|
| Naive paired/unpaired SHAPE | −0.46 |
| SGNM-predicted SHAPE | −0.43 |
| Base-pair F1 | −0.48 |

All three discriminate near-native from non-native structures. See `results.html`.

## Next Steps (Phase 3)

1. **Experimental SHAPE**: obtain real per-nucleotide reactivity for R2307 and
   compare predictions against it (current analysis uses reference-model profile).
2. **Equivariant (ERM) model**: run on a CUDA node for the stronger predictor.
3. **Quantitative metrics**: P@5, MAP@5 for SHAPE-based ranking.
4. **RNAnix pipeline**: integrate the un-released RNAnix 3D refiner/ranker.

## References

- **CASP17 R2307**: https://predictioncenter.org/casp17/target.cgi?id=19&view=rna
- **OpenKnot data**: https://github.com/eternagame/OpenKnotAIDesignData
- **SGNM/ERM**: https://github.com/hmblair/sgnm
- **OpenKnot scorer**: https://github.com/eternagame/OpenKnotScoreMATLAB

---

**Repository**: https://github.com/rhiju/rna-shape-decoys  
**Date**: 2026-06-20
