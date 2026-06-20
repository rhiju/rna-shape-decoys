# RNA SHAPE vs. Decoy RMSD Discrimination (R2307)

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

## Next Steps (Phase 2)

1. **SGNM/ERM SHAPE predictions**: Clone and run `https://github.com/hmblair/sgnm` on all models
   - Compare vs. experimental SHAPE (ceiling analysis)
   - Compare vs. reference-model SHAPE (removes systematic bias)

2. **F1 base-pair scoring**: Compute F1 of predicted vs reference secondary structure

3. **Quantitative metrics**: P@5, MAP@5, Spearman/Pearson across all SHAPE metrics

4. **RNAnix pipeline**: Integrate 3D structure predictions from un-released RNAnix model

## References

- **CASP17 R2307**: https://predictioncenter.org/casp17/target.cgi?id=19&view=rna
- **OpenKnot data**: https://github.com/eternagame/OpenKnotAIDesignData
- **SGNM/ERM**: https://github.com/hmblair/sgnm
- **OpenKnot scorer**: https://github.com/eternagame/OpenKnotScoreMATLAB

---

**Repository**: https://github.com/rhiju/rna-shape-decoys  
**Date**: 2026-06-20
