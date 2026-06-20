#!/usr/bin/env python3
"""
Run SGNM (GNM) SHAPE prediction on all structures.

Saves:
  - results/sgnm_profiles.csv : per-residue SGNM SHAPE for every model + reference
  - updates results/scores.csv with:
      shape_sgnm_vs_ref  : Pearson r of model SGNM SHAPE vs reference-model SGNM SHAPE (ceiling)
      shape_sgnm_vs_expt : Pearson r vs experimental SHAPE (if available, else blank)

Uses scripts/sgnm_predict.py which bypasses ciffy's mmCIF loader.
"""
import sys
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

sys.path.insert(0, 'scripts')
import sgnm
from sgnm_predict import predict

CHECKPOINT = '/Users/rhiju/src/sgnm/gnm-checkpoint.pth'
REF_PDB = 'data/farfar2/Mol9_reference_UtoG_buildloop.pdb'

print("Loading SGNM GNM model...")
model = sgnm.load(CHECKPOINT)
model.eval()

# Reference profile (ceiling reference)
print("Predicting reference profile...")
ref_pred, ref_kept = predict(model, REF_PDB)
ref_shape = np.array(ref_pred[:, 0].tolist())
N = len(ref_shape)
print(f"  Reference: {N} residues")

# Experimental SHAPE (may be all-NaN placeholder)
exp_df = pd.read_csv('data/experimental_shape.csv')
exp_shape = pd.to_numeric(exp_df['shape'], errors='coerce').values
has_expt = not np.all(np.isnan(exp_shape))
print(f"  Experimental SHAPE available: {has_expt}")

# Reference model's own SGNM-vs-experimental correlation (the ceiling / sanity check)
if has_expt and len(ref_shape) == len(exp_shape):
    _m = ~np.isnan(exp_shape)
    ref_vs_expt = pearsonr(ref_shape[_m], exp_shape[_m])[0]
    print(f"  Reference SGNM SHAPE vs experimental: r={ref_vs_expt:.3f}  (ceiling)")

# Gather all model PDBs
model_files = sorted(Path('data/farfar2').glob('Mol9.out.*.pdb'))
model_files += sorted(Path('data/casp17/R2307').glob('R2307TS*.pdb'))
print(f"Predicting SGNM SHAPE for {len(model_files)} models...")

profiles = {'reference': ref_shape}
rows = []
for i, pdb in enumerate(model_files):
    if i % 40 == 0:
        print(f"  {i}/{len(model_files)}")
    name = pdb.stem
    source = 'farfar2' if 'farfar2' in str(pdb) else 'casp17'
    try:
        pred, kept = predict(model, str(pdb))
        shape = np.array(pred[:, 0].tolist())
    except Exception as e:
        print(f"    skip {name}: {e}")
        continue

    # vs reference ceiling
    if len(shape) == len(ref_shape):
        r_ref = pearsonr(shape, ref_shape)[0]
    else:
        r_ref = np.nan
    # vs experimental
    if has_expt and len(shape) == len(exp_shape):
        mask = ~np.isnan(exp_shape)
        r_expt = pearsonr(shape[mask], exp_shape[mask])[0] if mask.sum() > 2 else np.nan
    else:
        r_expt = np.nan

    profiles[name] = shape
    rows.append({'model': name, 'source': source,
                 'shape_sgnm_vs_ref': r_ref, 'shape_sgnm_vs_expt': r_expt})

print(f"Predicted {len(rows)} models")

# Save per-residue profiles (rows=models, cols=positions)
prof_df = pd.DataFrame({k: v for k, v in profiles.items() if len(v) == N}).T
prof_df.columns = [f'pos_{i+1}' for i in range(N)]
prof_df.index.name = 'model'
prof_df.to_csv('results/sgnm_profiles.csv')
print(f"Saved results/sgnm_profiles.csv ({prof_df.shape[0]} profiles x {N} positions)")

# Merge correlations into scores.csv
sgnm_df = pd.DataFrame(rows).set_index('model')
scores = pd.read_csv('results/scores.csv')
scores['_key'] = scores['model_path'].apply(lambda p: Path(p).stem)
scores['shape_sgnm_vs_ref'] = scores['_key'].map(sgnm_df['shape_sgnm_vs_ref'])
scores['shape_sgnm_vs_expt'] = scores['_key'].map(sgnm_df['shape_sgnm_vs_expt'])
scores = scores.drop(columns='_key')
scores.to_csv('results/scores.csv', index=False)
print("Updated results/scores.csv with shape_sgnm_vs_ref, shape_sgnm_vs_expt")

# Summary
vr = sgnm_df['shape_sgnm_vs_ref'].dropna()
print(f"\nSGNM vs reference: mean r={vr.mean():.3f}, range [{vr.min():.3f}, {vr.max():.3f}]")
