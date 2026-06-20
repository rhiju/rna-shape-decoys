#!/usr/bin/env python3
"""
Naive SHAPE score: correlate each model's DSSR paired/unpaired vector (0/0.5)
against (a) the EXPERIMENTAL SHAPE profile and (b) the reference structure's
paired/unpaired vector.

Builds results/scores.csv with:
  model_path, source, group, shape_naive_vs_expt, shape_naive_vs_ref
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

sys.path.insert(0, 'scripts')
from dssr_util import dssr_ss, unpaired_vec

REF_PDB = 'data/farfar2/Mol9_reference_UtoG_buildloop.pdb'

# experimental SHAPE (100 nt, aligned to structure positions 1..100)
exp = pd.read_csv('data/experimental_shape.csv')
exp_shape = pd.to_numeric(exp['shape'], errors='coerce').values
N = len(exp_shape)
print(f"Experimental SHAPE: {N} positions")

# reference paired/unpaired
ref_vec = unpaired_vec(dssr_ss(REF_PDB))
print(f"Reference SS length: {len(ref_vec)}")

model_files = sorted(Path('data/farfar2').glob('Mol9.out.*.pdb'))
model_files += sorted(Path('data/casp17/R2307').glob('R2307TS*.pdb'))
print(f"Scoring {len(model_files)} models...")


def corr_vs_expt(vec):
    m = ~np.isnan(exp_shape)
    return pearsonr(vec[m], exp_shape[m])[0] if m.sum() > 2 else np.nan


rows = []
for i, pdb in enumerate(model_files):
    if i % 50 == 0:
        print(f"  {i}/{len(model_files)}")
    source = 'farfar2' if 'farfar2' in str(pdb) else 'casp17'
    group = pdb.stem[:10] if source == 'casp17' else 'FARFAR2'
    ss = dssr_ss(str(pdb))
    naive_expt = naive_ref = np.nan
    if ss is not None and len(ss) == N:
        vec = unpaired_vec(ss)
        naive_expt = corr_vs_expt(vec)
        naive_ref = pearsonr(vec, ref_vec)[0]
    rows.append({'model_path': str(pdb), 'source': source, 'group': group,
                 'shape_naive_vs_expt': naive_expt, 'shape_naive_vs_ref': naive_ref})

df = pd.DataFrame(rows)
df.to_csv('results/scores.csv', index=False)
print(f"Saved results/scores.csv ({len(df)} models)")
print(f"  naive_vs_expt computed for {df['shape_naive_vs_expt'].notna().sum()} models")

# reference self-scores (for plot reference point)
ref_naive_vs_expt = corr_vs_expt(ref_vec)
print(f"  reference naive vs expt: r={ref_naive_vs_expt:.3f}")
pd.DataFrame([{'metric': 'naive_vs_expt', 'ref_value': ref_naive_vs_expt},
              {'metric': 'naive_vs_ref', 'ref_value': 1.0}]).to_csv(
    'results/reference_point_naive.csv', index=False)
