#!/usr/bin/env python3
"""
Ingest ERM profiles produced on a GPU (results/erm_profiles.csv) and merge
ERM SHAPE correlations into results/scores.csv, then regenerate plots, metrics,
and the HTML report.

  shape_erm_vs_expt : Pearson r of ERM SHAPE vs experimental SHAPE
  shape_erm_vs_ref  : Pearson r of ERM SHAPE vs reference-model ERM SHAPE

Run locally after downloading erm_profiles.csv from Colab:
    python3 scripts/ingest_erm.py
"""
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import pearsonr

prof = pd.read_csv('results/erm_profiles.csv', index_col='model')
ref = prof.loc['reference'].values.astype(float)
exp = pd.read_csv('data/experimental_shape.csv')['shape'].values
m_exp = ~np.isnan(exp)

vs_expt, vs_ref = {}, {}
for name, row in prof.iterrows():
    if name == 'reference':
        continue
    v = row.values.astype(float)
    vs_expt[name] = pearsonr(v[m_exp], exp[m_exp])[0]
    vs_ref[name] = pearsonr(v, ref)[0]

scores = pd.read_csv('results/scores.csv')
key = scores['model_path'].apply(lambda p: Path(p).stem)
scores['shape_erm_vs_expt'] = key.map(vs_expt)
scores['shape_erm_vs_ref'] = key.map(vs_ref)
scores.to_csv('results/scores.csv', index=False)
print("Merged ERM columns into results/scores.csv")

# reference ERM vs experimental (ceiling sanity check)
print(f"Reference ERM vs experimental: r={pearsonr(ref[m_exp], exp[m_exp])[0]:.3f}")

# regenerate downstream artifacts
for script in ['scripts/08_metrics.py', 'scripts/04_plot.py', 'scripts/07_make_report.py']:
    print(f"--- running {script} ---")
    subprocess.run(['python3', script], check=False)
