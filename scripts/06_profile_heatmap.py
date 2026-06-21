#!/usr/bin/env python3
"""
Heatmap of predicted SHAPE profiles, sorted by Pearson correlation to the
reference-model profile. Works for either SGNM or ERM.

Raw output has a small dynamic range, so each profile is min-max normalized to
[0,1] for display (reveals the per-structure relative pattern; the correlation
used for sorting is on the raw values, scale-invariant). Low=white, high=red.

Usage:
    python3 scripts/06_profile_heatmap.py sgnm   # -> results/shape_heatmap_sgnm.png
    python3 scripts/06_profile_heatmap.py erm    # -> results/shape_heatmap_erm.png
"""
import sys
import matplotlib
matplotlib.rcParams['font.family'] = 'Helvetica'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

KIND = sys.argv[1] if len(sys.argv) > 1 else 'sgnm'
LABELS = {'sgnm': 'SGNM', 'erm': 'ERM'}
TAG = LABELS.get(KIND, KIND.upper())

CMAP = LinearSegmentedColormap.from_list('white_red', [(1, 1, 1), (1, 0, 0)])
prof = pd.read_csv(f'results/{KIND}_profiles.csv', index_col='model')
ref = prof.loc['reference'].values.astype(float)

corrs = {}
for name, row in prof.iterrows():
    if name == 'reference':
        continue
    corrs[name] = pearsonr(row.values.astype(float), ref)[0]
order = sorted(corrs, key=lambda k: corrs[k], reverse=True)


def norm(v):
    lo, hi = v.min(), v.max()
    return (v - lo) / (hi - lo) if hi > lo else v * 0


rows = [norm(ref)] + [norm(prof.loc[n].values.astype(float)) for n in order]
matrix = np.array(rows)
labels = ['REFERENCE'] + [f'{n}  r={corrs[n]:.2f}' for n in order]
N, n_rows = matrix.shape[1], matrix.shape[0]

fig, ax = plt.subplots(figsize=(16, max(8, n_rows * 0.12)))
im = ax.imshow(matrix, aspect='auto', cmap=CMAP, vmin=0, vmax=1, interpolation='nearest')
ax.set_yticks(range(n_rows))
ax.set_yticklabels(labels, fontsize=4)
ax.set_xlabel('Nucleotide position', fontsize=13, fontweight='bold')
ax.set_title(f'{TAG}-predicted SHAPE profiles (per-row min-max normalized)\n'
             'sorted by correlation to reference-model profile',
             fontsize=15, fontweight='bold')
ax.set_xticks([p - 1 for p in range(10, N + 1, 10)])
ax.set_xticklabels(range(10, N + 1, 10), fontsize=9)
plt.colorbar(im, ax=ax, label='normalized predicted reactivity', fraction=0.02)
plt.tight_layout()
out = f'results/shape_heatmap_{KIND}.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
print(f"Saved {out}  ({n_rows} rows x {N} positions)")
import subprocess
subprocess.run(['open', out])
