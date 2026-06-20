#!/usr/bin/env python3
"""
Heatmap of SGNM-predicted SHAPE profiles, sorted by Pearson correlation to the
reference-model SGNM profile.

SGNM raw output has a small dynamic range, so each profile is min-max
normalized to [0,1] for display (reveals the per-structure relative pattern;
correlation used for sorting is computed on the raw values, scale-invariant).

Color: low predicted reactivity = white, high = red.
Run after 03_shape_sgnm.py (needs results/sgnm_profiles.csv).
"""
import matplotlib
matplotlib.rcParams['font.family'] = 'Helvetica'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

CMAP = LinearSegmentedColormap.from_list('white_red', [(1, 1, 1), (1, 0, 0)])

prof = pd.read_csv('results/sgnm_profiles.csv', index_col='model')
ref = prof.loc['reference'].values.astype(float)

# correlation of each profile to reference (raw values)
corrs = {}
for name, row in prof.iterrows():
    if name == 'reference':
        continue
    v = row.values.astype(float)
    corrs[name] = pearsonr(v, ref)[0]

order = sorted(corrs, key=lambda k: corrs[k], reverse=True)

def norm(v):
    lo, hi = v.min(), v.max()
    return (v - lo) / (hi - lo) if hi > lo else v * 0

# reference on top, then sorted models
rows = [norm(ref)] + [norm(prof.loc[n].values.astype(float)) for n in order]
matrix = np.array(rows)
labels = ['REFERENCE'] + [f'{n}  r={corrs[n]:.2f}' for n in order]
N = matrix.shape[1]
n_rows = matrix.shape[0]

fig, ax = plt.subplots(figsize=(16, max(8, n_rows * 0.12)))
im = ax.imshow(matrix, aspect='auto', cmap=CMAP, vmin=0, vmax=1, interpolation='nearest')
ax.set_yticks(range(n_rows))
ax.set_yticklabels(labels, fontsize=4)
ax.set_xlabel('Nucleotide position', fontsize=13, fontweight='bold')
ax.set_title('SGNM-predicted SHAPE profiles (per-row min-max normalized)\n'
             'sorted by correlation to reference-model profile',
             fontsize=15, fontweight='bold')
ax.set_xticks(range(0, N, 10))
ax.set_xticklabels(range(1, N + 1, 10), fontsize=9)
plt.colorbar(im, ax=ax, label='normalized predicted reactivity', fraction=0.02)
plt.tight_layout()
out = 'results/shape_heatmap_sgnm.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
print(f"Saved {out}  ({n_rows} rows x {N} positions)")
import subprocess
subprocess.run(['open', out])
