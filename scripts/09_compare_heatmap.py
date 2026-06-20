#!/usr/bin/env python3
"""
Focused comparison heatmap: experimental SHAPE vs SGNM predictions for the
reference structure, the best-correlating GOOD (low-RMSD) FARFAR2 decoy, and the
best-correlating BAD (high-RMSD) FARFAR2 decoy. Dot-bracket secondary structure
is drawn on each structure row.

Each row is min-max normalized for display (white=low, red=high reactivity), so
the different scales (experimental vs predicted) are visually comparable.
"""
import matplotlib
matplotlib.rcParams['font.family'] = 'Helvetica'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import sys
from scipy.stats import pearsonr

sys.path.insert(0, 'scripts')
from dssr_util import dssr_ss

CMAP = LinearSegmentedColormap.from_list('white_red', [(1, 1, 1), (1, 0, 0)])
REF_PDB = 'data/farfar2/Mol9_reference_UtoG_buildloop.pdb'
GOOD = 'Mol9.out.rms.4'   # rmsd 5.3 A, best SGNM-vs-expt corr among good
BAD = 'Mol9.out.9'        # rmsd 20.9 A, best SGNM-vs-expt corr among bad

exp = pd.read_csv('data/experimental_shape.csv')['shape'].values.astype(float)
prof = pd.read_csv('results/sgnm_profiles.csv', index_col='model')
N = len(exp)
m = ~np.isnan(exp)

ref_ss = dssr_ss(REF_PDB)
good_ss = dssr_ss(f'data/farfar2/{GOOD}.pdb')
bad_ss = dssr_ss(f'data/farfar2/{BAD}.pdb')


def norm(v):
    lo, hi = np.nanmin(v), np.nanmax(v)
    return (v - lo) / (hi - lo) if hi > lo else v * 0


def corr(name):
    v = prof.loc[name].values.astype(float)
    return pearsonr(v[m], exp[m])[0]

# rows: (display values, dot-bracket SS, label)
rows = [
    (exp, ref_ss, 'Experimental 2A3 SHAPE  (target structure SS)'),
    (prof.loc['reference'].values.astype(float), ref_ss,
     f'Reference (true) structure  SGNM  r={corr("reference"):.2f}'),
    (prof.loc[GOOD].values.astype(float), good_ss,
     f'GOOD decoy {GOOD}  RMSD 5.3Å  SGNM  r={corr(GOOD):.2f}'),
    (prof.loc[BAD].values.astype(float), bad_ss,
     f'BAD decoy {BAD}  RMSD 20.9Å  SGNM  r={corr(BAD):.2f}'),
]
matrix = np.vstack([norm(r[0]) for r in rows])
ss_strings = [r[1] for r in rows]
labels = [r[2] for r in rows]

fig, ax = plt.subplots(figsize=(22, 4.5))
ax.imshow(matrix, aspect='auto', cmap=CMAP, vmin=0, vmax=1, interpolation='nearest')

# dot-bracket characters on each cell
for ri, ss in enumerate(ss_strings):
    if ss is None or len(ss) != N:
        continue
    for c in range(N):
        ax.text(c, ri, ss[c], ha='center', va='center', fontsize=7,
                fontfamily='monospace', color='black')

ax.set_yticks(range(len(labels)))
ax.set_yticklabels(labels, fontsize=12, fontfamily='monospace')
ax.set_xlabel('Nucleotide position', fontsize=16, fontweight='bold')
ax.set_title('Experimental vs SGNM-predicted SHAPE: reference, good & bad FARFAR2 decoys\n'
             '(per-row normalized; white=low, red=high; dot-bracket SS overlaid)',
             fontsize=16, fontweight='bold')
ax.set_xticks([p - 1 for p in range(10, N + 1, 10)])
ax.set_xticklabels(range(10, N + 1, 10), fontsize=12)

cbar = fig.colorbar(ax.images[0], ax=ax, orientation='vertical', fraction=0.015, pad=0.01,
                    ticks=[0, 1])
cbar.ax.set_yticklabels(['low', 'high'], fontsize=12)

plt.tight_layout()
out = 'results/compare_heatmap.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
print(f"Saved {out}")
import subprocess
subprocess.run(['open', out])
