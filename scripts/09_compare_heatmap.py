#!/usr/bin/env python3
"""
Comparison heatmap, 3 stacked panels (naive / SGNM / ERM). Each panel shows the
experimental SHAPE and the predictor's profile for the reference structure, the
best-correlating GOOD (low-RMSD) FARFAR2 decoy, and the best-correlating BAD
(high-RMSD) FARFAR2 decoy. Dot-bracket secondary structure is drawn on each
structure row. Each row is min-max normalized (white=low, red=high).
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
from dssr_util import dssr_ss, unpaired_vec

CMAP = LinearSegmentedColormap.from_list('white_red', [(1, 1, 1), (1, 0, 0)])
REF_PDB = 'data/farfar2/Mol9_reference_UtoG_buildloop.pdb'
GOOD = 'Mol9.out.rms.4'   # rmsd 5.3 A
BAD = 'Mol9.out.9'        # rmsd 20.9 A

exp = pd.read_csv('data/experimental_shape.csv')['shape'].values.astype(float)
sgnm = pd.read_csv('results/sgnm_profiles.csv', index_col='model')
erm = pd.read_csv('results/erm_profiles.csv', index_col='model')
N = len(exp)
m = ~np.isnan(exp)

# (model name, dot-bracket SS, short tag)
STRUCTS = [
    ('reference', dssr_ss(REF_PDB), 'Reference (true)'),
    (GOOD, dssr_ss(f'data/farfar2/{GOOD}.pdb'), f'GOOD {GOOD} 5.3Å'),
    (BAD, dssr_ss(f'data/farfar2/{BAD}.pdb'), f'BAD {BAD} 20.9Å'),
]
ref_ss = STRUCTS[0][1]


def norm(v):
    lo, hi = np.nanmin(v), np.nanmax(v)
    return (v - lo) / (hi - lo) if hi > lo else v * 0


def profile(predictor, name, ss):
    if predictor == 'naive':
        return unpaired_vec(ss)
    table = sgnm if predictor == 'sgnm' else erm
    return table.loc[name].values.astype(float)


def cc(v):
    return pearsonr(v[m], exp[m])[0]


PANELS = [('naive', 'Naive (paired/unpaired)'), ('sgnm', 'SGNM'), ('erm', 'ERM')]

fig, axes = plt.subplots(len(PANELS), 1, figsize=(22, 12))
for ax, (pred, title) in zip(axes, PANELS):
    # rows: experimental, then each structure's predicted profile
    rows = [(exp, ref_ss, 'Experimental 2A3 SHAPE')]
    for name, ss, tag in STRUCTS:
        v = profile(pred, name, ss)
        rows.append((v, ss, f'{tag}  r={cc(v):.2f}'))
    matrix = np.vstack([norm(r[0]) for r in rows])
    ax.imshow(matrix, aspect='auto', cmap=CMAP, vmin=0, vmax=1, interpolation='nearest')
    for ri, (_, ss, _) in enumerate(rows):
        if ss and len(ss) == N:
            for c in range(N):
                ax.text(c, ri, ss[c], ha='center', va='center', fontsize=6,
                        fontfamily='monospace', color='black')
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[2] for r in rows], fontsize=11, fontfamily='monospace')
    ax.set_title(title, fontsize=15, fontweight='bold', loc='left')
    ax.set_xticks([p - 1 for p in range(10, N + 1, 10)])
    ax.set_xticklabels(range(10, N + 1, 10), fontsize=11)

axes[-1].set_xlabel('Nucleotide position', fontsize=16, fontweight='bold')
fig.suptitle('Experimental vs predicted SHAPE (naive / SGNM / ERM) — reference, '
             'good & bad FARFAR2 decoys\n(per-row normalized: white=low, red=high; '
             'dot-bracket SS overlaid; r = corr to experiment)',
             fontsize=16, fontweight='bold')
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = 'results/compare_heatmap.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
print(f"Saved {out}")
import subprocess
subprocess.run(['open', out])
