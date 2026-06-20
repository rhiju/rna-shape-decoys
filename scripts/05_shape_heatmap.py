#!/usr/bin/env python3
"""
Heatmap of secondary-structure SHAPE proxies, with the dot-bracket character
drawn ON each cell.

- Color: paired = white (0.0), unpaired = light red (0.5)   [vmin=0, vmax=0.5]
- Cell text: the actual dot-bracket character ( . ) [ ] etc.
- Rows: experimental/reference at top, then each unique secondary structure,
        sorted by Pearson correlation (paired/unpaired vector) to the reference.
- Duplicate secondary structures are collapsed; the row label shows how many
  models share that structure and the best (lowest) RMSD among them.

Run after 01_compute_rmsd.py (needs results/scores.csv).
"""
import matplotlib
matplotlib.rcParams['font.family'] = 'Helvetica'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from pathlib import Path
import subprocess, sys
from scipy.stats import pearsonr

sys.path.insert(0, 'scripts')
from dssr_util import dssr_ss, unpaired_vec

CMAP = LinearSegmentedColormap.from_list('white_red', [(1, 1, 1), (1, 0, 0)])
REF_PDB = 'data/farfar2/Mol9_reference_UtoG_buildloop.pdb'


# --- reference ---
ref_ss = dssr_ss(REF_PDB)
ref_vec = unpaired_vec(ref_ss)
N = len(ref_ss)
print(f"Reference ({N} nt): {ref_ss}")

# --- all models with RMSD from scores.csv ---
scores = pd.read_csv('results/scores.csv')
scores['rmsd'] = pd.to_numeric(scores['rmsd'], errors='coerce')

print(f"Running DSSR on {len(scores)} models...")
records = []
for i, row in scores.iterrows():
    if i % 50 == 0:
        print(f"  {i}/{len(scores)}")
    ss = dssr_ss(row['model_path'])
    if ss is None or len(ss) != N:
        continue
    r = pearsonr(unpaired_vec(ss), ref_vec)[0]
    records.append({'name': Path(row['model_path']).stem, 'source': row['source'],
                    'ss': ss, 'rmsd': row['rmsd'], 'corr': r})

df = pd.DataFrame(records)
print(f"Collected {len(df)} models with valid SS")

# --- collapse duplicate secondary structures ---
groups = []
for ss, g in df.groupby('ss'):
    g_sorted = g.sort_values('rmsd', na_position='last')
    rep = g_sorted.iloc[0]  # representative = best (lowest) RMSD model
    groups.append({
        'ss': ss,
        'corr': g['corr'].iloc[0],
        'n': len(g),
        'best_rmsd': g['rmsd'].min(),
        'rep_name': rep['name'],
    })
gdf = pd.DataFrame(groups).sort_values('corr', ascending=False).reset_index(drop=True)
print(f"{len(gdf)} unique secondary structures")

# --- build matrix: reference on top, then unique SS ---
labels = ['Mol9_reference_UtoG_buildloop  (target)']
matrix = [ref_vec]
ss_strings = [ref_ss]
for _, r in gdf.iterrows():
    matrix.append(unpaired_vec(r['ss']))
    ss_strings.append(r['ss'])
    rmsd_str = f"{r['best_rmsd']:.1f}A" if pd.notna(r['best_rmsd']) else "NA"
    extra = f" (+{r['n']-1} more)" if r['n'] > 1 else ""
    labels.append(f"{r['rep_name']}{extra}  r={r['corr']:.2f} rmsd={rmsd_str}")
matrix = np.array(matrix)

n_rows = len(matrix)

# --- font sizes (1.5x previous) ---
FS_CELL, FS_YLAB, FS_XLAB, FS_TITLE, FS_XTICK = 8, 11, 20, 22, 14

# --- plot: cell size tuned so dot-bracket chars are legible ---
cell_w, cell_h = 0.16, 0.34
fig, ax = plt.subplots(figsize=(max(14, N * cell_w + 7), max(7, n_rows * cell_h + 2)))
im = ax.imshow(matrix, aspect='auto', cmap=CMAP, vmin=0, vmax=0.5, interpolation='nearest')

# draw dot-bracket characters on every cell
for r in range(n_rows):
    ss = ss_strings[r]
    for c in range(N):
        ax.text(c, r, ss[c], ha='center', va='center', fontsize=FS_CELL,
                fontfamily='monospace', color='black')

ax.set_yticks(range(n_rows))
ax.set_yticklabels(labels, fontsize=FS_YLAB, fontfamily='monospace')
ax.set_xlabel('Nucleotide position', fontsize=FS_XLAB, fontweight='bold')
ax.set_title('Secondary-structure SHAPE proxy (white=paired, red=unpaired)\n'
             'unique structures sorted by correlation to reference',
             fontsize=FS_TITLE, fontweight='bold')

# x ticks at 10, 20, 30, ... (nucleotide positions divisible by 10)
xticks = [p - 1 for p in range(10, N + 1, 10)]  # 0-indexed columns
ax.set_xticks(xticks)
ax.set_xticklabels(range(10, N + 1, 10), fontsize=FS_XTICK)

# horizontal colorbar on the bottom
cbar = fig.colorbar(im, ax=ax, orientation='horizontal', fraction=0.04, pad=0.06,
                    ticks=[0, 0.5])
cbar.ax.set_xticklabels(['paired (0)', 'unpaired (0.5)'], fontsize=FS_XTICK)

plt.tight_layout()
out = 'results/shape_heatmap_naive.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
print(f"Saved {out}")
subprocess.run(['open', out])
