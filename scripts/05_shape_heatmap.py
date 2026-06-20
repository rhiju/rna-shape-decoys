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
import subprocess, tempfile, os
from scipy.stats import pearsonr

CMAP = LinearSegmentedColormap.from_list('white_red', [(1, 1, 1), (1, 0, 0)])
REF_PDB = 'data/farfar2/Mol9_reference_UtoG_buildloop.pdb'


def dssr_ss(pdb_path):
    """Return dot-bracket secondary structure for a PDB via DSSR."""
    with tempfile.TemporaryDirectory() as tmp:
        ap = os.path.abspath(pdb_path)
        cwd = os.getcwd()
        try:
            os.chdir(tmp)
            subprocess.run(['x3dna-dssr', f'-i={ap}'],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            f = Path('dssr-2ndstrs.dbn')
            if f.exists():
                lines = f.read_text().splitlines()
                if len(lines) >= 3:
                    return lines[2].strip()
        finally:
            os.chdir(cwd)
    return None


def unpaired_vec(ss):
    """0.5 for unpaired (.), 0.0 for paired (everything else)."""
    return np.array([0.5 if c == '.' else 0.0 for c in ss])


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
    groups.append({
        'ss': ss,
        'corr': g['corr'].iloc[0],
        'n': len(g),
        'best_rmsd': g['rmsd'].min(),
        'sources': '/'.join(sorted(g['source'].unique())),
    })
gdf = pd.DataFrame(groups).sort_values('corr', ascending=False).reset_index(drop=True)
print(f"{len(gdf)} unique secondary structures")

# --- build matrix: experimental/reference on top, then unique SS ---
labels = [f'REFERENCE  (target)']
matrix = [ref_vec]
ss_strings = [ref_ss]
for _, r in gdf.iterrows():
    matrix.append(unpaired_vec(r['ss']))
    ss_strings.append(r['ss'])
    rmsd_str = f"{r['best_rmsd']:.1f}A" if pd.notna(r['best_rmsd']) else "NA"
    labels.append(f"r={r['corr']:.2f}  n={r['n']:<2d} best={rmsd_str} [{r['sources']}]")
matrix = np.array(matrix)

n_rows = len(matrix)

# --- plot: cell size tuned so dot-bracket chars are legible ---
cell_w, cell_h = 0.13, 0.30
fig, ax = plt.subplots(figsize=(max(12, N * cell_w + 5), max(6, n_rows * cell_h + 1)))
ax.imshow(matrix, aspect='auto', cmap=CMAP, vmin=0, vmax=0.5, interpolation='nearest')

# draw dot-bracket characters on every cell
for r in range(n_rows):
    ss = ss_strings[r]
    for c in range(N):
        ax.text(c, r, ss[c], ha='center', va='center', fontsize=5,
                fontfamily='monospace', color='black')

ax.set_yticks(range(n_rows))
ax.set_yticklabels(labels, fontsize=7, fontfamily='monospace')
ax.set_xlabel('Nucleotide position', fontsize=13, fontweight='bold')
ax.set_title('Secondary-structure SHAPE proxy (white=paired, red=unpaired)\n'
             'unique structures sorted by correlation to reference',
             fontsize=15, fontweight='bold')

# x ticks every 10 nt
ax.set_xticks(range(0, N, 10))
ax.set_xticklabels(range(1, N + 1, 10), fontsize=9)

plt.tight_layout()
out = 'results/shape_heatmap_naive.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
print(f"Saved {out}")
subprocess.run(['open', out])
