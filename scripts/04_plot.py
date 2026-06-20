#!/usr/bin/env python3
"""
Main deliverable plots: SHAPE-score vs RMSD (and vs base-pair F1).

Computes base-pair F1 (vs reference DSSR structure) for every model, then makes
a 2x2 panel:
  (0,0) naive paired/unpaired SHAPE corr vs RMSD
  (0,1) SGNM SHAPE corr (vs reference profile) vs RMSD
  (1,0) base-pair F1 vs RMSD
  (1,1) SGNM SHAPE corr vs base-pair F1

Run after 01_compute_rmsd.py and 03_shape_sgnm.py.
"""
import matplotlib
matplotlib.rcParams['font.family'] = 'Helvetica'
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import subprocess, tempfile, os, sys
from scipy.stats import spearmanr, pearsonr

sys.path.insert(0, 'scripts')
from f1_score import compute_f1

REF_PDB = 'data/farfar2/Mol9_reference_UtoG_buildloop.pdb'
COLORS = {'farfar2': '#1f77b4', 'casp17': '#ff7f0e'}


def dssr_ss(pdb_path):
    with tempfile.TemporaryDirectory() as tmp:
        ap = os.path.abspath(pdb_path); cwd = os.getcwd()
        try:
            os.chdir(tmp)
            subprocess.run(['x3dna-dssr', f'-i={ap}'],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            f = Path('dssr-2ndstrs.dbn')
            if f.exists():
                ls = f.read_text().splitlines()
                if len(ls) >= 3:
                    return ls[2].strip()
        finally:
            os.chdir(cwd)
    return None


df = pd.read_csv('results/scores.csv')
df['rmsd'] = pd.to_numeric(df['rmsd'], errors='coerce')
df['shape_naive'] = pd.to_numeric(df['shape_naive'], errors='coerce')
df['shape_sgnm_vs_ref'] = pd.to_numeric(df['shape_sgnm_vs_ref'], errors='coerce')

# --- compute base-pair F1 vs reference ---
ref_ss = dssr_ss(REF_PDB)
print(f"Computing base-pair F1 for {len(df)} models...")
f1s = []
for i, row in df.iterrows():
    if i % 50 == 0:
        print(f"  {i}/{len(df)}")
    ss = dssr_ss(row['model_path'])
    f1s.append(compute_f1(ref_ss, ss)[2] if ss else np.nan)
df['f1_bp'] = f1s
df.to_csv('results/scores.csv', index=False)
print("Updated results/scores.csv with f1_bp")


def panel(ax, x, y, xlabel, ylabel, title):
    for src in ['casp17', 'farfar2']:
        s = df[df['source'] == src]
        m = s[x].notna() & s[y].notna()
        ax.scatter(s[m][x], s[m][y], c=COLORS[src], label=src, s=40, alpha=0.6,
                   edgecolors='none')
    ax.set_xlabel(xlabel, fontsize=12, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    m = df[x].notna() & df[y].notna()
    if m.sum() > 2:
        rho = spearmanr(df[m][x], df[m][y])[0]
        ax.text(0.97, 0.04, f'Spearman ρ={rho:.3f}', transform=ax.transAxes,
                ha='right', va='bottom', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))


fig, axes = plt.subplots(2, 2, figsize=(15, 12))
panel(axes[0, 0], 'rmsd', 'shape_naive',
      'RMSD to reference (Å)', 'Naive SHAPE corr (Pearson r)',
      'Naive paired/unpaired SHAPE vs RMSD')
panel(axes[0, 1], 'rmsd', 'shape_sgnm_vs_ref',
      'RMSD to reference (Å)', 'SGNM SHAPE corr to ref (Pearson r)',
      'SGNM-predicted SHAPE vs RMSD')
panel(axes[1, 0], 'rmsd', 'f1_bp',
      'RMSD to reference (Å)', 'Base-pair F1 vs reference',
      'Secondary-structure F1 vs RMSD')
panel(axes[1, 1], 'f1_bp', 'shape_sgnm_vs_ref',
      'Base-pair F1 vs reference', 'SGNM SHAPE corr to ref (Pearson r)',
      'SGNM SHAPE vs base-pair F1')

fig.suptitle('SHAPE-based discrimination of RNA 3D decoys (CASP17 R2307)',
             fontsize=17, fontweight='bold')
fig.tight_layout()
out = 'results/shape_vs_rmsd.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
print(f"Saved {out}")
subprocess.run(['open', out])
