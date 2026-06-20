#!/usr/bin/env python3
"""
Main deliverable scatter plots. Produces TWO figures:

  results/shape_vs_rmsd_experimental.png
      SHAPE scores = correlation of predicted SHAPE to EXPERIMENTAL SHAPE.
      The reference (true) structure is a black star at its real (non-1.0)
      correlation to experiment — the ceiling / sanity check.

  results/shape_vs_rmsd_reference.png
      SHAPE scores = correlation of predicted SHAPE to the reference-MODEL
      simulated SHAPE (SGNM on the reference / reference paired-unpaired).
      The reference trivially sits at correlation = 1.0.

Each figure: 2x3 panels (naive & SGNM SHAPE vs RMSD and vs base-pair F1, plus
SGNM-vs-naive). Axis labels carry ↑/↓ arrows showing which direction is better.

Run after 01_compute_rmsd.py, 02_shape_naive.py, 03_shape_sgnm.py.
"""
import matplotlib
matplotlib.rcParams['font.family'] = 'Helvetica'
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import subprocess, sys
from scipy.stats import spearmanr, pearsonr

sys.path.insert(0, 'scripts')
from dssr_util import dssr_ss, unpaired_vec
from f1_score import compute_f1

REF_PDB = 'data/farfar2/Mol9_reference_UtoG_buildloop.pdb'
COLORS = {'farfar2': '#1f77b4', 'casp17': '#ff7f0e'}

df = pd.read_csv('results/scores.csv')
for c in df.columns:
    if c.startswith('shape_') or c in ('rmsd', 'f1_bp'):
        df[c] = pd.to_numeric(df[c], errors='coerce')

# --- base-pair F1 vs reference (compute if missing) ---
if 'f1_bp' not in df.columns or df['f1_bp'].notna().sum() == 0:
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

# --- reference self-scores (the black-star point in each figure) ---
exp = pd.read_csv('data/experimental_shape.csv')['shape'].values
prof = pd.read_csv('results/sgnm_profiles.csv', index_col='model')
ref_sgnm = prof.loc['reference'].values.astype(float)
ref_vec = unpaired_vec(dssr_ss(REF_PDB))
_m = ~np.isnan(exp)
REF = {
    'expt': {'rmsd': 0.0, 'f1_bp': 1.0,
             'naive': pearsonr(ref_vec[_m], exp[_m])[0],
             'sgnm': pearsonr(ref_sgnm[_m], exp[_m])[0]},
    'ref':  {'rmsd': 0.0, 'f1_bp': 1.0, 'naive': 1.0, 'sgnm': 1.0},
}


def panel(ax, xcol, ycol, xlabel, ylabel, title, refx, refy):
    for src in ['casp17', 'farfar2']:
        s = df[df['source'] == src]
        m = s[xcol].notna() & s[ycol].notna()
        ax.scatter(s[m][xcol], s[m][ycol], c=COLORS[src], label=src, s=40,
                   alpha=0.6, edgecolors='none')
    ax.scatter([refx], [refy], c='black', marker='*', s=340, edgecolors='white',
               linewidths=0.8, label='reference', zorder=5)
    ax.set_xlabel(xlabel, fontsize=12, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    m = df[xcol].notna() & df[ycol].notna()
    if m.sum() > 2:
        rho = spearmanr(df[m][xcol], df[m][ycol])[0]
        ax.text(0.97, 0.04, f'Spearman ρ={rho:.3f}', transform=ax.transAxes,
                ha='right', va='bottom', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))


def make_figure(kind, shape_label, out_png):
    """kind in {'expt','ref'}; uses the matching shape_*_vs_<kind> columns."""
    naive_col = f'shape_naive_vs_{kind}'
    sgnm_col = f'shape_sgnm_vs_{kind}'
    r = REF[kind]
    DOWN = r'  ($\downarrow$ better)'
    UP = r'  ($\uparrow$ better)'
    RMSD = 'RMSD to reference (Å)' + DOWN
    F1 = 'Base-pair F1 vs reference' + UP
    NAIVE = f'Naive SHAPE corr ({shape_label})' + UP
    SGNM = f'SGNM SHAPE corr ({shape_label})' + UP

    fig, ax = plt.subplots(2, 3, figsize=(20, 12))
    panel(ax[0, 0], 'rmsd', naive_col, RMSD, NAIVE, 'Naive SHAPE vs RMSD',
          r['rmsd'], r['naive'])
    panel(ax[0, 1], 'rmsd', sgnm_col, RMSD, SGNM, 'SGNM SHAPE vs RMSD',
          r['rmsd'], r['sgnm'])
    panel(ax[0, 2], 'rmsd', 'f1_bp', RMSD, F1, 'Base-pair F1 vs RMSD',
          r['rmsd'], r['f1_bp'])
    panel(ax[1, 0], 'f1_bp', naive_col, F1, NAIVE, 'Naive SHAPE vs F1',
          r['f1_bp'], r['naive'])
    panel(ax[1, 1], 'f1_bp', sgnm_col, F1, SGNM, 'SGNM SHAPE vs F1',
          r['f1_bp'], r['sgnm'])
    panel(ax[1, 2], naive_col, sgnm_col, NAIVE, SGNM, 'SGNM vs Naive SHAPE',
          r['naive'], r['sgnm'])

    title = ('SHAPE vs EXPERIMENTAL data' if kind == 'expt'
             else 'SHAPE vs REFERENCE-MODEL simulated data')
    fig.suptitle(f'{title} — RNA 3D decoy discrimination (CASP17 R2307)',
                 fontsize=18, fontweight='bold')
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches='tight')
    print(f"Saved {out_png}")
    subprocess.run(['open', out_png])


make_figure('expt', 'vs experiment', 'results/shape_vs_rmsd_experimental.png')
make_figure('ref', 'vs ref model', 'results/shape_vs_rmsd_reference.png')
