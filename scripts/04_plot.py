#!/usr/bin/env python3
"""
Plot SHAPE match vs RMSD and F1 metrics.
"""
import matplotlib
matplotlib.rcParams['font.family'] = 'Helvetica'
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import spearmanr

# Read scores
df = pd.read_csv('results/scores.csv')
df['rmsd'] = pd.to_numeric(df['rmsd'], errors='coerce')
df['shape_naive'] = pd.to_numeric(df['shape_naive'], errors='coerce')

# Create plot
fig, axes = plt.subplots(1, 1, figsize=(10, 6))

colors = {'farfar2': '#1f77b4', 'casp17': '#ff7f0e'}
for source in ['casp17', 'farfar2']:
    subset = df[df['source'] == source]
    subset_valid = subset[subset['rmsd'].notna() & subset['shape_naive'].notna()]

    axes.scatter(subset_valid['rmsd'], subset_valid['shape_naive'],
               c=colors[source], label=source, s=50, alpha=0.6, edgecolors='none')

axes.set_xlabel('RMSD to reference (Å)', fontsize=14, fontweight='bold')
axes.set_ylabel('Naive SHAPE correlation (Pearson r)', fontsize=14, fontweight='bold')
axes.set_title('SHAPE Score vs RMSD', fontsize=18, fontweight='bold')
axes.legend(fontsize=12)
axes.grid(True, alpha=0.3)

# Compute correlation
all_valid = df[df['rmsd'].notna() & df['shape_naive'].notna()]
if len(all_valid) > 1:
    rho, pval = spearmanr(all_valid['rmsd'], all_valid['shape_naive'])
    axes.text(0.95, 0.05, f'Spearman ρ={rho:.3f}\np={pval:.2e}',
            transform=axes.transAxes, ha='right', va='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
            fontsize=11)

plt.tight_layout()
out_png = 'results/shape_vs_rmsd.png'
fig.savefig(out_png, dpi=150, bbox_inches='tight')
print(f"Saved: {out_png}")

# Open in Preview
import subprocess
subprocess.run(['open', out_png])
