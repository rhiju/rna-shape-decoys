#!/usr/bin/env python3
"""
Create a heatmap of SHAPE predictions sorted by correlation to experimental.
Shows: experimental SHAPE at top, then naive SHAPE predictions
Rows are sorted by Pearson r to experimental.
Secondary structure (dot-bracket) displayed directly on heatmap as y-labels.
White (0) = paired, Light red (0.5) = unpaired.
"""
import matplotlib
matplotlib.rcParams['font.family'] = 'Helvetica'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from pathlib import Path
import sys
sys.path.insert(0, 'scripts')
from f1_score import get_base_pairs

# Custom colormap: white (0) to red (1)
cmap_white_red = LinearSegmentedColormap.from_list('white_red', [(1,1,1), (1,0,0)])

# Read scores
scores_df = pd.read_csv('results/scores.csv')

# Read experimental SHAPE (for now using reference SS as proxy)
exp_shape_df = pd.read_csv('data/experimental_shape.csv')
exp_shape_raw = exp_shape_df['shape'].values

# Convert experimental to 0/0.5 scale (0=paired, 0.5=unpaired)
# For now, use reference model structure
import tempfile, os
from pathlib import Path as PathlibPath

def get_ss_and_shape_naive(pdb_path):
    """Get secondary structure and naive SHAPE for a model."""
    import subprocess
    with tempfile.TemporaryDirectory() as tmpdir:
        abs_pdb = str(PathlibPath(pdb_path).resolve())
        old_cwd = os.getcwd()
        try:
            os.chdir(tmpdir)
            cmd = ['x3dna-dssr', f'-i={abs_pdb}']
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
            ss_file = PathlibPath('dssr-2ndstrs.dbn')
            if ss_file.exists():
                lines = ss_file.read_text().splitlines()
                if len(lines) >= 3:
                    ss = lines[2].strip()
                    # Convert SS to 0/0.5 scale (0=paired, 0.5=unpaired)
                    naive_shape = np.array([0.5 if c == '.' else 0.0 for c in ss])
                    return ss, naive_shape
        finally:
            os.chdir(old_cwd)
    return None, None

# Get reference model
ref_ss, ref_shape = get_ss_and_shape_naive('data/farfar2/Mol9_reference_UtoG_buildloop.pdb')
print(f"Reference SS: {ref_ss}")
print(f"Reference shape (0/0.5 scale): {ref_shape}")

# Use reference model SHAPE as experimental proxy
exp_shape = ref_shape
exp_ss = ref_ss

# Prepare data for heatmap
models_list = []
ss_list = []
naive_shapes = []

for idx, row in scores_df.iterrows():
    if idx % 50 == 0:
        print(f"Processing {idx}/{len(scores_df)}...")

    model_path = row['model_path']
    ss, naive_shape = get_ss_and_shape_naive(model_path)

    if naive_shape is not None and len(naive_shape) == len(exp_shape):
        models_list.append({
            'path': model_path,
            'name': PathlibPath(model_path).stem,
            'source': row['source']
        })
        ss_list.append(ss)
        naive_shapes.append(naive_shape)

print(f"Collected {len(models_list)} models with valid SHAPE data")

# Compute correlations to experimental SHAPE
from scipy.stats import pearsonr

correlations = []
for i, naive_shape in enumerate(naive_shapes):
    try:
        r, _ = pearsonr(naive_shape, exp_shape)
    except:
        r = np.nan
    correlations.append(r)

# Sort by correlation (descending)
sorted_indices = np.argsort([-c if not np.isnan(c) else -np.inf for c in correlations])

# Create heatmap data: experimental at top, then models
heatmap_data = np.vstack([
    exp_shape.reshape(1, -1),
    np.array(naive_shapes)[sorted_indices]
])

# Create y-labels with secondary structure and correlation
y_labels = [f"Experimental: {exp_ss}"]
for i in sorted_indices:
    corr_val = correlations[i]
    if np.isnan(corr_val):
        corr_str = "r=NA"
    else:
        corr_str = f"r={corr_val:.3f}"
    y_labels.append(f"{models_list[i]['name']} [{ss_list[i]}] {corr_str}")

# Create figure
fig, ax = plt.subplots(figsize=(20, max(10, len(sorted_indices) * 0.25)))

# Plot heatmap
im = ax.imshow(heatmap_data, aspect='auto', cmap=cmap_white_red, vmin=0, vmax=0.5, interpolation='nearest')

# Set ticks and labels
ax.set_yticks(range(len(y_labels)))
ax.set_yticklabels(y_labels, fontsize=7, fontfamily='monospace')
ax.set_xlabel('Nucleotide Position', fontsize=14, fontweight='bold')
ax.set_title('SHAPE Profiles: Naive DSSR Predictions Sorted by Correlation to Reference', fontsize=16, fontweight='bold')

# Colorbar
cbar = plt.colorbar(im, ax=ax, label='SHAPE Value (0=paired, 0.5=unpaired)')

# X-axis ticks at regular intervals
n_nt = len(exp_shape)
tick_interval = max(10, n_nt // 20)
ax.set_xticks(range(0, n_nt, tick_interval))
ax.set_xticklabels(range(1, n_nt + 1, tick_interval), fontsize=10)

plt.tight_layout()
out_png = 'results/shape_heatmap_naive.png'
fig.savefig(out_png, dpi=150, bbox_inches='tight')
print(f"Saved: {out_png}")

# Open in Preview
import subprocess
subprocess.run(['open', out_png])

# Save correlation data for reference
corr_df = pd.DataFrame({
    'model': [m['name'] for m in models_list],
    'source': [m['source'] for m in models_list],
    'secondary_structure': ss_list,
    'pearson_r_naive': correlations
})
corr_df = corr_df.sort_values('pearson_r_naive', ascending=False, na_position='last')
corr_df.to_csv('results/shape_correlations_naive.csv', index=False)
print(f"Saved correlation data to results/shape_correlations_naive.csv")

print(f"\nHeatmap includes {len(models_list)} models sorted by Pearson r to reference")

