#!/usr/bin/env python3
"""
Create a heatmap of SHAPE predictions sorted by correlation to experimental.
Shows: experimental SHAPE, naive SHAPE, SGNM SHAPE, ERM SHAPE (once available)
Rows are sorted by Pearson r to experimental.
Secondary structure (dot-bracket) shown below each profile.
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

# Read experimental SHAPE (placeholder for now - will be actual experimental once we have it)
exp_shape_df = pd.read_csv('data/experimental_shape.csv')
exp_shape = exp_shape_df['shape'].values

# Get reference secondary structure for pairing info
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
                    # Convert SS to binary (0=paired, 1=unpaired)
                    naive_shape = np.array([1.0 if c == '.' else 0.0 for c in ss])
                    return ss, naive_shape
        finally:
            os.chdir(old_cwd)
    return None, None

# Collect all SHAPE data
ref_ss, ref_naive_shape = get_ss_and_shape_naive('data/farfar2/Mol9_reference_UtoG_buildloop.pdb')

print(f"Reference SS: {ref_ss}")
print(f"Reference naive SHAPE shape: {ref_naive_shape.shape}")

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

# Compute correlations to experimental (naive SHAPE)
from scipy.stats import pearsonr

correlations = []
for i, naive_shape in enumerate(naive_shapes):
    try:
        r, _ = pearsonr(naive_shape, exp_shape)
    except:
        r = np.nan
    correlations.append(r)

# Sort by correlation
sorted_indices = np.argsort(correlations)[::-1]  # Descending

# Create heatmap
n_models = len(sorted_indices)
n_nt = len(exp_shape)

# Stack SHAPE profiles (rows = models, cols = nucleotides)
# Include experimental at the top
heatmap_data = np.vstack([
    exp_shape.reshape(1, -1),
    np.array(naive_shapes)[sorted_indices]
])

# Create figure with multiple rows
fig = plt.figure(figsize=(16, max(8, n_models * 0.3)))

# Main heatmap
ax_hm = plt.subplot(2, 1, 1)
im = ax_hm.imshow(heatmap_data, aspect='auto', cmap=cmap_white_red, vmin=0, vmax=1)

# Y labels: "Experimental" + model names sorted by correlation
y_labels = ['Experimental'] + [
    f"{models_list[i]['name']} (r={correlations[i]:.2f})"
    for i in sorted_indices
]

ax_hm.set_yticks(range(len(y_labels)))
ax_hm.set_yticklabels(y_labels, fontsize=8)
ax_hm.set_xlabel('Nucleotide Position', fontsize=12, fontweight='bold')
ax_hm.set_title('SHAPE Profiles: Naive DSSR Predictions Sorted by Correlation', fontsize=14, fontweight='bold')

# Colorbar
cbar = plt.colorbar(im, ax=ax_hm, label='SHAPE Value (0=paired, 1=unpaired)')

# Secondary structures subplot
ax_ss = plt.subplot(2, 1, 2)
ax_ss.axis('off')

# Display reference SS and all model SS
ss_text = "Secondary Structures (dot-bracket notation):\n\n"
ss_text += f"Reference: {ref_ss}\n\n"
for i, idx in enumerate(sorted_indices[:min(10, n_models)]):  # Show top 10
    ss_text += f"{models_list[idx]['name']}: {ss_list[idx]}\n"

if n_models > 10:
    ss_text += f"\n... and {n_models - 10} more structures"

ax_ss.text(0.05, 0.95, ss_text, transform=ax_ss.transAxes,
          fontfamily='monospace', fontsize=9, verticalalignment='top',
          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

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
    'pearson_r_naive': correlations
})
corr_df = corr_df.sort_values('pearson_r_naive', ascending=False)
corr_df.to_csv('results/shape_correlations_naive.csv', index=False)
print(f"Saved correlation data to results/shape_correlations_naive.csv")

