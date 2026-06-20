#!/usr/bin/env python3
"""
Predict SHAPE using SGNM model on all structures.
Requires: gnm-checkpoint.pth to be downloaded from https://github.com/hmblair/sgnm/releases
"""
import sys
import torch
import pandas as pd
from pathlib import Path
import numpy as np
from scipy.stats import pearsonr

# Import SGNM directly (avoid equivariant import issues for now)
from sgnm.models import SGNM
import ciffy

# Paths
model_checkpoint = Path('/Users/rhiju/src/sgnm/gnm-checkpoint.pth')
results_dir = Path('results')
scores_file = results_dir / 'scores.csv'

print("SGNM SHAPE Prediction")
print("=" * 50)

# Check for model checkpoint
if not model_checkpoint.exists():
    print(f"\n❌ ERROR: Model checkpoint not found at {model_checkpoint}")
    print("\nTo download the checkpoint:")
    print("  cd ~/src/sgnm")
    print("  wget https://github.com/hmblair/sgnm/releases/download/v1.0.0/gnm-checkpoint.pth")
    print("\nOr download manually from:")
    print("  https://github.com/hmblair/sgnm/releases")
    sys.exit(1)

print(f"✓ Found checkpoint: {model_checkpoint}")

# Load model
print("\nLoading SGNM model...")
try:
    model = SGNM.load(str(model_checkpoint))
    model.eval()
    print("✓ Model loaded successfully")
except Exception as e:
    print(f"❌ Failed to load model: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Read scores
scores_df = pd.read_csv(scores_file)

# Read experimental SHAPE for correlation
exp_shape_df = pd.read_csv('data/experimental_shape.csv')
exp_shape = exp_shape_df['shape'].values

# Reference model SHAPE (for ceiling analysis)
print("\nPredicting SHAPE for reference model...")
ref_pdb = 'data/farfar2/Mol9_reference_UtoG_buildloop.pdb'
try:
    ref_poly = ciffy.load(ref_pdb, backend='torch').strip()
    with torch.no_grad():
        ref_reactivity = model.ciffy(ref_poly)  # (N, 2) for [SHAPE, DMS]
    ref_shape = ref_reactivity[:, 0].cpu().numpy()
    print(f"✓ Reference SHAPE shape: {ref_shape.shape}")
except Exception as e:
    print(f"⚠ Could not load reference with ciffy: {e}")
    ref_shape = None

# Predict SHAPE for all models
print(f"\nPredicting SHAPE for {len(scores_df)} models...")

shape_sgnm_vs_expt = []
shape_sgnm_vs_ref = []

for idx, row in scores_df.iterrows():
    if idx % 50 == 0:
        print(f"  {idx}/{len(scores_df)}")

    model_path = row['model_path']

    try:
        # Load structure
        poly = ciffy.load(model_path, backend='torch').strip()

        # Predict SHAPE
        with torch.no_grad():
            reactivity = model.ciffy(poly)  # (N, 2)
        pred_shape = reactivity[:, 0].cpu().numpy()

        # Correlate with experimental
        if len(pred_shape) == len(exp_shape):
            try:
                r_expt, _ = pearsonr(pred_shape, exp_shape)
            except:
                r_expt = np.nan
        else:
            r_expt = np.nan

        # Correlate with reference model SHAPE
        if ref_shape is not None and len(pred_shape) == len(ref_shape):
            try:
                r_ref, _ = pearsonr(pred_shape, ref_shape)
            except:
                r_ref = np.nan
        else:
            r_ref = np.nan

    except Exception as e:
        print(f"    ⚠ Error processing {Path(model_path).name}: {e}")
        r_expt = np.nan
        r_ref = np.nan

    shape_sgnm_vs_expt.append(r_expt)
    shape_sgnm_vs_ref.append(r_ref)

print(f"✓ Predictions complete")

# Add to scores
scores_df['shape_sgnm_vs_expt'] = shape_sgnm_vs_expt
scores_df['shape_sgnm_vs_ref'] = shape_sgnm_vs_ref

# Save
scores_df.to_csv(scores_file, index=False)
print(f"\n✓ Updated {scores_file} with SGNM predictions")

# Summary
print("\nSummary:")
valid_expt = [r for r in shape_sgnm_vs_expt if not np.isnan(r)]
valid_ref = [r for r in shape_sgnm_vs_ref if not np.isnan(r)]

if valid_expt:
    print(f"  SGNM vs experimental: mean r={np.mean(valid_expt):.3f}, std={np.std(valid_expt):.3f}")
if valid_ref:
    print(f"  SGNM vs reference:   mean r={np.mean(valid_ref):.3f}, std={np.std(valid_ref):.3f}")

