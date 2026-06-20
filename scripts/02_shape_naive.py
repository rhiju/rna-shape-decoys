#!/usr/bin/env python3
"""
Naive SHAPE scoring: use binary paired/unpaired from DSSR as SHAPE proxy.
Compute correlation with experimental SHAPE (reference model predictions).
"""
import subprocess
import tempfile
import csv
import os
from pathlib import Path
import numpy as np
from scipy.stats import pearsonr
import sys

def run_dssr(pdb_file):
    """Run DSSR and return dot-bracket notation."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Get absolute path BEFORE changing directory
        abs_pdb = str(Path(pdb_file).resolve())
        old_cwd = os.getcwd()
        try:
            os.chdir(tmpdir)
            cmd = ['x3dna-dssr', f'-i={abs_pdb}']
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
            ss_file = Path('dssr-2ndstrs.dbn')
            if ss_file.exists():
                lines = ss_file.read_text().splitlines()
                if len(lines) >= 3:
                    return lines[2].strip()
        finally:
            os.chdir(old_cwd)
    return None

def ss_to_binary(dot_bracket):
    """Convert dot-bracket to binary: 0=paired, 1=unpaired."""
    binary = []
    for c in dot_bracket:
        if c == '.':
            binary.append(1)
        else:
            binary.append(0)
    return np.array(binary)

# Paths
results_dir = Path('results')
scores_file = results_dir / 'scores.csv'

# Read existing scores
scores = {}
if scores_file.exists():
    with open(scores_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            scores[row['model_path']] = dict(row)

# For now, create a synthetic SHAPE profile from the reference structure
# (We'll replace this with actual SGNM/ERM predictions in the next phase)
ref_pdb = 'data/farfar2/Mol9_reference_UtoG_buildloop.pdb'
ref_ss = run_dssr(ref_pdb)

if ref_ss:
    ref_binary = ss_to_binary(ref_ss)
    print(f"Reference SS: {ref_ss}")
    print(f"Paired/unpaired profile: {ref_binary}")

    # Use reference binary as "experimental SHAPE" for now
    # (This is temporary; will be replaced with actual SHAPE)
    exp_shape = ref_binary.astype(float)
else:
    print("Could not get reference secondary structure")
    sys.exit(1)

# Compute naive SHAPE scores for all models
model_files = []
for pdb in Path('data/farfar2').glob('Mol9.out.*.pdb'):
    model_files.append(str(pdb))
for pdb in Path('data/casp17/R2307').glob('R2307TS*.pdb'):
    model_files.append(str(pdb))

print(f"Computing naive SHAPE scores for {len(model_files)} models...")

for i, model_path in enumerate(model_files):
    if i % 50 == 0:
        print(f"  {i}/{len(model_files)}")

    # Get secondary structure
    pred_ss = run_dssr(model_path)
    if not pred_ss:
        continue

    pred_binary = ss_to_binary(pred_ss)

    # Compute correlation
    if len(pred_binary) == len(exp_shape):
        try:
            r, p = pearsonr(pred_binary, exp_shape)
        except:
            r = np.nan
    else:
        r = np.nan

    # Update score
    if model_path in scores:
        scores[model_path]['shape_naive'] = r
    else:
        scores[model_path] = {
            'model_path': model_path,
            'source': 'farfar2' if 'farfar2' in model_path else 'casp17',
            'group': Path(model_path).stem,
            'rmsd': None,
            'f1_bp': None,
            'shape_naive': r
        }

# Write updated scores
fieldnames = ['model_path', 'source', 'group', 'rmsd', 'f1_bp', 'shape_naive']
with open(scores_file, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, restval='')
    w.writeheader()
    for path in sorted(scores.keys()):
        w.writerow({k: scores[path].get(k, '') for k in fieldnames})

print(f"Saved: {scores_file}")
print(f"Total scores: {len(scores)}")
