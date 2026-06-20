#!/usr/bin/env python3
"""
Compute RMSD and base-pair F1 scores for all models vs reference.
"""
import subprocess
import csv
import re
import os
from pathlib import Path
from f1_score import compute_f1

def run_dssr(pdb_file):
    """Run DSSR and return dot-bracket notation."""
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        old_cwd = Path.cwd()
        try:
            import os
            os.chdir(tmpdir)
            cmd = f"x3dna-dssr -i={pdb_file} > /dev/null 2> /dev/null"
            subprocess.run(cmd, shell=True, check=False)
            ss_file = Path('dssr-2ndstrs.dbn')
            if ss_file.exists():
                lines = ss_file.read_text().splitlines()
                if len(lines) >= 3:
                    return lines[2].strip()
        finally:
            os.chdir(old_cwd)
    return None

def parse_tmscore_rmsd(tmscore_output):
    """Extract RMSD from TMscore output."""
    for line in tmscore_output.splitlines():
        if line.strip().startswith('RMSD'):
            m = re.search(r'([\d.]+)', line)
            if m:
                return float(m.group(1))
    return None

# Paths
ref_pdb = Path('data/farfar2/Mol9_reference_UtoG_buildloop.pdb')
farfar2_dir = Path('data/farfar2')
casp17_dir = Path('data/casp17/R2307')
results_dir = Path('results')

# Get reference SS
ref_ss = run_dssr(str(ref_pdb))
print(f"Reference SS: {ref_ss}")

# Build model list
models = []
for pdb in sorted(farfar2_dir.glob('Mol9.out.*.pdb')):
    models.append({
        'path': str(pdb),
        'source': 'farfar2',
        'group': 'FARFAR2'
    })

for pdb in sorted(casp17_dir.glob('R2307TS*.pdb')):
    group = pdb.stem[:8]  # Extract group ID like "R2307TS017"
    models.append({
        'path': str(pdb),
        'source': 'casp17',
        'group': group
    })

print(f"Total models: {len(models)}")

# Compute RMSD for each
results = []
for i, model in enumerate(models):
    if i % 20 == 0:
        print(f"Processing {i}/{len(models)}...")

    # RMSD
    cmd = ['TMscore', model['path'], str(ref_pdb), '-atom', " C1'"]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=False)
        rmsd = parse_tmscore_rmsd(res.stdout)
    except Exception as e:
        print(f"Error computing RMSD for {model['path']}: {e}")
        rmsd = None

    # F1 of base pairs
    pred_ss = run_dssr(model['path'])
    if pred_ss and ref_ss:
        p, r, f1 = compute_f1(ref_ss, pred_ss)
    else:
        f1 = None

    results.append({
        'model_path': model['path'],
        'source': model['source'],
        'group': model['group'],
        'rmsd': rmsd,
        'f1_bp': f1
    })

# Write results
output_csv = results_dir / 'scores.csv'
with open(output_csv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['model_path', 'source', 'group', 'rmsd', 'f1_bp'])
    w.writeheader()
    for r in results:
        w.writerow(r)

print(f"\nSaved: {output_csv}")

# Summary
casp17_rmsd = [r['rmsd'] for r in results if r['source'] == 'casp17' and r['rmsd']]
farfar2_rmsd = [r['rmsd'] for r in results if r['source'] == 'farfar2' and r['rmsd']]

if casp17_rmsd:
    print(f"CASP17: mean RMSD = {sum(casp17_rmsd)/len(casp17_rmsd):.2f}, min = {min(casp17_rmsd):.2f}, max = {max(casp17_rmsd):.2f}")
if farfar2_rmsd:
    print(f"FARFAR2: mean RMSD = {sum(farfar2_rmsd)/len(farfar2_rmsd):.2f}, min = {min(farfar2_rmsd):.2f}, max = {max(farfar2_rmsd):.2f}")
