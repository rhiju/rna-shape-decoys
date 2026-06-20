#!/usr/bin/env python3
"""
Run the SGNM Equivariant Reactivity Model (ERM) on all structures.

ERM is atom-level (ciffy PolymerEmbedding + k-NN) and needs flash-eq, which is
CUDA-only — run this on a GPU node / Colab, NOT on the Mac. It loads each PDB by
converting to mmCIF (scripts/pdb_to_cif.py, column-aligned so ciffy types atoms
correctly) and calling model.ciffy(poly).

Writes results/erm_profiles.csv (model, pos_1..pos_N) — same format as
results/sgnm_profiles.csv. Afterwards run scripts/ingest_erm.py locally.
"""
import sys
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, 'scripts')
import torch
import ciffy
import sgnm
from pdb_to_cif import pdb_to_cif

CKPT = sys.argv[1] if len(sys.argv) > 1 else 'equivariant-checkpoint.pth'

print(f"Loading ERM from {CKPT} ...")
model = sgnm.load(CKPT)
model.eval()
dev = 'cuda' if torch.cuda.is_available() else 'cpu'
model = model.to(dev)
print(f"  device: {dev}")


def erm_shape(pdb_path):
    with tempfile.NamedTemporaryFile(suffix='.cif') as tf:
        pdb_to_cif(pdb_path, tf.name)
        poly = ciffy.load(tf.name, backend='torch').strip()
    if dev == 'cuda':
        poly = poly.to('cuda')
    with torch.no_grad():
        pred = model.ciffy(poly)
    return pred[:, 0].cpu().numpy()


ref_pdb = 'data/farfar2/Mol9_reference_UtoG_buildloop.pdb'
model_files = [ref_pdb]
model_files += sorted(str(p) for p in Path('data/farfar2').glob('Mol9.out.*.pdb'))
model_files += sorted(str(p) for p in Path('data/casp17/R2307').glob('R2307TS*.pdb'))

profiles = {}
for i, pdb in enumerate(model_files):
    if i % 20 == 0:
        print(f"  {i}/{len(model_files)}")
    name = 'reference' if pdb == ref_pdb else Path(pdb).stem
    try:
        profiles[name] = erm_shape(pdb)
    except Exception as e:
        print(f"    skip {name}: {e}")

N = len(profiles['reference'])
prof = pd.DataFrame({k: v for k, v in profiles.items() if len(v) == N}).T
prof.columns = [f'pos_{i+1}' for i in range(N)]
prof.index.name = 'model'
Path('results').mkdir(exist_ok=True)
prof.to_csv('results/erm_profiles.csv')
print(f"Saved results/erm_profiles.csv ({prof.shape[0]} profiles x {N} positions)")
