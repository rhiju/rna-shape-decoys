#!/usr/bin/env python3
"""
Convert CASP17 TS format to PDB format.
The TS format is nearly identical to PDB, just with different header lines.
"""
import sys
import os
from pathlib import Path

def convert_ts_to_pdb(ts_file, pdb_file):
    """Convert a single TS file to PDB format."""
    with open(ts_file, 'r') as f:
        lines = f.readlines()

    pdb_lines = []
    for line in lines:
        if line.startswith('PFRMAT') or line.startswith('TARGET') or line.startswith('PARENT') or line.startswith('MODEL'):
            continue
        if line.startswith('ATOM'):
            pdb_lines.append(line)
        elif line.startswith('END'):
            pdb_lines.append('END\n')
            break

    with open(pdb_file, 'w') as f:
        f.writelines(pdb_lines)

# Convert all TS files in casp17/R2307/ to PDB
casp17_dir = Path('data/casp17/R2307')
if casp17_dir.exists():
    for ts_file in sorted(casp17_dir.glob('*')):
        if ts_file.is_file() and not ts_file.suffix:  # No suffix = TS format
            pdb_file = ts_file.with_suffix('.pdb')
            print(f"Converting {ts_file.name} -> {pdb_file.name}")
            convert_ts_to_pdb(ts_file, pdb_file)
    print(f"Converted {len(list(casp17_dir.glob('*.pdb')))} files")
else:
    print(f"Directory not found: {casp17_dir}")
