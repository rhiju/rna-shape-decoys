#!/usr/bin/env python3
"""Shared DSSR helper. Returns dot-bracket secondary structure for a PDB.

DSSR inserts '&' at backbone breaks (chain discontinuities); these are stripped
so the returned string aligns 1:1 with the nucleotides.
"""
import subprocess
import tempfile
import os
from pathlib import Path
import numpy as np


def dssr_ss(pdb_path):
    """Return dot-bracket secondary structure (no '&' separators), or None."""
    with tempfile.TemporaryDirectory() as tmp:
        ap = os.path.abspath(pdb_path)
        cwd = os.getcwd()
        try:
            os.chdir(tmp)
            subprocess.run(['x3dna-dssr', f'-i={ap}'],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            f = Path('dssr-2ndstrs.dbn')
            if f.exists():
                lines = f.read_text().splitlines()
                if len(lines) >= 3:
                    return lines[2].strip().replace('&', '')
        finally:
            os.chdir(cwd)
    return None


def unpaired_vec(ss):
    """0.5 for unpaired (.), 0.0 for paired (everything else)."""
    return np.array([0.5 if c == '.' else 0.0 for c in ss])
