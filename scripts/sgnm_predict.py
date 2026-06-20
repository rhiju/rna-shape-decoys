#!/usr/bin/env python3
"""
Predict SHAPE/DMS reactivity from an RNA PDB using the SGNM (GNM) model,
bypassing ciffy's mmCIF file loader.

SGNM's forward(coords, frames) needs only:
  - coords:  per-residue heavy-atom centroid, shape (N, 3)
  - frames:  per-residue C2/C4/C6 nucleobase frame, shape (N, 3, 3)

These are computed directly from the PDB here. This was verified to
reproduce `model.ciffy(poly)` exactly (max abs diff 0.0) on the SGNM
example structures, so it is faithful to the model's intended input.

Hydrogens are stripped so FARFAR2 (H-containing) and CASP17 (no H)
structures get a consistent heavy-atom representation.
"""
import torch
from ciffy.geometry.transforms import frame_from_positions


def parse_pdb_residues(pdb_path):
    """Parse first model of a PDB into ordered residues of heavy atoms.

    Returns a list of dicts: {chain, res_seq, res_name, atoms: {name: (x,y,z)}}
    """
    residues = []
    index = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('ENDMDL'):
                break
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            atom_name = line[12:16].strip()
            elem = line[76:78].strip() or atom_name[0]
            if elem == 'H' or (atom_name and atom_name[0] == 'H'):
                continue  # skip hydrogens
            res_name = line[17:20].strip()
            chain = line[21].strip() or 'A'
            res_seq = int(line[22:26])
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            key = (chain, res_seq)
            if key not in index:
                index[key] = len(residues)
                residues.append({'chain': chain, 'res_seq': res_seq,
                                 'res_name': res_name, 'atoms': {}})
            residues[index[key]]['atoms'][atom_name] = (x, y, z)
    return residues


def coords_and_frames(pdb_path):
    """Compute per-residue centers and C2/C4/C6 frame positions from a PDB.

    Returns (centers (N,3), positions (N,3,3), kept_residue_indices).
    Residues missing any of C2/C4/C6 are skipped (with their index recorded).
    """
    residues = parse_pdb_residues(pdb_path)
    centers = []
    positions = []
    kept = []
    for i, res in enumerate(residues):
        atoms = res['atoms']
        if not all(a in atoms for a in ('C2', 'C4', 'C6')):
            continue
        all_xyz = torch.tensor(list(atoms.values()), dtype=torch.float32)
        centers.append(all_xyz.mean(dim=0))
        positions.append(torch.tensor([atoms['C2'], atoms['C4'], atoms['C6']],
                                      dtype=torch.float32))
        kept.append(i)
    if not centers:
        raise ValueError(f'No usable RNA residues (with C2/C4/C6) in {pdb_path}')
    return torch.stack(centers), torch.stack(positions), kept


def predict(model, pdb_path):
    """Predict (N, 2) [SHAPE, DMS] reactivity for a PDB structure."""
    centers, positions, kept = coords_and_frames(pdb_path)
    _, R = frame_from_positions(positions)
    with torch.no_grad():
        pred = model(centers, R)
    return pred, kept


if __name__ == '__main__':
    import sys
    import sgnm
    model = sgnm.load(sys.argv[1]); model.eval()
    pred, kept = predict(model, sys.argv[2])
    shape = pred[:, 0].tolist()
    print('SHAPE:', ','.join(f'{x:.4f}' for x in shape))
