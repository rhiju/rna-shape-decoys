#!/usr/bin/env python3
"""
Convert RNA PDB files to mmCIF format that ciffy can parse.

ciffy's C parser requires both an _atom_site loop (21 columns) and a
_pdbx_poly_seq_scheme block. BioPython and gemmi don't reliably produce
the poly_seq_scheme block for chain-less FARFAR2 / CASP TS structures, so
we emit a minimal, complete mmCIF directly.

Usage:
    pdb_to_cif.py in.pdb out.cif
    from pdb_to_cif import pdb_to_cif; pdb_to_cif('in.pdb', 'out.cif')
"""
import sys

# mmCIF _atom_site columns (matches RCSB / ciffy example format)
_ATOM_SITE_COLS = [
    'group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id',
    'label_comp_id', 'label_asym_id', 'label_entity_id', 'label_seq_id',
    'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy',
    'B_iso_or_equiv', 'pdbx_formal_charge', 'auth_seq_id', 'auth_comp_id',
    'auth_asym_id', 'auth_atom_id', 'pdbx_PDB_model_num',
]

_POLY_SEQ_COLS = [
    'asym_id', 'entity_id', 'seq_id', 'mon_id', 'ndb_seq_num',
    'pdb_seq_num', 'auth_seq_num', 'pdb_mon_id', 'auth_mon_id',
    'pdb_strand_id', 'pdb_ins_code', 'hetero',
]


def _parse_pdb_atoms(pdb_path):
    """Parse ATOM/HETATM records from a PDB file. Reads first model only."""
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('ENDMDL'):
                break  # only first model
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain = line[21].strip() or 'A'
            res_seq = int(line[22:26])
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            occ = line[54:60].strip() or '1.00'
            bfac = line[60:66].strip() or '0.00'
            elem = line[76:78].strip() or atom_name[0]
            atoms.append({
                'atom_name': atom_name, 'res_name': res_name, 'chain': chain,
                'res_seq': res_seq, 'x': x, 'y': y, 'z': z,
                'occ': occ, 'bfac': bfac, 'elem': elem,
            })
    return atoms


def pdb_to_cif(pdb_path, cif_path, data_name='structure'):
    """Convert a PDB file to a ciffy-parseable mmCIF file."""
    atoms = _parse_pdb_atoms(pdb_path)
    if not atoms:
        raise ValueError(f'No ATOM records found in {pdb_path}')

    # Assign sequential label_seq_id per chain (in order of first appearance)
    label_seq = {}
    chain_counters = {}
    residues = []  # ordered list of (chain, res_seq, res_name, label_seq_id)
    for a in atoms:
        key = (a['chain'], a['res_seq'])
        if key not in label_seq:
            chain_counters.setdefault(a['chain'], 0)
            chain_counters[a['chain']] += 1
            label_seq[key] = chain_counters[a['chain']]
            residues.append((a['chain'], a['res_seq'], a['res_name'], chain_counters[a['chain']]))

    def q(name):
        """Quote atom names containing a prime."""
        return f'"{name}"' if "'" in name else name

    # Single-chain assumption: use the first chain as the polymer entity.
    # (All FARFAR2 / CASP R2307 structures are single-chain RNA.)
    chain_id = residues[0][0]
    one_letter = ''.join(r[2] if r[2] in 'ACGU' else 'N' for r in residues)

    with open(cif_path, 'w') as out:
        out.write(f'data_{data_name}\n#\n')

        # _entry
        out.write(f'_entry.id   {data_name}\n#\n')

        # _entity
        out.write('_entity.id                       1\n')
        out.write('_entity.type                     polymer\n')
        out.write('_entity.src_method               syn\n')
        out.write('_entity.pdbx_description         RNA\n')
        out.write(f'_entity.pdbx_number_of_molecules 1\n#\n')

        # _entity_poly
        out.write('_entity_poly.entity_id                    1\n')
        out.write('_entity_poly.type                         polyribonucleotide\n')
        out.write('_entity_poly.nstd_linkage                 no\n')
        out.write('_entity_poly.nstd_monomer                 no\n')
        out.write(f'_entity_poly.pdbx_seq_one_letter_code     {one_letter}\n')
        out.write(f'_entity_poly.pdbx_seq_one_letter_code_can {one_letter}\n')
        out.write(f'_entity_poly.pdbx_strand_id               {chain_id}\n#\n')

        # _entity_poly_seq loop
        out.write('loop_\n')
        for c in ['entity_id', 'num', 'mon_id', 'hetero']:
            out.write(f'_entity_poly_seq.{c}\n')
        for chain, res_seq, res_name, lseq in residues:
            out.write(f'1 {lseq} {res_name} n\n')
        out.write('#\n')

        # _struct_asym
        out.write(f'_struct_asym.id                          {chain_id}\n')
        out.write('_struct_asym.pdbx_blank_PDB_chainid_flag N\n')
        out.write('_struct_asym.pdbx_modified               N\n')
        out.write('_struct_asym.entity_id                   1\n')
        out.write('_struct_asym.details                     ?\n#\n')

        # _atom_site loop
        out.write('loop_\n')
        for c in _ATOM_SITE_COLS:
            out.write(f'_atom_site.{c}\n')
        for i, a in enumerate(atoms, 1):
            lseq = label_seq[(a['chain'], a['res_seq'])]
            an = q(a['atom_name'])
            out.write(
                f"ATOM {i} {a['elem']} {an} . {a['res_name']} {a['chain']} 1 {lseq} ? "
                f"{a['x']:.3f} {a['y']:.3f} {a['z']:.3f} {a['occ']} {a['bfac']} ? "
                f"{a['res_seq']} {a['res_name']} {a['chain']} {an} 1\n"
            )
        out.write('#\n')

        # _pdbx_poly_seq_scheme loop
        out.write('loop_\n')
        for c in _POLY_SEQ_COLS:
            out.write(f'_pdbx_poly_seq_scheme.{c}\n')
        for chain, res_seq, res_name, lseq in residues:
            out.write(
                f'{chain} 1 {lseq} {res_name} {lseq} {res_seq} {res_seq} '
                f'{res_name} {res_name} {chain} . n\n'
            )
        out.write('#\n')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: pdb_to_cif.py in.pdb out.cif')
        sys.exit(1)
    pdb_to_cif(sys.argv[1], sys.argv[2])
    print(f'Wrote {sys.argv[2]}')
