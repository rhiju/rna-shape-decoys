#!/usr/bin/env python3
"""
Convert CASP17 TS-format predictions to standard-column PDB.

Some groups (e.g. TS278) write atom names with non-standard column
justification, which breaks fixed-column parsers (TMscore finds 0 residues in
common -> RMSD defaults to 0; DSSR fails). We parse each ATOM record by
whitespace (all 197 R2307 files have a consistent 12-field layout) and re-emit
strict standard PDB columns with proper atom-name justification.
"""
from pathlib import Path


def format_atom_name(name):
    """Right-justify per PDB convention: single-char-element names start at col 14."""
    if len(name) >= 4:
        return name[:4]
    return ' ' + name.ljust(3)


def convert_ts_to_pdb(ts_file, pdb_file):
    out = []
    with open(ts_file) as f:
        for line in f:
            if line.startswith('ENDMDL') or line.startswith('END'):
                break
            if not line.startswith('ATOM'):
                continue
            parts = line.split()
            # ATOM serial name resName chain resSeq x y z occ bfac element
            _, serial, name, res_name, chain, res_seq, x, y, z, occ, bfac, elem = parts[:12]
            chain = 'A' if chain in ('0', '.', '') else chain[0]
            an = format_atom_name(name)
            out.append(
                f"ATOM  {int(serial):>5} {an}{'':1}{res_name:>3} {chain}{int(res_seq):>4}    "
                f"{float(x):8.3f}{float(y):8.3f}{float(z):8.3f}"
                f"{float(occ):6.2f}{float(bfac):6.2f}          {elem:>2}\n"
            )
    out.append('END\n')
    Path(pdb_file).write_text(''.join(out))


if __name__ == '__main__':
    casp17_dir = Path('data/casp17/R2307')
    n = 0
    for ts_file in sorted(casp17_dir.glob('R2307TS*')):
        if ts_file.suffix:  # skip already-.pdb
            continue
        convert_ts_to_pdb(ts_file, ts_file.with_suffix('.pdb'))
        n += 1
    print(f"Converted {n} CASP17 TS files to standard PDB")
