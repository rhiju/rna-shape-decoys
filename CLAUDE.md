# SHAPE decoy discrimination — project notes

Repo: github.com/rhiju/rna-shape-decoys. Goal: test whether SHAPE reactivity
discriminates near-native RNA 3D structures (CASP17 target R2307 / "Mol9").

## Pipeline (run in order, from the project root)

```bash
python3 scripts/get_openknot_shape.py      # experimental 2A3 SHAPE (gRNAde P20)
python3 scripts/convert_casp17_to_pdb.py   # CASP TS -> standard-column PDB
# RMSD vs reference:
python3 ~/src/daslab_tools/structure/rmsd.py -refpdb data/farfar2/Mol9_reference_UtoG_buildloop.pdb \
    -pdb data/farfar2/Mol9.out.*.pdb data/casp17/R2307/R2307TS*.pdb -t > results/rmsd_raw.csv
python3 scripts/02_shape_naive.py          # DSSR paired/unpaired vs expt & ref
# (merge rmsd into results/scores.csv)
python3 scripts/03_shape_sgnm.py           # SGNM GNM SHAPE vs expt & ref
python3 scripts/04_plot.py                 # 4 figures (2 main 2x2 + f1_vs_rmsd + sgnm_vs_naive)
python3 scripts/05_shape_heatmap.py        # secondary-structure heatmap
python3 scripts/06_sgnm_heatmap.py         # SGNM profile heatmap
python3 scripts/07_make_report.py          # results.html
```

Key helpers: `scripts/dssr_util.py` (DSSR SS, strips `&` backbone breaks),
`scripts/sgnm_predict.py` (SGNM SHAPE from PDB, bypasses ciffy's mmCIF loader).

## SGNM install (this Mac) — see ~/.claude memory "sgnm-install-gotchas"

Regular `pip install` (NOT editable) for ~/src/{dlu,ciffy,sgnm}; checkpoints via
`gh release download v2.0.2 -p "*checkpoint*"`. GNM runs on CPU; ERM needs CUDA.

## TODO — next session (deferred)

1. **Run ERM (equivariant) model on a GPU node.** It's the stronger predictor
   (val corr +0.63 SHAPE / +0.75 DMS vs GNM +0.39/+0.35) but needs `flash-eq`
   (CUDA-only). Install `pip install 'sgnm[equivariant]'` on a CUDA box, load
   `equivariant-checkpoint.pth`, add an `shape_erm_vs_expt` / `shape_erm_vs_ref`
   column, and regenerate the scatter figures.
2. **Sanity-check: re-run GNM/SGNM on the GPU node** and confirm predictions
   match the Mac CPU results (rules out platform/numerical issues).
3. **After ERM is in:** draft PRs into the upstream repos to reduce install/usage
   pain points for future users:
   - `hmblair/ciffy`: duplicate `fallback_version` in pyproject (already fixed
     locally); consider native PDB-format input (ciffy currently reads mmCIF only —
     this is the biggest pain point; see pdb_to_cif note below).
   - `hmblair/sgnm`: pyproject requires `dlu` but the dist is `dlu-torch`
     (use `dlu-torch` or document `--no-deps`); editable-install namespace clash.
   - `hmblair/dlu`: document the `dlu-torch` dist-name vs `dlu` import-name.

## pdb_to_cif status (for the ciffy PR)

`scripts/pdb_to_cif.py` writes mmCIF that ciffy will *parse* (all required blocks:
atom_site, entity, entity_poly, entity_poly_seq, struct_asym, pdbx_poly_seq_scheme)
BUT ciffy then mis-types most atoms (code 0) for our PDBs — not yet a clean,
proposable pipeline. The proven path is `sgnm_predict.py` (bypass). Before
proposing PDB support to ciffy, root-cause the atom-typing failure (atom ordering
within residue / hydrogen naming) so a converter or a native reader types atoms
correctly. Do NOT propose pdb_to_cif.py to ciffy as-is.
