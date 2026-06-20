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

## Running ERM (deferred to a GPU — Colab)

ERM is the stronger predictor (val corr +0.63 SHAPE / +0.75 DMS vs GNM
+0.39/+0.35) but is atom-level and needs `flash-eq` (CUDA-only; no CPU/MPS path,
and it can't be validated on the Mac). Run it on a free Colab T4:

1. Open `notebooks/run_erm_colab.ipynb` in Colab, set GPU runtime, Run all.
   It installs the stack, clones this repo, downloads `equivariant-checkpoint.pth`,
   verifies atom typing, runs `scripts/run_erm.py`, and downloads
   `results/erm_profiles.csv`.
2. Locally: drop `erm_profiles.csv` into `results/`, then
   `python3 scripts/ingest_erm.py` — merges `shape_erm_vs_expt` /
   `shape_erm_vs_ref` into scores.csv and regenerates plots, metrics, report.
   (04_plot / 08_metrics auto-include ERM once the columns exist.)
3. **Sanity check on the GPU**: the notebook can also run `scripts/03_shape_sgnm.py`
   to confirm GNM predictions match the Mac CPU results.

## pdb_to_cif — SOLVED (and a ciffy PR candidate)

`scripts/pdb_to_cif.py` now produces mmCIF that ciffy parses **and types
correctly** (code-0 = 0 on FARFAR2 + CASP17). Verified: SGNM via the converted
CIF matches the `sgnm_predict.py` bypass (corr 0.987). This is what makes ERM
possible (ERM needs a real ciffy Polymer; the bypass only works for SGNM).

Root cause of the earlier mis-typing: **ciffy's C mmCIF parser is whitespace-
column sensitive** — collapsing aligned `_atom_site` columns to single spaces
makes it mis-read atom names (even on the official 9FO9 example: 972/1053 atoms
go untyped). The converter now emits per-column left-justified (aligned) rows.

## TODO — PR phase (after ERM is in)

Draft PRs to reduce install/usage pain for future users:
- `hmblair/ciffy`: (a) **whitespace-tolerant mmCIF parsing** (the alignment bug
  above — strong, well-isolated PR with a repro); (b) duplicate `fallback_version`
  in pyproject (fixed locally); (c) optionally native PDB input.
- `hmblair/sgnm`: pyproject requires `dlu` but the dist is `dlu-torch`
  (use `dlu-torch` or document `--no-deps`); editable-install namespace clash.
- `hmblair/dlu`: document the `dlu-torch` dist-name vs `dlu` import-name.
