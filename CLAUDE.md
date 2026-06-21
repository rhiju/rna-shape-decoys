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
python3 scripts/08_metrics.py              # decoy-retrieval metrics (best-in-top-k, AUPRC, EF)
python3 scripts/09_compare_heatmap.py      # 3-panel naive/SGNM/ERM comparison heatmap
python3 scripts/07_make_report.py          # results.html + results.md + results.pdf
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

## Upstream PRs / issues (submitted 2026-06-20)

- **hmblair/sgnm#2** (PR): dependency `dlu` -> `dlu-torch` (the published dist name;
  bare `dlu` doesn't exist on PyPI so `pip install sgnm` failed).
- **hmblair/ciffy#1** (issue): mmCIF parser is whitespace-column-sensitive —
  single-spaced `_atom_site` rows fail atom typing (repro on their 9FO9 example:
  0 untyped aligned vs 972 single-spaced). Highest-value fix; in the C parser.
- **hmblair/ciffy#2** (PR): remove duplicate `fallback_version` in pyproject
  (source/sdist builds only; PyPI wheels unaffected).

Skipped (decided not worth a PR):
- dlu dist-vs-import doc: README already shows `pip install dlu-torch` + `import dlu`.
- flash-eq as an sgnm dep: already declared as the `equivariant` optional extra.
  (Its real issue — not on PyPI + Colab wheel coverage — is a flash-eq matter.)
- ciffy native PDB input: large; depends on the parser issue being resolved first.
