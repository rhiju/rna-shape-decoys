#!/usr/bin/env python3
"""Generate results.html summarizing the SHAPE decoy-discrimination analysis."""
import base64
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import spearmanr

def img_b64(path):
    return base64.b64encode(Path(path).read_bytes()).decode()

df = pd.read_csv('results/scores.csv')
for c in ['rmsd', 'f1_bp', 'shape_naive', 'shape_sgnm_vs_ref']:
    df[c] = pd.to_numeric(df[c], errors='coerce')

def rho(x, y):
    m = df[x].notna() & df[y].notna()
    return spearmanr(df[m][x], df[m][y])[0] if m.sum() > 2 else float('nan')

metrics = {
    'Naive paired/unpaired SHAPE vs RMSD': rho('shape_naive', 'rmsd'),
    'SGNM-predicted SHAPE vs RMSD': rho('shape_sgnm_vs_ref', 'rmsd'),
    'Base-pair F1 vs RMSD': rho('f1_bp', 'rmsd'),
}

n_farfar = (df['source'] == 'farfar2').sum()
n_casp = (df['source'] == 'casp17').sum()

rows = ''.join(
    f"<tr><td>{k}</td><td>{v:+.3f}</td></tr>" for k, v in metrics.items()
)

html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>SHAPE Decoy Discrimination — R2307</title>
<style>
body {{ font-family: Helvetica, Arial, sans-serif; max-width: 1100px; margin: 2em auto; padding: 0 1em; color: #222; }}
h1 {{ font-size: 28px; }} h2 {{ font-size: 22px; border-bottom: 2px solid #eee; padding-bottom: 4px; margin-top: 1.5em; }}
table {{ border-collapse: collapse; margin: 1em 0; }}
td, th {{ border: 1px solid #ccc; padding: 6px 14px; text-align: left; }}
th {{ background: #f5f5f5; }}
img {{ max-width: 100%; border: 1px solid #ddd; margin: 0.5em 0; }}
.cap {{ color: #666; font-size: 14px; margin-bottom: 1.5em; }}
code {{ background: #f0f0f0; padding: 1px 5px; border-radius: 3px; }}
</style></head><body>
<h1>Can SHAPE data discriminate RNA 3D structures?</h1>
<p>CASP17 target <b>R2307</b> (reference: FARFAR2 <code>Mol9_reference_UtoG_buildloop</code>).
Decoys: <b>{n_farfar}</b> FARFAR2 + <b>{n_casp}</b> CASP17 submitted models.</p>

<h2>Discrimination summary (Spearman ρ vs RMSD)</h2>
<p class="cap">More negative ρ = better discrimination (lower RMSD → higher SHAPE agreement).</p>
<table><tr><th>Metric</th><th>Spearman ρ</th></tr>{rows}</table>

<h2>SHAPE score vs RMSD &amp; base-pair F1</h2>
<img src="data:image/png;base64,{img_b64('results/shape_vs_rmsd.png')}">
<p class="cap">Top-right: SGNM-predicted SHAPE correlation to the reference-model profile vs RMSD —
good models cluster near r≈0.8 and fan down toward negative correlation as RMSD grows.</p>

<h2>Secondary-structure SHAPE proxy heatmap</h2>
<img src="data:image/png;base64,{img_b64('results/shape_heatmap_naive.png')}">
<p class="cap">Unique DSSR secondary structures (dot-bracket on each cell), white=paired, red=unpaired,
sorted by correlation to the reference.</p>

<h2>SGNM-predicted SHAPE profile heatmap</h2>
<img src="data:image/png;base64,{img_b64('results/shape_heatmap_sgnm.png')}">
<p class="cap">Continuous SGNM-predicted reactivity (per-row min-max normalized), sorted by
correlation to the reference-model profile.</p>

<h2>Methods</h2>
<ul>
<li><b>RMSD</b>: C1' RMSD to reference via TMscore.</li>
<li><b>Naive SHAPE</b>: DSSR paired/unpaired (0/0.5), Pearson r to reference pattern.</li>
<li><b>SGNM SHAPE</b>: hmblair/sgnm GNM model predicts per-residue reactivity from 3D
structure; Pearson r of each model's profile to the reference-model profile.
(The stronger equivariant model needs CUDA; run on a GPU node for phase 3.)</li>
<li><b>Base-pair F1</b>: F1 of DSSR base pairs vs reference.</li>
</ul>
<p class="cap">Note: experimental SHAPE comparison pending real per-nucleotide data;
current "vs reference" analysis establishes the ceiling on discrimination power.</p>
</body></html>"""

Path('results.html').write_text(html)
print("Saved results.html")
import subprocess
subprocess.run(['open', 'results.html'])
