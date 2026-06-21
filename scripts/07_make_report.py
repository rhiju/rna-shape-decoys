#!/usr/bin/env python3
"""Generate results.html, results.md (GitHub-friendly), and results.pdf."""
import base64
import subprocess
from pathlib import Path
import pandas as pd

ROOT = Path('.').resolve()

# ---- metrics table (from 08_metrics.py) ----
md_metrics = pd.read_csv('results/discrimination_metrics.csv')
SHOW = ['predictor', 'best_top5', 'best_top10', 'auprc', 'auprc_ci95',
        'auprc_p', 'EF_top10pct']
n_good = int(md_metrics['n_good'].iloc[0])
n_models = int(md_metrics['n_models'].iloc[0])
rand_auprc = md_metrics['random_auprc'].iloc[0]

# ---- figures: (title, path, caption) ----
FIGURES = [
    ('SHAPE vs EXPERIMENTAL data', 'results/shape_vs_rmsd_experimental.png',
     'Predicted SHAPE (naive / SGNM / ERM) correlated to experimental 2A3 SHAPE '
     '(OpenKnot gRNAde P20). Black star = reference (true) structure at its real '
     'correlation (~0.6), not 1.0.'),
    ('SHAPE vs REFERENCE-MODEL simulated data', 'results/shape_vs_rmsd_reference.png',
     'Predicted SHAPE vs the reference-model simulated SHAPE (ceiling analysis); '
     'reference trivially sits at 1.0.'),
    ('Structural-accuracy consistency: base-pair F1 vs RMSD', 'results/f1_vs_rmsd.png',
     'The two structural-accuracy measures agree (lower RMSD -> higher F1).'),
    ('Do the two SHAPE predictors agree?', 'results/sgnm_vs_naive.png',
     'SGNM vs naive SHAPE correlation, for experimental and reference targets.'),
    ('Focused comparison: experimental vs naive / SGNM / ERM',
     'results/compare_heatmap.png',
     'Experimental 2A3 plus each predictor for the reference, the best GOOD '
     '(5.3 A) and best BAD (20.9 A) FARFAR2 decoys; dot-bracket SS overlaid.'),
    ('Secondary-structure SHAPE proxy heatmap', 'results/shape_heatmap_naive.png',
     'Unique DSSR secondary structures (dot-bracket on each cell), white=paired, '
     'light-red=unpaired (0.5), sorted by correlation to the reference.'),
    ('SGNM-predicted SHAPE profile heatmap', 'results/shape_heatmap_sgnm.png',
     'Continuous SGNM-predicted reactivity (per-row normalized), sorted by '
     'correlation to the reference-model profile.'),
    ('ERM-predicted SHAPE profile heatmap', 'results/shape_heatmap_erm.png',
     'Continuous ERM-predicted reactivity (per-row normalized), sorted by '
     'correlation to the reference-model profile.'),
]

INTRO = (f"CASP17 target **R2307** (reference: FARFAR2 `Mol9_reference_UtoG_buildloop`). "
         f"Decoys: FARFAR2 + CASP17 submitted models. Can SHAPE reactivity fish "
         f"low-RMSD models from a mostly-bad pool?")

FINDINGS = [
    f"This is a **retrieval** problem: positive = RMSD < 6 A ({n_good} of ~{n_models} "
    f"models). Metrics: best RMSD recovered in the top-k SHAPE-ranked models, AUPRC "
    f"(random ~= {rand_auprc}), and enrichment factor at top 10%.",
    "**vs experiment, SGNM is the best discriminator** (AUPRC 0.14, enrichment 2.8x, "
    "best-in-top-10 = 5.2 A). Naive barely beats random; ERM, despite being the "
    "stronger sequence-validated model, does not beat SGNM here (AUPRC ~ random).",
    "**The ceiling is modest**: even the true structure's predicted SHAPE only "
    "correlates ~0.57-0.61 with experiment, so no model can discriminate sharply.",
    "Only ~11 positives -> wide CIs; this is **one target**. GPU-run GNM matched "
    "the Mac CPU run exactly (max diff 0.0).",
]

METHODS = [
    "**RMSD**: C1' RMSD to reference via TMscore.",
    "**Naive SHAPE**: DSSR paired/unpaired (0 / 0.5), Pearson r vs target.",
    "**SGNM / ERM SHAPE**: hmblair/sgnm GNM (CPU) and equivariant (GPU) models "
    "predict per-residue reactivity from 3D structure; Pearson r vs target.",
    "**Base-pair F1**: F1 of DSSR base pairs vs reference.",
    "Targets: 'vs experiment' = OpenKnot 2A3 data; 'vs reference' = SHAPE predicted "
    "for the reference structure (ceiling).",
]


def metrics_rows_html():
    body = ''
    for _, r in md_metrics.iterrows():
        body += '<tr>' + ''.join(f'<td>{r[c]}</td>' for c in SHOW) + '</tr>'
    hdr = ''.join(f'<th>{c}</th>' for c in SHOW)
    return f'<table><tr>{hdr}</tr>{body}</table>'


def metrics_rows_md():
    hdr = '| ' + ' | '.join(SHOW) + ' |'
    sep = '| ' + ' | '.join('---' for _ in SHOW) + ' |'
    lines = [hdr, sep]
    for _, r in md_metrics.iterrows():
        lines.append('| ' + ' | '.join(str(r[c]) for c in SHOW) + ' |')
    return '\n'.join(lines)


# ---- HTML (self-contained, base64 images) ----
def b64(p):
    return base64.b64encode(Path(p).read_bytes()).decode()

html = [f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<title>SHAPE decoy discrimination — R2307</title><style>
body {{ font-family: Helvetica, Arial, sans-serif; max-width: 1100px; margin: 2em auto; padding: 0 1em; color: #222; }}
h1 {{ font-size: 28px; }} h2 {{ font-size: 21px; border-bottom: 2px solid #eee; padding-bottom: 4px; margin-top: 1.5em; }}
table {{ border-collapse: collapse; margin: 1em 0; }} td, th {{ border: 1px solid #ccc; padding: 6px 12px; }}
th {{ background: #f5f5f5; }} img {{ max-width: 100%; border: 1px solid #ddd; }}
.cap {{ color: #666; font-size: 14px; margin-bottom: 1.5em; }}</style></head><body>
<h1>Can SHAPE data discriminate RNA 3D structures? (CASP17 R2307)</h1>
<p>{INTRO}</p>
<h2>Bottom line — decoy-retrieval metrics</h2>
{metrics_rows_html()}
<ul>{''.join(f'<li>{f}</li>' for f in FINDINGS)}</ul>"""]
for title, path, cap in FIGURES:
    html.append(f'<h2>{title}</h2>\n<img src="data:image/png;base64,{b64(path)}">\n'
                f'<p class="cap">{cap}</p>')
html.append('<h2>Methods</h2><ul>' + ''.join(f'<li>{m}</li>' for m in METHODS) +
            '</ul></body></html>')
Path('results.html').write_text('\n'.join(html))
print("Saved results.html")

# ---- Markdown (relative image paths; renders on GitHub) ----
md = [f"# Can SHAPE data discriminate RNA 3D structures? (CASP17 R2307)\n",
      INTRO + "\n",
      "## Bottom line — decoy-retrieval metrics\n",
      metrics_rows_md() + "\n",
      *[f"- {f}" for f in FINDINGS], ""]
for title, path, cap in FIGURES:
    md += [f"## {title}\n", f"![{title}]({path})\n", f"*{cap}*\n"]
md += ["## Methods\n", *[f"- {m}" for m in METHODS]]
# strip markdown bold markers? keep them — GitHub renders **.
Path('results.md').write_text('\n'.join(md).replace('`', '`'))
print("Saved results.md")

# ---- PDF via headless Chrome (weasyprint/pandoc unavailable) ----
chrome = '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'
if Path(chrome).exists():
    subprocess.run([chrome, '--headless=new', '--disable-gpu', '--no-pdf-header-footer',
                    f'--print-to-pdf={ROOT}/results.pdf', f'file://{ROOT}/results.html'],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if Path('results.pdf').exists():
        print("Saved results.pdf")
    else:
        print("PDF generation failed (Chrome headless)")
else:
    print("Chrome not found; skipped PDF")
