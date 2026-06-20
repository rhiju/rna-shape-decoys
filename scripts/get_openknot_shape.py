#!/usr/bin/env python3
"""
Extract experimental SHAPE (2A3) reactivity for the gRNAde P20 design from the
OpenKnot benchmark CSV (highest OpenKnot-score design), and save it aligned to
the R2307 / Mol9 structure positions (1-100).

The design sequence matches the Mol9 modeling sequence exactly except position 1
(the U->G "buildloop" mutation). The 100 experimental reactivity values map
directly to structure positions 1-100.

Requires a local checkout of the OpenKnot data:
  git clone https://github.com/eternagame/OpenKnotAIDesignData /tmp/OpenKnotAIDesignData
"""
import sys
import pandas as pd
import numpy as np

CSV = '/tmp/OpenKnotAIDesignData/Data/OpenKnotBench_data.v4.5.1.csv'

df = pd.read_csv(CSV, low_memory=False)
p20 = df[(df['method'] == 'gRNAde') & (df['puzzle'] == 'P20')]
p20 = p20.sort_values('target_openknot_score', ascending=False)
best = p20.iloc[0]
print(f"Selected gRNAde P20 design, OpenKnot score = {best['target_openknot_score']}",
      file=sys.stderr)

design_seq = str(best['design_sequence'])
assert len(design_seq) == 100, f"expected 100-nt design, got {len(design_seq)}"

# reactivity columns, in order; take the 100 non-null values for this design
rcols = sorted([c for c in df.columns if c.startswith('reactivity_') and 'error' not in c],
               key=lambda x: int(x.split('_')[1]))
vals = [best[c] for c in rcols]
nonnull = [float(v) for v in vals if pd.notna(v)]
assert len(nonnull) == 100, f"expected 100 reactivity values, got {len(nonnull)}"

out = pd.DataFrame({'position': range(1, 101), 'shape': nonnull})
out.to_csv('data/experimental_shape.csv', index=False)
print(f"Wrote data/experimental_shape.csv (100 positions)", file=sys.stderr)
print(f"  design seq: {design_seq}", file=sys.stderr)
print(f"  SHAPE range: [{min(nonnull):.3f}, {max(nonnull):.3f}]", file=sys.stderr)
