#!/usr/bin/env python3
"""
Fetch SHAPE data from OpenKnot CSV for method=gRNAde, target=P20, max OpenKnot score.
Save per-nt SHAPE and canonical sequence.
"""
import sys
import pandas as pd
import csv

# Download CSV
csv_url = "https://raw.githubusercontent.com/eternagame/OpenKnotAIDesignData/main/Data/OpenKnotBench_data.v4.5.1.csv"
print(f"Downloading {csv_url}...", file=sys.stderr)
df = pd.read_csv(csv_url)

print(f"CSV shape: {df.shape}", file=sys.stderr)
print(f"Columns: {df.columns.tolist()}", file=sys.stderr)

# Filter method=gRNAde, target=P20
df_filtered = df[(df['method'] == 'gRNAde') & (df['target'] == 'P20')]
print(f"Filtered rows: {len(df_filtered)}", file=sys.stderr)

if df_filtered.empty:
    print("ERROR: No rows with method=gRNAde and target=P20", file=sys.stderr)
    sys.exit(1)

# Sort by OpenKnot score (descending) and pick first
if 'openknot_score' in df_filtered.columns:
    score_col = 'openknot_score'
elif 'score' in df_filtered.columns:
    score_col = 'score'
else:
    score_col = df_filtered.columns[-1]
    print(f"Using column '{score_col}' as score", file=sys.stderr)

df_filtered = df_filtered.sort_values(score_col, ascending=False)
best = df_filtered.iloc[0]

print(f"Best row: {best['method']} target={best['target']} score={best[score_col]}", file=sys.stderr)
print(f"Sequence: {best['sequence'][:50]}...", file=sys.stderr)

# Extract SHAPE column
shape_col = None
for col in ['SHAPE', 'shape', '2A3', '2a3', 'reactivity', 'SHAPE_2A3', 'shape_2a3']:
    if col in best.index:
        shape_col = col
        break

if shape_col is None:
    print("Available columns:", best.index.tolist(), file=sys.stderr)
    print("Trying to parse 'SHAPE' or similar...", file=sys.stderr)
    # Check if SHAPE is a string of comma-separated values
    for col in best.index:
        if isinstance(best[col], str) and ',' in str(best[col])[:100]:
            shape_col = col
            print(f"Using column '{col}' as SHAPE data", file=sys.stderr)
            break

if shape_col is None:
    print("ERROR: Could not find SHAPE column", file=sys.stderr)
    sys.exit(1)

shape_data = best[shape_col]
if isinstance(shape_data, str):
    try:
        shape_values = [float(x.strip()) for x in shape_data.strip().split(',')]
    except ValueError:
        print(f"ERROR: Could not parse SHAPE column as floats: {shape_data[:100]}", file=sys.stderr)
        sys.exit(1)
else:
    print(f"SHAPE column is {type(shape_data)}, not string", file=sys.stderr)
    sys.exit(1)

seq = best['sequence'].strip()
print(f"Sequence length: {len(seq)}", file=sys.stderr)
print(f"SHAPE values count: {len(shape_values)}", file=sys.stderr)

# Save experimental_shape.csv (position, shape)
out_shape_csv = 'data/experimental_shape.csv'
with open(out_shape_csv, 'w') as f:
    w = csv.writer(f)
    w.writerow(['position', 'shape'])
    for i, val in enumerate(shape_values, start=1):
        w.writerow([i, val])

print(f"Saved: {out_shape_csv}", file=sys.stderr)

# Save sequence.fasta
out_fasta = 'data/sequence.fasta'
with open(out_fasta, 'w') as f:
    f.write(f">R2307_gRNAde_P20\n{seq}\n")

print(f"Saved: {out_fasta}", file=sys.stderr)
print("Done!", file=sys.stderr)
