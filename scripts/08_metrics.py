#!/usr/bin/env python3
"""
Decoy-discrimination metrics for each SHAPE predictor.

This is a RETRIEVAL problem (fish low-RMSD models from a mostly-bad pool), not a
correlation problem, so we report:

  - best RMSD recovered in the top-k SHAPE-ranked models (k=1,5,10), vs a random
    top-k baseline and the pool's achievable best; with a permutation p-value
    (how often random top-k beats the SHAPE selection).
  - AUPRC (average precision, positive = RMSD < threshold) -- preferred over
    AUROC because positives are rare; with bootstrap 95% CI and permutation p.
  - Enrichment factor at the top 10% (fraction-good-in-top / fraction-good-all).

Models are ranked by SHAPE score descending (higher correlation = predicted
better structure). Writes results/discrimination_metrics.csv.
"""
import numpy as np
import pandas as pd

RMSD_GOOD = 6.0      # "good" = RMSD < 6 A (target-specific; see CLAUDE.md)
KS = [1, 5, 10]
N_BOOT = 2000
SEED = 0

PREDICTORS = {
    'naive (vs expt)': 'shape_naive_vs_expt',
    'SGNM (vs expt)':  'shape_sgnm_vs_expt',
    'ERM (vs expt)':   'shape_erm_vs_expt',
    'naive (vs ref)':  'shape_naive_vs_ref',
    'SGNM (vs ref)':   'shape_sgnm_vs_ref',
    'ERM (vs ref)':    'shape_erm_vs_ref',
}


def average_precision(scores, labels):
    """AUPRC / average precision; higher score = more likely positive."""
    order = np.argsort(-scores)
    lab = labels[order]
    if lab.sum() == 0:
        return np.nan
    cum_tp = np.cumsum(lab)
    precision_at_k = cum_tp / np.arange(1, len(lab) + 1)
    return float((precision_at_k * lab).sum() / lab.sum())


def best_in_topk(scores, rmsd, k):
    order = np.argsort(-scores)
    return float(np.min(rmsd[order[:k]]))


def enrichment_factor(scores, labels, frac=0.10):
    order = np.argsort(-scores)
    n_top = max(1, int(round(frac * len(scores))))
    top_rate = labels[order[:n_top]].mean()
    base_rate = labels.mean()
    return float(top_rate / base_rate) if base_rate > 0 else np.nan


def main():
    rng = np.random.default_rng(SEED)
    df = pd.read_csv('results/scores.csv')
    df['rmsd'] = pd.to_numeric(df['rmsd'], errors='coerce')

    rows = []
    for name, col in PREDICTORS.items():
        if col not in df.columns:
            continue
        sub = df[df[col].notna() & df['rmsd'].notna()].copy()
        scores = pd.to_numeric(sub[col], errors='coerce').values
        rmsd = sub['rmsd'].values
        labels = (rmsd < RMSD_GOOD).astype(int)
        n, n_pos = len(scores), int(labels.sum())
        pool_best = float(rmsd.min())

        rec = {'predictor': name, 'n_models': n, 'n_good': n_pos,
               'pool_best_rmsd': round(pool_best, 2)}

        # best-in-top-k + random baseline + permutation p
        for k in KS:
            achieved = best_in_topk(scores, rmsd, k)
            rand_best = np.array([np.min(rng.choice(rmsd, size=k, replace=False))
                                  for _ in range(N_BOOT)])
            p_rand_better = float((rand_best <= achieved).mean())
            rec[f'best_top{k}'] = round(achieved, 2)
            rec[f'rand_top{k}_mean'] = round(float(rand_best.mean()), 2)
            rec[f'best_top{k}_p'] = round(p_rand_better, 4)

        # AUPRC + bootstrap CI + permutation p
        auprc = average_precision(scores, labels)
        boot = []
        for _ in range(N_BOOT):
            idx = rng.integers(0, n, n)
            if labels[idx].sum() > 0:
                boot.append(average_precision(scores[idx], labels[idx]))
        lo, hi = np.percentile(boot, [2.5, 97.5])
        perm = np.array([average_precision(rng.permutation(scores), labels)
                         for _ in range(N_BOOT)])
        rec['auprc'] = round(auprc, 3)
        rec['auprc_ci95'] = f'[{lo:.3f}, {hi:.3f}]'
        rec['auprc_p'] = round(float((perm >= auprc).mean()), 4)
        rec['random_auprc'] = round(n_pos / n, 3)  # AUPRC of random ranking

        # enrichment factor at top 10%
        rec['EF_top10pct'] = round(enrichment_factor(scores, labels, 0.10), 2)

        rows.append(rec)

    out = pd.DataFrame(rows)
    out.to_csv('results/discrimination_metrics.csv', index=False)
    pd.set_option('display.width', 200, 'display.max_columns', 50)
    print(out.to_string(index=False))
    print(f"\nPositive class: RMSD < {RMSD_GOOD} A. Lower best_topk_p / auprc_p = "
          f"more significant. EF_top10pct: x-fold over random.")
    print("Saved results/discrimination_metrics.csv")


if __name__ == '__main__':
    main()
