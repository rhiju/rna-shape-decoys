#!/usr/bin/env python3
"""
Compute F1 score for base pairs between two secondary structures (dot-bracket notation).
Based on: https://github.com/eternagame/OpenKnotScoreMATLAB/blob/main/matlab/score/get_F1_scores.m
"""

def get_base_pairs(dot_bracket):
    """Extract base-pair set from dot-bracket notation.

    Returns: set of (i, j) tuples with i < j, using 0-based indexing.
    """
    stack = []
    pairs = set()

    for i, c in enumerate(dot_bracket):
        if c in '([{':
            stack.append((i, c))
        elif c in ')]}':
            if not stack:
                continue
            j, open_c = stack.pop()
            # Verify matching bracket types
            if (open_c == '(' and c == ')') or \
               (open_c == '[' and c == ']') or \
               (open_c == '{' and c == '}'):
                pairs.add((j, i))

    return pairs

def compute_f1(reference_ss, predicted_ss):
    """
    Compute F1 score for base pairs.

    Args:
        reference_ss: reference dot-bracket string
        predicted_ss: predicted dot-bracket string

    Returns: (precision, recall, f1)
    """
    ref_pairs = get_base_pairs(reference_ss)
    pred_pairs = get_base_pairs(predicted_ss)

    if len(pred_pairs) == 0 and len(ref_pairs) == 0:
        return 1.0, 1.0, 1.0

    if len(pred_pairs) == 0:
        return 0.0, 0.0, 0.0

    tp = len(ref_pairs & pred_pairs)
    fp = len(pred_pairs - ref_pairs)
    fn = len(ref_pairs - pred_pairs)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0

    return precision, recall, f1

if __name__ == '__main__':
    import sys
    if len(sys.argv) == 3:
        ref = sys.argv[1]
        pred = sys.argv[2]
        p, r, f = compute_f1(ref, pred)
        print(f"Precision: {p:.4f}, Recall: {r:.4f}, F1: {f:.4f}")
    else:
        print("Usage: f1_score.py <reference_ss> <predicted_ss>")
