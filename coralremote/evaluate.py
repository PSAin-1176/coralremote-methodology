"""Accuracy metrics and report formatting (port of ``computeMetrics3`` /
``writeMetricBlock``)."""
from __future__ import annotations

import numpy as np

IDS = (1, 2, 3, 4)
NAMES = ("HPY", "PLE", "BLD", "UNK")


def compute_metrics3(y_true, y_pred):
    """4x4 confusion matrix (rows true, cols predicted, UNK last); OA and
    kappa on the 3x3 core; per-class precision / recall / F1 / support."""
    y_true = np.asarray(y_true).ravel()
    y_pred = np.asarray(y_pred).ravel()
    C = np.array([[np.count_nonzero((y_true == r) & (y_pred == c)) for c in IDS] for r in IDS],
                 dtype=np.int64)
    core = C[:3, :3]
    tot = max(core.sum(), 1)
    oa = np.trace(core) / tot
    pe = np.sum((core.sum(1) / tot) * (core.sum(0) / tot))
    kappa = (oa - pe) / max(1 - pe, 1e-12)
    tp = np.diag(C).astype(float)
    total = C.sum(1)
    prec = tp / np.maximum(C.sum(0), 1)
    rec = tp / np.maximum(total, 1)
    f1 = 2 * prec * rec / np.maximum(prec + rec, 1e-12)
    return dict(C=C, OA=oa, kappa=kappa, prec=prec, rec=rec, f1=f1, tp=tp, total=total)


def format_metric_block(m):
    lines = [f"  OA    : {m['OA']:.6f}", f"  Kappa : {m['kappa']:.6f}"]
    for i, n in enumerate(NAMES):
        lines.append(f"  {n}: Prec {m['prec'][i]:.4f}  Rec {m['rec'][i]:.4f}  "
                     f"F1 {m['f1'][i]:.4f}  Support {int(m['total'][i])}")
    lines.append("  Confusion matrix:")
    for r in m["C"]:
        lines.append("    " + "".join(f"{int(x):8d}" for x in r))
    return "\n".join(lines)


def parse_report_confusion(text):
    """Extract the 4x4 confusion matrix from a MATLAB Stage 2 report."""
    lines = text.splitlines()
    i = next(k for k, l in enumerate(lines) if l.strip() == "Confusion matrix:")
    return np.array([[int(x) for x in lines[i + 1 + r].split()] for r in range(4)], dtype=np.int64)


def compare_date_summary(rows_py, rows_ref, tol=0.0):
    """Cell-by-cell comparison of two date summaries (Table 5).

    Returns (n_cells, n_mismatch, max_abs_diff, details)."""
    ref = {r["date"]: r for r in rows_ref}
    details, maxd, n, bad = [], 0.0, 0, 0
    for r in rows_py:
        q = ref[r["date"]]
        for k in ("HPY", "PLE", "BLD", "UNK"):
            d = abs(float(r[k]) - float(q[k]))
            n += 1
            maxd = max(maxd, d)
            if d > tol:
                bad += 1
                details.append((r["date"], k, r[k], q[k], d))
    return n, bad, maxd, details
