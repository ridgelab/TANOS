#!/usr/bin/env python3
"""
analyze_jackknife_convergence.py

Compares jackknife support values across different replicate counts to determine
how many replicates are needed for stable (converged) support values.

Usage:
    python3 analyze_jackknife_convergence.py <results_dir> [--ref <n_replicates>]

Arguments:
    results_dir   Path to the TANOS_results directory containing subdirectories
                  like one_replicate/, five_replicates/, 10_replicates/, etc.
    --ref N       Use N replicates as the reference (default: highest found)

Each subdirectory is expected to contain a file named either:
    primary_with_support.treefile
    jackknifeSupport/primary_with_jackknife.suptree

The script extracts support values from Newick-format trees (numbers inside
parentheses before branch lengths, e.g. )99.7:0.012) and computes:
  - Mean absolute difference from the reference run
  - Max absolute difference from the reference run
  - RMS difference from the reference run
  - Number of nodes differing by more than 1% from the reference

Output is printed to stdout as a table; optionally writes a CSV.
An optional PNG plots differences from the reference and between consecutive counts.
"""

import re
import sys
import os
import argparse


# ── helpers ──────────────────────────────────────────────────────────────────

TREEFILE_NAMES = [
    "primary_with_support.treefile",
    os.path.join("jackknifeSupport", "primary_with_jackknife.suptree"),
]

# Map common directory name patterns → integer replicate count
NAME_PATTERNS = [
    (r"^one_replicate$",        1),
    (r"^(\d+)_replicates?$",    None),   # e.g. 35_replicates
    (r"^(\w+)_replicates?$",    None),   # e.g. five_replicates
]

WORD_NUMBERS = {
    "one": 1, "two": 2, "three": 3, "four": 4, "five": 5,
    "six": 6, "seven": 7, "eight": 8, "nine": 9, "ten": 10,
    "eleven": 11, "twelve": 12, "thirteen": 13, "fourteen": 14,
    "fifteen": 15, "sixteen": 16, "seventeen": 17, "eighteen": 18,
    "nineteen": 19, "twenty": 20, "twentyfive": 25, "thirty": 30,
    "thirtyfive": 35, "forty": 40, "fortyfive": 45, "fifty": 50,
}


def parse_replicate_count(dirname):
    """Return integer replicate count from a directory name, or None if unrecognised."""
    d = dirname.lower().rstrip("/")
    # numeric prefix: 35_replicates, 5_replicates, etc.
    m = re.match(r"^(\d+)_replicates?$", d)
    if m:
        return int(m.group(1))
    # word prefix: fifteen_replicates, one_replicate, etc.
    m = re.match(r"^([a-z]+)_replicates?$", d)
    if m:
        word = m.group(1)
        if word in WORD_NUMBERS:
            return WORD_NUMBERS[word]
    return None


def find_treefile(dirpath):
    """Return path to the support treefile inside dirpath, or None."""
    for name in TREEFILE_NAMES:
        p = os.path.join(dirpath, name)
        if os.path.isfile(p):
            return p
    return None


def extract_support_values(treefile):
    """Parse a Newick tree and return a list of internal-node support values."""
    with open(treefile) as f:
        content = f.read()
    # Support values appear as )VALUE: in Newick format
    values = re.findall(r"\)(\d+(?:\.\d+)?)", content)
    return [float(v) for v in values]


def mean(lst):
    return sum(lst) / len(lst)


def rms(lst):
    return (sum(x**2 for x in lst) / len(lst)) ** 0.5


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Convergence analysis of jackknife support values across replicate counts."
    )
    parser.add_argument("results_dir", help="Directory containing per-replicate result subdirectories")
    parser.add_argument("--ref", type=int, default=None,
                        help="Replicate count to use as reference (default: highest found)")
    parser.add_argument("--csv", metavar="FILE", default=None,
                        help="Write results table to a CSV file")
    parser.add_argument("--png", metavar="FILE", default=None,
                        help="Write a convergence plot to a PNG file")
    parser.add_argument("--threshold", type=float, default=1.0,
                        help="Support difference threshold for counting 'unstable' nodes (default: 1.0)")
    args = parser.parse_args()

    results_dir = os.path.abspath(args.results_dir)
    if not os.path.isdir(results_dir):
        sys.exit(f"ERROR: {results_dir} is not a directory")

    # ── discover result subdirectories ───────────────────────────────────────
    data = {}   # {n_reps: [support_values]}
    skipped = []

    for entry in sorted(os.listdir(results_dir)):
        fullpath = os.path.join(results_dir, entry)
        if not os.path.isdir(fullpath):
            continue
        n = parse_replicate_count(entry)
        if n is None:
            skipped.append(entry)
            continue
        treefile = find_treefile(fullpath)
        if treefile is None:
            print(f"  [skip] {entry}/ — no treefile found", file=sys.stderr)
            continue
        values = extract_support_values(treefile)
        if not values:
            print(f"  [skip] {entry}/ — no support values parsed", file=sys.stderr)
            continue
        data[n] = values
        print(f"  [ok]   {entry}/ → {len(values)} nodes  ({treefile.split(results_dir)[-1]})",
              file=sys.stderr)

    if skipped:
        print(f"\n  [ignored dirs] {', '.join(skipped)}", file=sys.stderr)

    if not data:
        sys.exit("ERROR: No usable result directories found.")

    reps_sorted = sorted(data.keys())

    # ── choose reference ──────────────────────────────────────────────────────
    ref_n = args.ref if args.ref is not None else max(reps_sorted)
    if ref_n not in data:
        sys.exit(f"ERROR: Reference replicate count {ref_n} not found in data. "
                 f"Available: {reps_sorted}")
    ref = data[ref_n]
    n_nodes = len(ref)

    # Warn if node counts differ
    for r, vals in data.items():
        if len(vals) != n_nodes:
            print(f"  [warn] {r} replicates has {len(vals)} nodes vs {n_nodes} in reference",
                  file=sys.stderr)

    # ── compute statistics ────────────────────────────────────────────────────
    rows = []
    consecutive = {}
    for i, r in enumerate(reps_sorted):
        vals = data[r]
        n = min(len(vals), n_nodes)
        diffs = [abs(vals[i] - ref[i]) for i in range(n)]
        if i == 0:
            step_mean = step_max = None
        else:
            previous = data[reps_sorted[i - 1]]
            step_n = min(len(previous), len(vals))
            step_diffs = [abs(previous[j] - vals[j]) for j in range(step_n)]
            step_mean = round(mean(step_diffs), 4)
            step_max = round(max(step_diffs), 2)
            consecutive[r] = (reps_sorted[i - 1], step_mean, step_max)
        rows.append({
            "reps":          r,
            "n_nodes":       n,
            "mean_support":  round(mean(vals[:n]), 2),
            "pct_ge95":      round(100 * sum(1 for v in vals[:n] if v >= 95) / n, 1),
            "pct_lt50":      round(100 * sum(1 for v in vals[:n] if v <  50) / n, 1),
            "mean_diff":     round(mean(diffs), 4),
            "max_diff":      round(max(diffs),  2),
            "rms_diff":      round(rms(diffs),  4),
            "nodes_gt_thr":  sum(1 for d in diffs if d > args.threshold),
            "step_mean_diff": step_mean,
            "step_max_diff":  step_max,
        })

    # ── print table ───────────────────────────────────────────────────────────
    thr = args.threshold
    header = (
        f"{'Reps':>6}  {'Nodes':>5}  {'MeanSupp':>9}  "
        f"{'>=95%':>6}  {'<50%':>5}  "
        f"{'MeanDiff':>9}  {'MaxDiff':>8}  {'RMSDiff':>8}  "
        f"{f'Nodes>{thr:.0f}%':>10}"
    )
    separator = "-" * len(header)

    print()
    print(f"Reference: {ref_n} replicates  |  threshold: {thr}%  |  {results_dir}")
    print(separator)
    print(header)
    print(separator)
    for row in rows:
        ref_marker = " <-- ref" if row["reps"] == ref_n else ""
        print(
            f"{row['reps']:>6}  {row['n_nodes']:>5}  {row['mean_support']:>9.2f}  "
            f"{row['pct_ge95']:>6.1f}  {row['pct_lt50']:>5.1f}  "
            f"{row['mean_diff']:>9.4f}  {row['max_diff']:>8.2f}  {row['rms_diff']:>8.4f}  "
            f"{row['nodes_gt_thr']:>10}{ref_marker}"
        )
    print(separator)

    # ── convergence summary ───────────────────────────────────────────────────
    print()
    print("Convergence summary (consecutive replicate counts):")
    print(f"  {'Comparison':>12}  {'MeanDiff':>9}  {'MaxDiff':>8}")
    print(f"  {'-'*35}")
    for r2 in reps_sorted[1:]:
        r1, step_mean, step_max = consecutive[r2]
        print(f"  {r1:>5} -> {r2:<5}  {step_mean:>9.4f}  {step_max:>8.2f}")

    # First replicate count where no node exceeds threshold vs reference
    converged_at = next(
        (row["reps"] for row in rows if row["nodes_gt_thr"] == 0 and row["reps"] != ref_n),
        None
    )
    print()
    if converged_at:
        print(f"  ✓ Values converge to within {thr}% of the {ref_n}-replicate reference at: {converged_at} replicates")
    else:
        print(f"  ✗ Values never fully converge to within {thr}% — consider running more replicates")

    # ── optional CSV output ───────────────────────────────────────────────────
    if args.csv:
        import csv
        fieldnames = list(rows[0].keys())
        with open(args.csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        print(f"\n  CSV written to: {args.csv}")

    if args.png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        counts = [row["reps"] for row in rows]
        fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        axes[0].plot(counts, [row["mean_diff"] for row in rows], label="Mean difference")
        axes[0].plot(counts, [row["max_diff"] for row in rows], label="Maximum difference")
        axes[0].axhline(thr, color="gray", linestyle="--", linewidth=1,
                        label=f"Threshold ({thr:g}%)")
        axes[0].set_title(f"Difference from {ref_n}-replicate reference")
        axes[0].set_ylabel("Absolute support change (%)")
        axes[0].legend()
        axes[0].grid(True, alpha=0.25)

        step_counts = [r for r in reps_sorted[1:]]
        axes[1].plot(step_counts, [consecutive[r][1] for r in step_counts],
                     label="Mean change")
        axes[1].plot(step_counts, [consecutive[r][2] for r in step_counts],
                     label="Maximum change")
        axes[1].set_title("Change between consecutive replicate counts")
        axes[1].set_xlabel("Number of replicates")
        axes[1].set_ylabel("Absolute support change (%)")
        axes[1].legend()
        axes[1].grid(True, alpha=0.25)
        fig.tight_layout()
        fig.savefig(args.png, dpi=160)
        plt.close(fig)
        print(f"PNG written to: {args.png}")


if __name__ == "__main__":
    main()
