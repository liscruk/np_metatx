#!/usr/bin/env python3
"""Merge per-sample oarfish .quant files into a read-count and a TPM table (one column per sample),
and write a MultiQC custom-content table from the oarfish meta_info.json files.

TPM follows the usual definition: reads per transcript length, scaled to sum to 1e6 per sample.
oarfish itself only reports read counts.
"""
import argparse
import json

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--quants", nargs="+", required=True, help="<sample>.quant files")
ap.add_argument("--meta", nargs="+", required=True, help="<sample>.meta_info.json files, same order")
ap.add_argument("--samples", nargs="+", required=True, help="sample names, same order")
ap.add_argument("--counts", required=True)
ap.add_argument("--tpm", required=True)
ap.add_argument("--mqc", required=True)
args = ap.parse_args()

if not len(args.quants) == len(args.meta) == len(args.samples):
    raise SystemExit("--quants, --meta and --samples need the same number of entries")

# ---------- counts and TPM ----------
counts, lengths = {}, None
for sample, path in zip(args.samples, args.quants):
    q = pd.read_csv(path, sep="\t", index_col="tname")
    if lengths is None:
        lengths = q["len"]
    elif not q.index.equals(lengths.index):
        raise SystemExit(f"{path}: transcripts differ from the first sample (different reference?)")
    counts[sample] = q["num_reads"]

counts = pd.DataFrame(counts)
counts.index.name = "Name"
rpk = counts.div(lengths, axis=0)
tpm = rpk.div(rpk.sum(axis=0).where(lambda s: s > 0), axis=1) * 1e6

counts.to_csv(args.counts, sep="\t", float_format="%.6g")
tpm.fillna(0).to_csv(args.tpm, sep="\t", float_format="%.6g")

# ---------- MultiQC table ----------
rows = []
for sample, path in zip(args.samples, args.meta):
    with open(path) as fh:
        meta = json.load(fh)
    d, em = meta["discard_table"], meta.get("em", {})
    rows.append({
        "Sample": sample,
        "reads_assigned": d["valid_best_aln"],
        "aln_dropped_score": d["discard_score"],
        "aln_dropped_aligned_fraction": d["discard_aln_frac"],
        "aln_dropped_aligned_length": d["discard_aln_len"],
        "em_iterations": em.get("iterations", ""),
        "em_converged": em.get("converged", ""),
        "model": meta.get("prob_model", ""),
    })

header = """# id: 'oarfish'
# section_name: 'oarfish'
# description: 'Transcript quantification by oarfish (alignment mode, coverage model). Reads assigned = reads with at least one alignment passing the filters; dropped counts are alignments.'
# plot_type: 'table'
"""
with open(args.mqc, "w") as fh:
    fh.write(header)
    pd.DataFrame(rows).to_csv(fh, sep="\t", index=False)

print(f"Merged {len(args.samples)} samples, {len(counts)} transcripts")
