#!/usr/bin/env python3
"""Split merged transcript tables into one read-count and one TPM table per species.

Read counts are copied unchanged. TPM is rescaled so each sample sums to 1e6 within each species,
which keeps values comparable between samples with different host/xeno proportions.
"""
import argparse
import csv
import os

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--counts", required=True, help="merged numreads.tsv")
ap.add_argument("--tpm", required=True, help="merged tpm.tsv")
ap.add_argument("--txmap", required=True, help="TSV: transcript -> species tag")
ap.add_argument("--species", nargs=2, required=True, metavar=("HOST", "XENO"))
ap.add_argument("--outdir", required=True)
args = ap.parse_args()


def read_table(path, **kwargs):
    # QUOTE_NONE: transcript names are used exactly as written
    return pd.read_csv(path, sep="\t", index_col=0, quoting=csv.QUOTE_NONE, **kwargs)


def write_table(table, path, **kwargs):
    table.to_csv(path, sep="\t", quoting=csv.QUOTE_NONE, **kwargs)


# Counts stay text so their values are not altered by a float round trip
counts = read_table(args.counts, dtype=str)
tpm = read_table(args.tpm)

if not counts.index.equals(tpm.index) or list(counts.columns) != list(tpm.columns):
    raise SystemExit("numreads and TPM tables differ in transcripts, order or samples")
if not counts.index.is_unique:
    raise SystemExit("merged tables contain duplicate transcript names")

txmap = pd.read_csv(args.txmap, sep="\t", header=None, names=["name", "species"],
                    index_col=0, quoting=csv.QUOTE_NONE)["species"]
if not txmap.index.is_unique:
    raise SystemExit(f"{args.txmap}: duplicate transcript names")

species = counts.index.map(txmap)
unmapped = int(pd.isna(species).sum())
unknown = set(species.dropna()) - set(args.species)
if unmapped or unknown:
    raise SystemExit(f"Cannot assign species: {unmapped} unmapped transcripts, unknown tags {sorted(unknown)}")

os.makedirs(args.outdir, exist_ok=True)
for sp in args.species:
    keep = species == sp
    write_table(counts.loc[keep], os.path.join(args.outdir, f"{sp}.numreads.tsv"))

    sp_tpm = tpm.loc[keep]
    # A sample with no reads in this species stays at 0
    totals = sp_tpm.sum(axis=0).where(lambda s: s > 0)
    rescaled = (sp_tpm.div(totals, axis=1) * 1e6).fillna(0)
    write_table(rescaled, os.path.join(args.outdir, f"{sp}.tpm.tsv"), float_format="%.8g")
    print(f"{sp}: {int(keep.sum())} transcripts")
