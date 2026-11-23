#!/usr/bin/env python3
"""Collect Snakemake's per-job benchmark files into one table: <bench>/benchmarks.tsv.

Snakemake needs one file per job during a run. This script runs when the pipeline ends, adds each
file as a row (a rerun job replaces its old row) and deletes the files it has safely stored.

File names:  <rule>.tsv, <rule>.<space>.tsv, <rule>/<sample>.tsv, <rule>/<space>.<sample>.tsv
Uses only the standard library because it runs outside the conda environments.
"""
import argparse
import csv
import os
import sys

TABLE = "benchmarks.tsv"
ID_COLUMNS = ["rule", "space", "sample", "run_id", "commit"]
MEMORY_COLUMNS = ("max_rss", "max_vms", "max_uss", "max_pss")
UNKNOWN = "unknown"


def job_from_path(rel):
    """Return (rule, space, sample) from a path relative to the benchmark directory."""
    parts = rel.split(os.sep)
    stem = parts[-1][:-len(".tsv")]
    if len(parts) == 1:
        rule, _, space = stem.partition(".")
        return rule, space or "NA", "NA"
    if len(parts) == 2:
        space, dot, sample = stem.partition(".")
        return (parts[0], space, sample) if dot else (parts[0], "NA", stem)
    raise ValueError("unexpected directory depth")


def read_benchmark(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows or None in rows[-1]:
        raise ValueError("no complete benchmark row")
    row = rows[-1]
    # Jobs shorter than Snakemake's memory sampling interval report 0
    for col in MEMORY_COLUMNS:
        if col in row and row[col] in ("0", "0.0", "0.00", "-", ""):
            row[col] = "NA"
    return row


def read_state(path):
    """Current run as (run_id, commit, start time), or None outside a tracked run."""
    try:
        with open(path) as fh:
            run_id, commit, started = fh.readline().rstrip("\n").split("\t")
        return run_id, commit, float(started)
    except (OSError, ValueError):
        return None


ap = argparse.ArgumentParser()
ap.add_argument("--bench", required=True, help="benchmark directory")
ap.add_argument("--state", required=True, help="current-run file written by run_info.py")
args = ap.parse_args()

if not os.path.isdir(args.bench):
    sys.exit(0)

table_path = os.path.join(args.bench, TABLE)
state = read_state(args.state)

table, columns = {}, list(ID_COLUMNS)
if os.path.exists(table_path):
    with open(table_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        columns += [c for c in reader.fieldnames or [] if c not in columns]
        for row in reader:
            table[(row["rule"], row["space"], row["sample"])] = row

stored = []
for root, _, files in os.walk(args.bench):
    for name in sorted(files):
        path = os.path.join(root, name)
        rel = os.path.relpath(path, args.bench)
        if rel == TABLE or not name.endswith(".tsv"):
            continue
        try:
            rule, space, sample = job_from_path(rel)
            row = read_benchmark(path)
        except (OSError, ValueError, csv.Error) as err:
            print(f"WARNING: kept {path}: {err}", file=sys.stderr)
            continue
        # Files older than this run came from an earlier, untracked run
        current = state is not None and os.path.getmtime(path) >= state[2]
        row.update(rule=rule, space=space, sample=sample,
                   run_id=state[0] if current else UNKNOWN,
                   commit=state[1] if current else UNKNOWN)
        columns += [c for c in row if c not in columns]
        table[(rule, space, sample)] = row
        stored.append(path)

if not table:
    sys.exit(0)

tmp = f"{table_path}.tmp"
with open(tmp, "w", newline="") as fh:
    writer = csv.DictWriter(fh, columns, delimiter="\t", restval="NA", lineterminator="\n")
    writer.writeheader()
    writer.writerows(table[key] for key in sorted(table))
os.replace(tmp, table_path)

# Only remove per-job files once the table is safely in place
for path in stored:
    os.remove(path)
for root, _, _ in os.walk(args.bench, topdown=False):
    if root != args.bench:
        try:
            os.rmdir(root)
        except OSError:
            pass

print(f"Stored {len(stored)} benchmark file(s) in {table_path}")
