#!/usr/bin/env python3
"""Collect Snakemake's per-job benchmark files into one table: <bench>/benchmarks.tsv.

Runs when the pipeline ends. The table has one row per job (a rerun job replaces its old row) plus
the run and commit that wrote it. The per-job files stay in place: Snakemake counts a benchmark
file as a job output, so a missing one makes that job look unfinished.

Each new or rerun job is credited to the run whose time window (run_info/runs.tsv) contains the
benchmark file's modification time. Jobs already in the table keep their run, so copying the
results folder without keeping file times does not change the record.

File names:  <rule>.tsv, <rule>.<space>.tsv, <rule>/<sample>.tsv, <rule>/<space>.<sample>.tsv
Uses only the standard library because it runs outside the conda environments.
"""
import argparse
import csv
import datetime
import os
import sys

TABLE = "benchmarks.tsv"
ID_COLUMNS = ["rule", "space", "sample", "run_id", "commit"]
MEMORY_COLUMNS = ("max_rss", "max_vms", "max_uss", "max_pss")
UNKNOWN = "unknown"
TIME_FORMAT = "%Y-%m-%dT%H:%M:%S%z"


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


def path_from_job(rule, space, sample):
    """The reverse of job_from_path."""
    if sample == "NA":
        return f"{rule}.tsv" if space == "NA" else f"{rule}.{space}.tsv"
    return os.path.join(rule, f"{sample}.tsv" if space == "NA" else f"{space}.{sample}.tsv")


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


def parse_time(text):
    try:
        return datetime.datetime.strptime(text, TIME_FORMAT).timestamp()
    except (TypeError, ValueError):
        return None


def run_windows(path):
    """Recorded runs as (start, end, run_id, commit), oldest first.

    A run ends at its 'finished' time, rounded up to the next second because times are recorded to
    the second. A killed run never recorded an end, so it lasts until the next run started. The
    current run has no end yet.
    """
    try:
        with open(path, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
    except OSError:
        return []
    runs = sorted((parse_time(r.get("started")), parse_time(r.get("finished")), r["run_id"], r["commit"])
                  for r in rows if parse_time(r.get("started")) is not None)
    windows = []
    for i, (start, finished, run_id, commit) in enumerate(runs):
        if finished is not None:
            end = finished + 1
        else:
            end = runs[i + 1][0] if i + 1 < len(runs) else float("inf")
        windows.append((start, end, run_id, commit))
    return windows


def run_of(mtime, windows):
    """(run_id, commit) of the run that was active at mtime; the latest one if two overlap."""
    for start, end, run_id, commit in reversed(windows):
        if start <= mtime < end:
            return run_id, commit
    # Before the first recorded run, or a run started with --no-hooks
    return UNKNOWN, UNKNOWN


ap = argparse.ArgumentParser()
ap.add_argument("--bench", required=True, help="benchmark directory")
ap.add_argument("--runs", required=True, help="runs.tsv written by run_info.py")
args = ap.parse_args()

if not os.path.isdir(args.bench):
    sys.exit(0)

table_path = os.path.join(args.bench, TABLE)
windows = run_windows(args.runs)

table, columns = {}, list(ID_COLUMNS)
if os.path.exists(table_path):
    with open(table_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        columns += [c for c in reader.fieldnames or [] if c not in columns]
        for row in reader:
            table[(row["rule"], row["space"], row["sample"])] = row

seen, updated = set(), 0
for root, _, files in os.walk(args.bench):
    for name in sorted(files):
        path = os.path.join(root, name)
        rel = os.path.relpath(path, args.bench)
        if rel == TABLE or not name.endswith(".tsv"):
            continue
        try:
            key = job_from_path(rel)
            row = read_benchmark(path)
        except (OSError, ValueError, csv.Error) as err:
            print(f"WARNING: skipped {path}: {err}", file=sys.stderr)
            continue
        seen.add(key)
        stored = table.get(key)
        # Same measurements as the stored row: the job has not run again since it was added
        if stored is not None and all(stored.get(c, "NA") == v for c, v in row.items()):
            continue
        run_id, commit = run_of(os.path.getmtime(path), windows)
        row.update(rule=key[0], space=key[1], sample=key[2], run_id=run_id, commit=commit)
        columns += [c for c in row if c not in columns]
        table[key] = row
        updated += 1

# Rows without a file: an earlier version of this script deleted the per-job files after storing
# them, which made those jobs look unfinished. Write them back from the table.
restored = 0
bench_columns = [c for c in columns if c not in ID_COLUMNS]
for key in sorted(set(table) - seen):
    path = os.path.join(args.bench, path_from_job(*key))
    if os.path.exists(path):
        continue
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, bench_columns, delimiter="\t", restval="NA",
                                    extrasaction="ignore", lineterminator="\n")
            writer.writeheader()
            writer.writerow(table[key])
        restored += 1
    except OSError as err:
        print(f"WARNING: could not restore {path}: {err}", file=sys.stderr)

if updated:
    tmp = f"{table_path}.tmp"
    with open(tmp, "w", newline="") as fh:
        writer = csv.DictWriter(fh, columns, delimiter="\t", restval="NA", lineterminator="\n")
        writer.writeheader()
        writer.writerows(table[key] for key in sorted(table))
    os.replace(tmp, table_path)

if updated or restored:
    print(f"Benchmark table {table_path}: {updated} job(s) added or updated, {restored} per-job file(s) restored")
