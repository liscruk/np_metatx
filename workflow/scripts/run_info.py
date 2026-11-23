#!/usr/bin/env python3
"""Record which code and config produced a pipeline run, in <results>/run_info/.

  runs.tsv            one row per run
  config.<hash>.yml   the resolved config (run configs are not in git), stored once per version
  <run_id>.patch      uncommitted changes under workflow/, only when there are any
  .current            the active run; read by the benchmark table and the MultiQC header
"""
import argparse
import csv
import hashlib
import json
import os
import secrets
import subprocess
import sys
import time

import yaml

COLUMNS = ["run_id", "started", "finished", "status", "commit", "branch", "dirty", "config", "snakemake", "mode"]


def git(repo, *args):
    """Git output as bytes, or None if git is missing or the command fails."""
    try:
        res = subprocess.run(["git", "-C", repo, *args], capture_output=True)
    except OSError:
        return None
    return res.stdout if res.returncode == 0 else None


def git_text(repo, *args):
    out = git(repo, *args)
    return None if out is None else out.decode().strip()


def workflow_patch(top, rel):
    """Diff of tracked changes plus untracked files under the workflow directory."""
    patch = git(top, "diff", "--no-ext-diff", "--binary", "HEAD", "--", rel) or b""
    untracked = git(top, "ls-files", "--others", "--exclude-standard", "-z", "--", rel) or b""
    for path in filter(None, untracked.decode().split("\0")):
        # `git diff --no-index` exits 1 when the files differ, so this does not use git()
        res = subprocess.run(["git", "-C", top, "diff", "--no-ext-diff", "--no-index", "--binary",
                              "--", os.devnull, path], capture_output=True)
        patch += res.stdout
    return patch


def code_version(workflow_dir, run_dir, run_id):
    """Return (commit label, branch, dirty). A plain copy of the pipeline gives 'unknown'."""
    workflow_dir = os.path.realpath(workflow_dir)
    # Guards against a copy of the pipeline sitting inside some other git repository
    if git(workflow_dir, "ls-files", "--error-unmatch", "snakefile") is None:
        return "unknown", "unknown", "unknown"
    top = git_text(workflow_dir, "rev-parse", "--show-toplevel")
    commit = git_text(workflow_dir, "rev-parse", "--short", "HEAD")
    branch = git_text(workflow_dir, "rev-parse", "--abbrev-ref", "HEAD")
    if not top or not commit:
        return "unknown", "unknown", "unknown"
    rel = os.path.relpath(workflow_dir, os.path.realpath(top))
    if not git_text(top, "status", "--porcelain", "--", rel):
        return commit, branch or "unknown", "no"
    with open(os.path.join(run_dir, f"{run_id}.patch"), "wb") as fh:
        fh.write(workflow_patch(top, rel))
    return f"{commit}-dirty", branch or "unknown", "yes"


def store_config(run_dir, config):
    canonical = json.dumps(config, sort_keys=True, default=str)
    name = f"config.{hashlib.sha256(canonical.encode()).hexdigest()[:8]}.yml"
    path = os.path.join(run_dir, name)
    if not os.path.exists(path):
        with open(path, "w") as fh:
            yaml.safe_dump(json.loads(canonical), fh, sort_keys=False)
    return name


def read_runs(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_runs(path, rows):
    tmp = f"{path}.tmp"
    with open(tmp, "w", newline="") as fh:
        writer = csv.DictWriter(fh, COLUMNS, delimiter="\t", restval="", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def now():
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def start(args):
    os.makedirs(args.dir, exist_ok=True)
    started = time.time()
    run_id = f"{time.strftime('%Y%m%dT%H%M%S', time.localtime(started))}-{secrets.token_hex(2)}"
    commit, branch, dirty = code_version(args.repo, args.dir, run_id)
    config_name = store_config(args.dir, json.load(sys.stdin))

    runs_path = os.path.join(args.dir, "runs.tsv")
    rows = read_runs(runs_path)
    # A run still marked as running was killed before it could finish
    for row in rows:
        if row["status"] == "running":
            row["status"] = "interrupted"
    rows.append(dict(zip(COLUMNS, [run_id, now(), "", "running", commit, branch, dirty,
                                   config_name, args.snakemake, args.mode])))
    write_runs(runs_path, rows)
    with open(os.path.join(args.dir, ".current"), "w") as fh:
        fh.write(f"{run_id}\t{commit}\t{started}\n")


def finish(args):
    state = os.path.join(args.dir, ".current")
    if not os.path.exists(state):
        return
    with open(state) as fh:
        run_id = fh.readline().split("\t")[0]
    runs_path = os.path.join(args.dir, "runs.tsv")
    rows = read_runs(runs_path)
    for row in rows:
        if row["run_id"] == run_id:
            row.update(finished=now(), status=args.status)
    write_runs(runs_path, rows)
    os.remove(state)


ap = argparse.ArgumentParser()
sub = ap.add_subparsers(dest="command", required=True)
p = sub.add_parser("start", help="record a new run; reads the resolved config as JSON on stdin")
p.add_argument("--dir", required=True)
p.add_argument("--repo", required=True, help="workflow directory")
p.add_argument("--snakemake", required=True)
p.add_argument("--mode", required=True)
p.set_defaults(func=start)
p = sub.add_parser("finish", help="close the current run")
p.add_argument("--dir", required=True)
p.add_argument("--status", required=True, choices=["success", "error"])
p.set_defaults(func=finish)

args = ap.parse_args()
args.func(args)
