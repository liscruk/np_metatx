### Run tracking: records the commit and config of each run (run_info/) and collects the per-job
### benchmark files into benchmarks/benchmarks.tsv when the run ends.
### The commit is deliberately not a rule param: that would make every new commit trigger reruns.

onstart:
	_track("run_info.py", ["start", "--dir", RUN_INFO, "--repo", workflow.basedir,
		"--snakemake", SNAKEMAKE_VERSION, "--mode", config["mode"]],
		stdin = json.dumps(config, default = str))

onsuccess:
	_finish_run("success")

onerror:
	_finish_run("error")
