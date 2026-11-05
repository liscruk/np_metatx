### Shared settings, paths and input checks. Included first by the snakefile.

import gzip

### Settings
RESULTS = config["results"]
SAMPLES = config["samples"]
SPECIES = config["species"]
THREADS = config["threads"]
MODE = config["mode"].lower()
REFNAME = f"{SPECIES[0]}.{SPECIES[1]}"

# Mapping spaces that run: tx = hybrid transcriptome (canonical), gn = hybrid genome (novel)
SPACES = (["tx"] if MODE in ("canonical", "both") else []) + (["gn"] if MODE in ("novel", "both") else [])

LOGS = f"{RESULTS}logs/"
BENCH = f"{RESULTS}benchmarks/"
SCRIPTS = os.path.join(workflow.basedir, "scripts")

### Source references
HOST_SEQ = config["hostSeq"]
HOST_GTF = config["hostGtf"]
HOST_TX = config["hostTx"]
XENO_SEQ = config["XenoSeq"]
XENO_GTF = config["XenoGtf"]
XENO_TX = config["XenoTx"]

### Hybrid references (built by makeRef) and their minimap2 indexes
HYB_GENOME = f"{RESULTS}refs/Hybrid_genome.{REFNAME}.fa"
HYB_FAI = f"{HYB_GENOME}.fai"
HYB_GTF = f"{RESULTS}refs/Hybrid_annotation.{REFNAME}.gtf"
HYB_TX = f"{RESULTS}refs/Hybrid_transcriptome.{REFNAME}.fa"
HYB_GENOME_MMI = f"{RESULTS}refs/Hybrid_genome.{REFNAME}.splice.mmi"
HYB_TX_MMI = f"{RESULTS}refs/Hybrid_transcriptome.{REFNAME}.map-ont.mmi"

# {sample} may only match sample names from the config, {space} only tx/gn
wildcard_constraints:
	sample = "|".join(re.escape(s) for s in SAMPLES),
	space = "tx|gn"


### Input checks: catch a reference given under the wrong key (e.g. a FASTA as GTF) before
### anything runs. Only the first records are read. Missing files are left to Snakemake.

def _first_records(path, n):
	opener = gzip.open if path.endswith(".gz") else open
	records = []
	with opener(path, "rt") as fh:
		for line in fh:
			if line.strip() and not line.startswith("#"):
				records.append(line.rstrip("\n"))
				if len(records) == n:
					break
	return records


def _bad_input(key, path, expected):
	sys.exit(f"Config '{key}': {path} does not look like {expected}.")


def check_fasta(key):
	path = config[key]
	if not os.path.exists(path):
		return
	try:
		records = _first_records(path, 1)
	except (OSError, UnicodeDecodeError):
		records = []
	if not records or not records[0].startswith(">"):
		_bad_input(key, path, "a FASTA file (first line should start with '>')")


def check_gtf(key):
	path = config[key]
	if not os.path.exists(path):
		return
	try:
		records = _first_records(path, 5)
	except (OSError, UnicodeDecodeError):
		records = []
	expected = "a GTF file (9 tab-separated columns, gene_id \"...\" attributes; GFF3 is not supported)"
	if not records:
		_bad_input(key, path, expected)
	for rec in records:
		cols = rec.split("\t")
		if len(cols) != 9 or not (cols[3].isdigit() and cols[4].isdigit()) or 'gene_id "' not in cols[8]:
			_bad_input(key, path, expected)


for key in ("hostSeq", "hostTx", "XenoSeq", "XenoTx"):
	check_fasta(key)
for key in ("hostGtf", "XenoGtf"):
	check_gtf(key)
