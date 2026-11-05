### Species composition: reads per species and the ratio between them, per mapping space.
###   tx = hybrid transcriptome BAMs (canonical branch), also reports the salmon EM split
###   gn = hybrid genome BAMs (novel branch)
### Each read is assigned by its primary alignment; reads with primary/secondary hits in both
### species are additionally counted as cross-species ambiguous (needs secondaries in the BAM:
### canonical maps with -N 100, the novel branch uses --secondary=no, so ambiguity is 0 there).

SPECIES_TARGETS = expand(f"{RESULTS}species/{{space}}/species_summary.tsv", space = SPACES) + \
	expand(f"{RESULTS}species/{{space}}/species_composition.{{ext}}", space = SPACES, ext = ["pdf", "png"])

### Reference name -> species tables
rule txSpeciesMap:
	input:
		hx = HOST_TX,
		xx = XENO_TX
	output:
		f"{RESULTS}refs/tx2species.{REFNAME}.tsv"
	params:
		hs = SPECIES[0],
		xs = SPECIES[1]
	log:
		f"{LOGS}txSpeciesMap.log"
	conda:
		"../envs/species.yaml"
	shell:"""
		dc() {{ case "$1" in *.gz) gzip -dc "$1" ;; *) cat "$1" ;; esac; }}
		ids() {{ awk -v s="$1" 'substr($0, 1, 1) == ">" {{print substr($1, 2) "\\t" s}}'; }}
		{{ dc {input.hx} | ids {params.hs}; dc {input.xx} | ids {params.xs}; }} > {output} 2> {log}
	"""

rule gnSpeciesMap:
	input:
		HYB_FAI
	output:
		f"{RESULTS}refs/gn2species.{REFNAME}.tsv"
	params:
		hs = SPECIES[0],
		xs = SPECIES[1]
	log:
		f"{LOGS}gnSpeciesMap.log"
	conda:
		"../envs/species.yaml"
	shell:"""
		awk -F'\\t' -v a={params.hs} -v b={params.xs} \
		'index($1, a "_") == 1 {{print $1 "\\t" a; next}} index($1, b "_") == 1 {{print $1 "\\t" b}}' \
		{input} > {output} 2> {log}
	"""

### Per-sample read counts per species.
### tx BAMs are already grouped by read name (minimap2 order); gn BAMs are coordinate-sorted -> collate.
rule speciesCounts:
	input:
		bam = f"{RESULTS}bam/{{space}}/{{sample}}.bam",
		smap = f"{RESULTS}refs/{{space}}2species.{REFNAME}.tsv"
	output:
		f"{RESULTS}species/{{space}}/{{sample}}.counts.tsv"
	params:
		script = os.path.join(SCRIPTS, "species_counts.py"),
		tmp = f"{RESULTS}species/{{space}}/collate_{{sample}}",
		species = " ".join(SPECIES)
	threads:
		THREADS
	log:
		f"{LOGS}speciesCounts/{{space}}.{{sample}}.log"
	benchmark:
		f"{BENCH}speciesCounts/{{space}}.{{sample}}.tsv"
	conda:
		"../envs/species.yaml"
	shell:"""
		{{
		if [ "{wildcards.space}" = "gn" ]; then
			samtools collate -O -u -@ {threads} -T {params.tmp} {input.bam}
		else
			samtools view -u {input.bam}
		fi | samtools view -F 0x800 - | \
		python {params.script} --map {input.smap} --species {params.species} \
			--sample {wildcards.sample} --space {wildcards.space} --out {output}
		}} 2> {log}
	"""

### Merge samples, compute ratios, plot
rule speciesSummary:
	input:
		counts = expand(f"{RESULTS}species/{{space}}/{{sample}}.counts.tsv", sample = SAMPLES, allow_missing = True),
		salmon = lambda w: expand(f"{RESULTS}counts/{{sample}}/quant.sf", sample = SAMPLES) if w.space == "tx" else [],
		txmap = lambda w: f"{RESULTS}refs/tx2species.{REFNAME}.tsv" if w.space == "tx" else []
	output:
		tsv = f"{RESULTS}species/{{space}}/species_summary.tsv",
		pdf = f"{RESULTS}species/{{space}}/species_composition.pdf",
		png = f"{RESULTS}species/{{space}}/species_composition.png"
	params:
		script = os.path.join(SCRIPTS, "species_summary.py"),
		species = " ".join(SPECIES),
		salmon = lambda w, input: f"--salmon {' '.join(input.salmon)} --txmap {input.txmap}" if w.space == "tx" else "",
		title = lambda w: "Species composition - " + ("transcriptome (tx)" if w.space == "tx" else "genome (gn)"),
		prefix = f"{RESULTS}species/{{space}}/species_composition"
	log:
		f"{LOGS}speciesSummary/{{space}}.log"
	conda:
		"../envs/species.yaml"
	shell:"""
		python {params.script} --counts {input.counts} --species {params.species} {params.salmon} \
			--title "{params.title}" --tsv {output.tsv} --plot {params.prefix} > {log} 2>&1
	"""
