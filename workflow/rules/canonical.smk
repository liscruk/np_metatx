### Canonical branch: competitive mapping against the hybrid transcriptome, quantified in
### parallel by salmon and oarfish from the same alignments. Output per quantifier:
###   counts/salmon/   <sample>/quant.sf, numreads.tsv, tpm.tsv, by_species/
###   counts/oarfish/  <sample>/<sample>.quant (+ .meta_info.json, .ambig_info.tsv), numreads.tsv, tpm.tsv, by_species/

CANONICAL_TARGETS = [
	f"{SALMON_DIR}numreads.tsv",
	f"{SALMON_DIR}tpm.tsv",
	f"{OARFISH_DIR}numreads.tsv",
	f"{OARFISH_DIR}tpm.tsv",
] + expand(f"{RESULTS}counts/{{quant}}/by_species/{{species}}.{{table}}.tsv",
	quant = ["salmon", "oarfish"], species = SPECIES, table = ["numreads", "tpm"])

### minimap2 index of the hybrid transcriptome, built once instead of in every mapping job.
### The index stores the preset's k-mer settings, so it is built with the same preset (map-ont)
### used for mapping. -I 16G keeps it in one part: multi-part indexes give wrong secondaries/MAPQ.
rule indexTranscriptome:
	input:
		HYB_TX
	output:
		HYB_TX_MMI
	threads:
		THREADS
	log:
		f"{LOGS}indexTranscriptome.log"
	benchmark:
		f"{BENCH}indexTranscriptome.tsv"
	conda:
		"../envs/mapping.yaml"
	shell:
		"minimap2 -t {threads} -x map-ont -I 16G -d {output} {input} 2> {log}"

### Competitive mapping of each sample against the hybrid transcriptome.
### If you stick to the model it might make sense to keep the reference stable.
rule mapping:
	input:
		reads = lambda wildcards: SAMPLES[wildcards.sample],
		idx = HYB_TX_MMI
	output:
		bam = f"{RESULTS}bam/tx/{{sample}}.bam"
	threads:
		THREADS
	log:
		f"{LOGS}mapping/{{sample}}.log"
	benchmark:
		f"{BENCH}mapping/{{sample}}.tsv"
	conda:
		"../envs/mapping.yaml"

	# Output must stay in minimap2 read order (NOT coordinate-sorted): salmon and oarfish group
	# consecutive records with the same read name into one alignment group. -N 100 keeps all competing
	# isoforms/orthologues as candidates (default -N 5 truncates them).
	shell:"""
		{{
		minimap2 -t {threads} -ax map-ont -N 100 --eqx {input.idx} {input.reads} | \
		samtools view -b -o {output.bam} -
		}} 2> {log}
	"""

### Quantification by salmon for each sample
rule salmonQuant:
	input:
		bam = f"{RESULTS}bam/tx/{{sample}}.bam",
		tx = HYB_TX
	output:
		quant = f"{SALMON_DIR}{{sample}}/quant.sf"
	params:
		outdir = lambda wildcards, output: os.path.dirname(output.quant)
	threads:
		THREADS
	log:
		f"{LOGS}salmonQuant/{{sample}}.log"
	benchmark:
		f"{BENCH}salmonQuant/{{sample}}.tsv"
	conda:
		"../envs/salmon.yaml"
	# --ont: salmon's ONT alignment model (undocumented on readthedocs, listed in
	# `salmon quant --help-alignment`); uses alignment errors/clipping instead of ignoring them.
	shell:
		"salmon quant --ont -p {threads} -l U -t {input.tx} -a {input.bam} -o {params.outdir} > {log} 2>&1"

### Merge per-sample salmon quantifications: estimated read counts (for DE) and TPM
rule salmonMerge:
	input:
		quants = expand(f"{SALMON_DIR}{{sample}}/quant.sf", sample = SAMPLES)
	output:
		counts = f"{SALMON_DIR}numreads.tsv",
		tpm = f"{SALMON_DIR}tpm.tsv"
	params:
		dirs = lambda wildcards, input: " ".join(os.path.dirname(q) for q in input.quants)
	log:
		f"{LOGS}salmonMerge.log"
	benchmark:
		f"{BENCH}salmonMerge.tsv"
	conda:
		"../envs/salmon.yaml"
	shell:"""
		{{
		salmon quantmerge --quants {params.dirs} --column numreads --output {output.counts}
		salmon quantmerge --quants {params.dirs} --column tpm --output {output.tpm}
		}} > {log} 2>&1
	"""

### Quantification by oarfish for each sample (alignment mode, same BAM as salmon).
### --model-coverage models where reads land on a transcript, which helps with 5'-truncated ONT reads.
### Filters are pinned to oarfish 0.13 defaults so results don't shift if the defaults change:
### alignments need >= 95% of the read's best score (-s), >= 50% of the read (-m) and >= 50 nt (-l) aligned.
rule oarfishQuant:
	input:
		bam = f"{RESULTS}bam/tx/{{sample}}.bam"
	output:
		quant = f"{OARFISH_DIR}{{sample}}/{{sample}}.quant",
		meta = f"{OARFISH_DIR}{{sample}}/{{sample}}.meta_info.json",
		ambig = f"{OARFISH_DIR}{{sample}}/{{sample}}.ambig_info.tsv"
	params:
		prefix = f"{OARFISH_DIR}{{sample}}/{{sample}}"
	threads:
		THREADS
	log:
		f"{LOGS}oarfishQuant/{{sample}}.log"
	benchmark:
		f"{BENCH}oarfishQuant/{{sample}}.tsv"
	conda:
		"../envs/oarfish.yaml"
	shell:
		"oarfish -a {input.bam} -o {params.prefix} -j {threads} --model-coverage -s 0.95 -m 0.5 -l 50 > {log} 2>&1"

### Merge per-sample oarfish quantifications into the same layout as salmon (numreads.tsv, tpm.tsv)
### and write a MultiQC table from oarfish's meta_info.json (MultiQC has no oarfish module).
rule oarfishMerge:
	input:
		quants = expand(f"{OARFISH_DIR}{{sample}}/{{sample}}.quant", sample = SAMPLES),
		metas = expand(f"{OARFISH_DIR}{{sample}}/{{sample}}.meta_info.json", sample = SAMPLES)
	output:
		counts = f"{OARFISH_DIR}numreads.tsv",
		tpm = f"{OARFISH_DIR}tpm.tsv",
		mqc = f"{RESULTS}qc/oarfish/oarfish_mqc.tsv"
	params:
		script = os.path.join(SCRIPTS, "oarfish_merge.py"),
		samples = " ".join(SAMPLES)
	log:
		f"{LOGS}oarfishMerge.log"
	benchmark:
		f"{BENCH}oarfishMerge.tsv"
	conda:
		"../envs/species.yaml"
	shell:"""
		python {params.script} --quants {input.quants} --meta {input.metas} --samples {params.samples} \
			--counts {output.counts} --tpm {output.tpm} --mqc {output.mqc} > {log} 2>&1
	"""

### Split the merged tables by species. Counts are unchanged; TPM is rescaled to 1e6 within each
### species so values stay comparable between samples with different host/xeno proportions.
rule splitBySpecies:
	input:
		counts = f"{RESULTS}counts/{{quant}}/numreads.tsv",
		tpm = f"{RESULTS}counts/{{quant}}/tpm.tsv",
		txmap = f"{RESULTS}refs/tx2species.{REFNAME}.tsv"
	output:
		host_counts = f"{RESULTS}counts/{{quant}}/by_species/{SPECIES[0]}.numreads.tsv",
		host_tpm = f"{RESULTS}counts/{{quant}}/by_species/{SPECIES[0]}.tpm.tsv",
		xeno_counts = f"{RESULTS}counts/{{quant}}/by_species/{SPECIES[1]}.numreads.tsv",
		xeno_tpm = f"{RESULTS}counts/{{quant}}/by_species/{SPECIES[1]}.tpm.tsv"
	params:
		script = os.path.join(SCRIPTS, "split_by_species.py"),
		species = " ".join(SPECIES),
		outdir = lambda wildcards, output: os.path.dirname(output.host_tpm)
	wildcard_constraints:
		quant = "salmon|oarfish"
	log:
		f"{LOGS}splitBySpecies/{{quant}}.log"
	benchmark:
		f"{BENCH}splitBySpecies_{{quant}}.tsv"
	conda:
		"../envs/species.yaml"
	shell:"""
		python {params.script} --counts {input.counts} --tpm {input.tpm} --txmap {input.txmap} \
			--species {params.species} --outdir {params.outdir} > {log} 2>&1
	"""
