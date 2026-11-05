### Canonical branch: competitive mapping against the hybrid transcriptome + salmon quantification

CANONICAL_TARGETS = [
	f"{RESULTS}counts/tx/numreads.tsv",
	f"{RESULTS}counts/tx/tpm.tsv",
]

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

	# Output must stay in minimap2 read order (NOT coordinate-sorted): salmon groups consecutive
	# records with the same read name into one alignment group. -N 100 keeps all competing
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
		quant = f"{RESULTS}counts/{{sample}}/quant.sf"
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

### Merge per-sample quantifications: estimated read counts (for DE) and TPM
rule salmonMerge:
	input:
		quants = expand(f"{RESULTS}counts/{{sample}}/quant.sf", sample = SAMPLES)
	output:
		counts = f"{RESULTS}counts/tx/numreads.tsv",
		tpm = f"{RESULTS}counts/tx/tpm.tsv"
	params:
		dirs = lambda wildcards, input: " ".join(os.path.dirname(q) for q in input.quants)
	log:
		f"{LOGS}salmonMerge.log"
	conda:
		"../envs/salmon.yaml"
	shell:"""
		{{
		salmon quantmerge --quants {params.dirs} --column numreads --output {output.counts}
		salmon quantmerge --quants {params.dirs} --column tpm --output {output.tpm}
		}} > {log} 2>&1
	"""
