### Novel branch: spliced genomic mapping against the hybrid genome + FLAIR isoform discovery
### 1. Map each sample to the hybrid genome (spliced)
### 2. Pool all samples and build one isoform set with `flair transcriptome` (FLAIR >= 3.0)
### 3. Quantify every sample against that common isoform set with `flair quantify`

FLAIR_PREFIX = f"{RESULTS}flair/{REFNAME}"

NOVEL_TARGETS = [
	f"{FLAIR_PREFIX}.isoforms.gtf",
	f"{FLAIR_PREFIX}.quant.counts.tsv",
	f"{FLAIR_PREFIX}.quant.tpm.tsv",
]

### minimap2 index of the hybrid genome, built once instead of in every mapping job.
### Built with the splice preset so its k-mer settings match the spliced mapping below.
### -I 16G keeps the human + xeno genome in one index part.
rule indexGenome:
	input:
		HYB_GENOME
	output:
		HYB_GENOME_MMI
	threads:
		THREADS
	log:
		f"{LOGS}indexGenome.log"
	benchmark:
		f"{BENCH}indexGenome.tsv"
	conda:
		"../envs/mapping.yaml"
	shell:
		"minimap2 -t {threads} -x splice -I 16G -d {output} {input} 2> {log}"

### Spliced alignment against the hybrid genome (parameters as recommended by FLAIR)
rule minimap2splice:
	input:
		reads = lambda wildcards: SAMPLES[wildcards.sample],
		idx = HYB_GENOME_MMI
	output:
		bam = f"{RESULTS}bam/gn/{{sample}}.bam",
		bai = f"{RESULTS}bam/gn/{{sample}}.bam.bai"
	threads:
		THREADS
	log:
		f"{LOGS}minimap2splice/{{sample}}.log"
	benchmark:
		f"{BENCH}minimap2splice/{{sample}}.tsv"
	conda:
		"../envs/mapping.yaml"
	shell:"""
		{{
		minimap2 -t {threads} -ax splice -s 80 -G 200k --secondary=no {input.idx} {input.reads} | \
		samtools sort -@ {threads} -o {output.bam} -
		samtools index {output.bam}
		}} 2> {log}
	"""

### Pool all samples so isoforms are discovered once across the whole experiment
rule mergeGenomic:
	input:
		bams = expand(f"{RESULTS}bam/gn/{{sample}}.bam", sample = SAMPLES)
	output:
		bam = f"{RESULTS}flair/merged.bam",
		bai = f"{RESULTS}flair/merged.bam.bai"
	threads:
		THREADS
	log:
		f"{LOGS}mergeGenomic.log"
	benchmark:
		f"{BENCH}mergeGenomic.tsv"
	conda:
		"../envs/mapping.yaml"
	shell:"""
		{{
		samtools merge -f -@ {threads} -o {output.bam} {input.bams}
		samtools index {output.bam}
		}} 2> {log}
	"""

### Isoform discovery (correct + collapse in one step)
rule flairTranscriptome:
	input:
		bam = f"{RESULTS}flair/merged.bam",
		bai = f"{RESULTS}flair/merged.bam.bai",
		genome = HYB_GENOME,
		gtf = HYB_GTF
	output:
		bed = f"{FLAIR_PREFIX}.isoforms.bed",
		gtf = f"{FLAIR_PREFIX}.isoforms.gtf",
		fa = f"{FLAIR_PREFIX}.isoforms.fa"
	params:
		prefix = FLAIR_PREFIX
	threads:
		THREADS
	log:
		f"{LOGS}flairTranscriptome.log"
	benchmark:
		f"{BENCH}flairTranscriptome.tsv"
	conda:
		"../envs/flair.yaml"
	shell:"""
		flair transcriptome -b {input.bam} \
		-g {input.genome} \
		-f {input.gtf} \
		-t {threads} \
		--stringent \
		-o {params.prefix} > {log} 2>&1
	"""

### Reads manifest for flair quantify: sample, condition, batch, reads
### Condition/batch are placeholders; --sample_id_only keeps the sample names as column headers.
rule flairManifest:
	output:
		manifest = f"{RESULTS}flair/reads_manifest.tsv"
	log:
		f"{LOGS}flairManifest.log"
	run:
		with open(output.manifest, "w") as fh:
			for s, fq in SAMPLES.items():
				fh.write(f"{s}\tcond\tbatch1\t{fq}\n")
		with open(log[0], "w") as fh:
			fh.write(f"Wrote {len(SAMPLES)} samples to {output.manifest}\n")

### Quantify each sample against the common isoform set
rule flairQuantify:
	input:
		manifest = f"{RESULTS}flair/reads_manifest.tsv",
		fa = f"{FLAIR_PREFIX}.isoforms.fa",
		bed = f"{FLAIR_PREFIX}.isoforms.bed"
	output:
		counts = f"{FLAIR_PREFIX}.quant.counts.tsv",
		tpm = f"{FLAIR_PREFIX}.quant.tpm.tsv"
	params:
		prefix = f"{FLAIR_PREFIX}.quant",
		tmp = f"{RESULTS}flair/tmp"
	threads:
		THREADS
	log:
		f"{LOGS}flairQuantify.log"
	benchmark:
		f"{BENCH}flairQuantify.tsv"
	conda:
		"../envs/flair.yaml"
	shell:"""
		flair quantify -r {input.manifest} \
		--isoform_bed {input.bed} \
		-i {input.fa} \
		--sample_id_only \
		--tpm \
		--temp_dir {params.tmp} \
		-t {threads} \
		-o {params.prefix} > {log} 2>&1
	"""
