### QC: read statistics (seqkit), mapping statistics (samtools flagstat) and one MultiQC
### report over seqkit, flagstat and salmon.
### MultiQC names samples from these files: seqkit rows carry the sample name, flagstat files are
### named <sample>.<tx|gn> so transcriptome and genome mapping stay apart, salmon uses the
### counts/<sample>/ directory name.

QC_TARGETS = [f"{RESULTS}qc/multiqc/multiqc_report.html"]

### Read statistics per sample: count, length (incl. N50) and quality (Q20/Q30, average Q).
### The file column is replaced by the sample name so MultiQC rows match the other tools.
rule seqkitStats:
	input:
		reads = lambda wildcards: SAMPLES[wildcards.sample]
	output:
		f"{RESULTS}qc/seqkit/{{sample}}.tsv"
	log:
		f"{LOGS}seqkitStats/{{sample}}.log"
	benchmark:
		f"{BENCH}seqkitStats/{{sample}}.tsv"
	conda:
		"../envs/qc.yaml"
	shell:"""
		seqkit stats --all --tabular {input.reads} 2> {log} | \
		awk -v s={wildcards.sample} 'BEGIN{{FS=OFS="\\t"}} NR > 1 {{$1 = s}} 1' > {output}
	"""

### Mapping statistics. With -N 100 (canonical), use the 'primary mapped' line for map rates.
rule flagstat:
	input:
		bam = f"{RESULTS}bam/{{space}}/{{sample}}.bam"
	output:
		f"{RESULTS}qc/flagstat/{{sample}}.{{space}}.flagstat"
	log:
		f"{LOGS}flagstat/{{sample}}.{{space}}.log"
	conda:
		"../envs/mapping.yaml"
	shell:
		"samtools flagstat {input.bam} > {output} 2> {log}"

rule multiqc:
	input:
		seqkit = expand(f"{RESULTS}qc/seqkit/{{sample}}.tsv", sample = SAMPLES),
		flagstat = expand(f"{RESULTS}qc/flagstat/{{sample}}.{{space}}.flagstat", sample = SAMPLES, space = SPACES),
		salmon = expand(f"{RESULTS}counts/{{sample}}/quant.sf", sample = SAMPLES) if "tx" in SPACES else []
	output:
		html = f"{RESULTS}qc/multiqc/multiqc_report.html"
	params:
		outdir = lambda wildcards, output: os.path.dirname(output.html),
		# salmon's MultiQC files (aux_info/, lib_format_counts.json) sit next to quant.sf
		salmon = lambda wildcards, input: " ".join(os.path.dirname(q) for q in input.salmon)
	log:
		f"{LOGS}multiqc.log"
	conda:
		"../envs/qc.yaml"
	shell:"""
		multiqc --force --outdir {params.outdir} --filename multiqc_report.html \
		{input.seqkit} {input.flagstat} {params.salmon} > {log} 2>&1
	"""
