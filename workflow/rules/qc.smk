### QC: read statistics (seqkit), mapping statistics (samtools flagstat) and one MultiQC
### report over seqkit, flagstat, salmon and oarfish.
### MultiQC names samples from these files: seqkit rows carry the sample name, flagstat files are
### named <sample>.<tx|gn> so transcriptome and genome mapping stay apart, salmon uses the
### counts/salmon/<sample>/ directory name, oarfish comes in as a custom table (oarfish_mqc.tsv).

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
	benchmark:
		f"{BENCH}flagstat/{{space}}.{{sample}}.tsv"
	conda:
		"../envs/mapping.yaml"
	shell:
		"samtools flagstat {input.bam} > {output} 2> {log}"

rule multiqc:
	input:
		seqkit = expand(f"{RESULTS}qc/seqkit/{{sample}}.tsv", sample = SAMPLES),
		flagstat = expand(f"{RESULTS}qc/flagstat/{{sample}}.{{space}}.flagstat", sample = SAMPLES, space = SPACES),
		salmon = expand(f"{SALMON_DIR}{{sample}}/quant.sf", sample = SAMPLES) if "tx" in SPACES else [],
		oarfish = [f"{RESULTS}qc/oarfish/oarfish_mqc.tsv"] if "tx" in SPACES else []
	output:
		html = f"{RESULTS}qc/multiqc/multiqc_report.html"
	params:
		outdir = lambda wildcards, output: os.path.dirname(output.html),
		# salmon's MultiQC files (aux_info/, lib_format_counts.json) sit next to quant.sf
		salmon = lambda wildcards, input: " ".join(os.path.dirname(q) for q in input.salmon),
		run_state = RUN_STATE
	log:
		f"{LOGS}multiqc.log"
	benchmark:
		f"{BENCH}multiqc.tsv"
	conda:
		"../envs/qc.yaml"
	shell:"""
		label=$(awk -F'\\t' '{{print "np_metatx: report built by run " $1 " (commit " $2 ")"}}' {params.run_state} 2> /dev/null || true)
		multiqc --force --comment "$label" --outdir {params.outdir} --filename multiqc_report.html \
		{input.seqkit} {input.flagstat} {params.salmon} {input.oarfish} > {log} 2>&1
	"""
