### Canonical branch: competitive mapping against the hybrid transcriptome + salmon quantification

CANONICAL_TARGETS = [RESULTS+"counts/tx/numreads.tsv", RESULTS+"counts/tx/tpm.tsv"] + \
	expand(RESULTS+"bam/tx/{sample}.txt", sample = SAMPLES)

### Canonical workflow perform competetive hybrid mapping on a hybrid reference generated per run.
### If you stick to the model it might make sense to keep the reference stable. 
rule mapping:
	input:
		reads = lambda wildcards: config["samples"][wildcards.sample],
		ref = RESULTS + "refs/Hybrid_transcriptome." + REFNAME + ".fa"
	output:
		bam = RESULTS+"bam/tx/{sample}.bam"
	threads:
		config["threads"]
	params:
		ref = RESULTS + "refs/Hybrid_transcriptome." + REFNAME + ".fa"
	conda:
		"../envs/mapping.yaml"

	# Output must stay in minimap2 read order (NOT coordinate-sorted): salmon groups consecutive
	# records with the same read name into one alignment group. -N 100 keeps all competing
	# isoforms/orthologues as candidates (default -N 5 truncates them).
	shell:"""
		minimap2 -t {threads} -ax map-ont -N 100 --eqx -I 16G {params.ref} {input.reads} | \
		samtools view -b -o {output.bam} -
	"""
### Simple mapping statistics gathering. With -N 100, use the 'primary mapped' line for map rates.
rule flagstat:
	input:
		bam = RESULTS+"bam/tx/{sample}.bam"
	output:
		stat = RESULTS+"bam/tx/{sample}.txt"
	conda:
		"../envs/mapping.yaml"
	
	shell:"""
		samtools flagstat {input.bam} > {output.stat};
	"""
### Generation quantifcation by salmon for seperate mapping
rule salmonQuant:
	input:
		bam = RESULTS+"bam/tx/{sample}.bam",
		trp = RESULTS + "refs/Hybrid_transcriptome." + REFNAME + ".fa"
	output:
		trp_quant = RESULTS+"counts/{sample}/quant.sf",
		trp_dir = directory(RESULTS+"counts/{sample}/")
	threads:
		config["threads"]
	conda:
		"../envs/salmon.yaml"
	# --ont: salmon's ONT alignment model (undocumented on readthedocs, listed in
	# `salmon quant --help-alignment`); uses alignment errors/clipping instead of ignoring them.
	shell:
		"salmon quant --ont -p {threads} -l U -t {input.trp} -a {input.bam} -o {output.trp_dir}"

### Merge per-sample quantifications: estimated read counts (for DE) and TPM
rule salmonMerge:
	input:
		dirs = expand(RESULTS+"counts/{sample}/", sample=SAMPLES)
	output:
		counts = RESULTS+"counts/tx/numreads.tsv",
		tpm = RESULTS+"counts/tx/tpm.tsv"
	conda:
		"../envs/salmon.yaml"
	shell:"""
		salmon quantmerge --quants {input.dirs} --column numreads --output {output.counts};
		salmon quantmerge --quants {input.dirs} --column tpm --output {output.tpm};
	"""


