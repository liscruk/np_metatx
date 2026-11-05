# np_metatx

Snakemake pipeline for meta-/xeno-transcriptomics on Oxford Nanopore long reads.
Reads from mixed samples (e.g. human host + mouse xenograft, or human + *Aspergillus fumigatus*)
are mapped **competitively** against a hybrid host + xeno reference and quantified per transcript,
so species assignment follows from the transcript of origin.

## Modes

| Mode | What it does |
|------|--------------|
| `Canonical` (default) | Build hybrid transcriptome → `minimap2 -ax map-ont` → `samtools flagstat` → `salmon quant` (alignment mode) → `salmon quantmerge` |
| `Novel` | Spliced genomic mapping (`minimap2 -ax splice`) per sample → pool samples → `flair transcriptome` (one isoform set for the experiment) → `flair quantify` per sample |
| `Both` | Both of the above |

Every mode also writes a **species composition** report per mapping type (`tx` = transcriptome, `gn` = genome):
reads per species (by primary alignment), cross-species ambiguous reads (secondary hits in the other species),
unmapped reads, the ratio between the species and, for `tx`, the salmon (EM) based split.

## Layout

```
config/config.default.yml      template config (copy and adapt; project configs are git-ignored)
workflow/snakefile             entry point: config, mode switch, rule all, hybrid reference (makeRef)
workflow/rules/canonical.smk   canonical branch: transcriptome mapping + salmon quantification
workflow/rules/novel.smk       novel branch: spliced genomic mapping + FLAIR isoform discovery/quantification
workflow/rules/species.smk     reads per species + species ratio, table and plots
workflow/scripts/              helper scripts (species counting / summary plots)
workflow/envs/*.yaml           conda environments (mapping, salmon, flair, species)
resources/sequentialMapping.py experimental sequential host→xeno mapping (not wired in)
```

## Requirements

- Snakemake ≥ 7 with conda/mamba
- Host and xeno references: genome FASTA, GTF (not GFF3) and transcriptome FASTA, plain or gzipped.
  Transcript IDs must be distinct between species (e.g. `ENST…` vs `ENSMUST…`).
- Genome contigs and GTF seqnames are prefixed with the species tag in the hybrid reference
  (`hs_chr1`, `mm_chr1`, …), since GENCODE human and mouse share contig names.
- Novel mode needs FLAIR ≥ 3.0 (`flair transcriptome`); `workflow/envs/flair.yaml` pins 3.0.1.

## Usage

```bash
cp config/config.default.yml config/config_myproject.yml
# edit paths, samples, species and mode

snakemake -s workflow/snakefile \
  --configfile config/config_myproject.yml \
  --use-conda --cores 10
```

## Outputs (in `results`)

```
refs/Hybrid_genome.<hs.xx>.fa(.fai)    concatenated genomes, contigs prefixed with species tag
refs/Hybrid_annotation.<hs.xx>.gtf     concatenated annotations (same prefixes)
refs/Hybrid_transcriptome.<hs.xx>.fa    concatenated transcriptomes

# Canonical
bam/tx/<sample>.bam / .txt             transcriptome alignments + flagstat
counts/<sample>/quant.sf               per-sample salmon quantification
counts/tx/numreads.tsv                 merged estimated read counts (use for DE)
counts/tx/tpm.tsv                      merged TPM

# Novel
bam/gn/<sample>.bam(.bai) / .txt       spliced genome alignments + flagstat
flair/merged.bam                       pooled alignments used for isoform discovery
flair/<hs.xx>.isoforms.{bed,gtf,fa}    FLAIR isoform set (annotated + novel)
flair/<hs.xx>.quant.counts.tsv         isoform counts per sample
flair/<hs.xx>.quant.tpm.tsv            isoform TPM per sample

# Species composition (per mapping type: tx and/or gn)
species/<tx|gn>/<sample>.counts.tsv           per-sample counts: total, unmapped, per species, ambiguous
species/<tx|gn>/species_summary.tsv           all samples + percentages and <xeno>/<host> ratios
                                              (primary, unique-only, and salmon for tx)
species/<tx|gn>/species_composition.pdf/.png  read composition, reads per species, species ratio
```

## License

See `LICENSE`.
