# Beneficial fungi are the major driver of root fungal microbiome assembly and network robustness

**Saveria Mosca, Edda Francomano, Meriem Miyassa Aci, Nesma Zakaria Mohamed, Leonardo Schena, Antonino Malacrinò**

## Abstract

# Disclaimer

This repository contains the main components used to process the raw data and to analyze it. Raw data is available at NCBI SRA under the BioProject number `PRJNA1080585`.

Our pipeline included:
* nf-core/ampliseq v2.7.1 [https://github.com/nf-core/ampliseq/](https://github.com/nf-core/ampliseq/)
* MAFFT [https://academic.oup.com/nar/article/30/14/3059/2904316](https://academic.oup.com/nar/article/30/14/3059/2904316)
* FastTree [https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0009490)
* R v4.4.1 [https://www.R-project.org/](https://www.R-project.org/)

# Data processing

Run ampliseq

```bash
nextflow run nf-core/ampliseq -r 2.7.1 -profile singularity \
--input_folder $INDIR \
--FW_primer "CTTGGTCATTTAGAGGAAGTAA" \
--RV_primer "GCTGCGTTCTTCATCGATG" \
--outdir $OUTDIR \
--extension "/*_{1,2}.fastq.gz" \
--illumina_pe_its \
--dada_ref_tax_custom $taxDB/sh_general_release_dynamic_04.04.2024.fasta \
--skip_dada_addspecies \
--skip_qiime \
--skip_barplot \
--skip_abundance_tables \
--skip_alpha_rarefaction \
--skip_diversity_indices \
--skip_ancom \
--ignore_empty_input_files \
--ignore_failed_trimming \
--ignore_failed_filtering \
--max_cpus 16 \
--max_memory '64.GB'

mafft --thread $NTHREADS ASV_seqs.fasta > asv_aligned.fasta

FastTree -gtr -nt < asv_aligned.fasta > tree.tre
```

Build a phyloseq object

```r

```

## Data analysis

### Load libraries

```r

```

