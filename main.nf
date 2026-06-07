#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CHIKFLOW } from './workflows/chikflow'

workflow {
    if (params.help) {
        log.info """
        CHIK-FLOW

        Usage:
          nextflow run . --input samplesheet.csv --outdir results -profile docker

        Required arguments:
          --input    CSV samplesheet with columns: sample,fastq_1,fastq_2
          --outdir   Output directory

        Reference arguments:
          --reference_fasta  Optional FASTA reference panel; defaults to bundled CHIKV panel
          --reference_gff    Optional GFF annotation; defaults to bundled NC_004162.2 annotation
          --genotype_references  Optional genotype/source FASTA; defaults to bundled CHIKV panel
          --reference_min_kmer_identity  Minimum selected-reference k-mer identity for high confidence
          --reference_min_matched_read_fraction  Minimum read fraction with a matching k-mer for high confidence
          --reference_min_score_margin  Minimum selected-vs-runner-up k-mer margin for high confidence

        Optional arguments:
          --skip_fastqc   Skip FastQC before and after trimming
          --skip_fastp    Skip fastp trimming
          --skip_multiqc  Skip MultiQC aggregation
          --skip_reference_prep  Skip reference panel validation
          --skip_alignment  Skip BWA-MEM alignment and BAM processing
          --skip_coverage   Skip per-base depth and coverage summary
          --skip_consensus  Skip variant calling and consensus FASTA generation
          --skip_genotyping  Skip nearest-reference genotype/lineage assignment
          --skip_phylogeny   Skip batch distance tree generation
          --skip_report      Skip final batch HTML/PDF report
          --min_consensus_depth  Minimum depth for unmasked consensus bases
        """.stripIndent()
        return
    }

    if (!params.input) {
        error "Missing required parameter: --input"
    }
    if (!params.outdir) {
        error "Missing required parameter: --outdir"
    }
    if (!params.skip_reference_prep && !params.reference_fasta) {
        error "Missing required parameter: --reference_fasta"
    }

    CHIKFLOW()
}
