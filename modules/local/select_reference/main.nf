process SELECT_REFERENCE {
    tag "$meta.id"
    label 'process_low'

    conda "python=3.12"
    container "python:3.12"

    publishDir "${params.outdir}/${meta.id}/reference_selection", mode: params.publish_dir_mode, pattern: "*.reference_selection.csv"

    input:
    tuple val(meta), path(reads)
    path reference_fasta

    output:
    tuple val(meta), path(reads), path("${meta.id}.selected_reference.fasta"), emit: reads_reference
    tuple val(meta), path("${meta.id}.selected_reference.fasta"), emit: selected_reference
    tuple val(meta), path("${meta.id}.reference_selection.csv"), emit: report
    path "versions.yml", emit: versions

    script:
    def readArgs = reads.join(' ')

    """
    select_reference.py \
        --sample-id "${meta.id}" \
        --reads $readArgs \
        --references "$reference_fasta" \
        --selected-fasta ${meta.id}.selected_reference.fasta \
        --output ${meta.id}.reference_selection.csv \
        --min-kmer-identity ${params.reference_min_kmer_identity} \
        --min-matched-read-fraction ${params.reference_min_matched_read_fraction} \
        --min-score-margin ${params.reference_min_score_margin}

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
