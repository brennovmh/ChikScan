process BATCH_SUMMARY {
    tag "batch_summary"
    label 'process_low'

    conda "python=3.12"
    container "python:3.12"

    publishDir "${params.outdir}/batch_reports", mode: params.publish_dir_mode

    input:
    path summaries

    output:
    path "sample_summary.csv", emit: sample_summary
    path "versions.yml", emit: versions

    script:
    """
    aggregate_sample_summaries.py \
        --summary ${summaries} \
        --output sample_summary.csv

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
