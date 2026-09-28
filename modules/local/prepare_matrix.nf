process PREPARE_MATRIX {
    tag "${sample}"
    label 'process_low'
    publishDir "${params.outdir}/prepared_matrices", mode: 'copy'

    input:
    tuple val(sample), path(raw_tpm)

    output:
    tuple val(sample), path("${sample}.matrix.tsv"), emit: matrix
    path  "${sample}.prepare.log",                   emit: log

    script:
    """
    set -o pipefail
    prepare_matrix.R \\
        --input ${raw_tpm} \\
        --sample ${sample} \\
        --species ${params.species} \\
        --out ${sample}.matrix.tsv \\
        2>&1 | tee ${sample}.prepare.log
    """
}
