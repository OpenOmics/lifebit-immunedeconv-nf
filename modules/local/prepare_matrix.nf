process PREPARE_MATRIX {
    tag "${sample}"
    label 'process_low'
    publishDir "${params.outdir}/prepared_matrices", mode: 'copy'

    container 'ghcr.io/openomics/lifebit-immunedeconv-nf:latest'

    input:
    tuple val(sample), path(raw_tpm)

    output:
    tuple val(sample), path("${sample}.matrix.tsv"), emit: matrix
    path  "${sample}.prepare.log",                   emit: log

    script:
    """
    Rscript ${projectDir}/bin/prepare_matrix.R \\
        --input ${raw_tpm} \\
        --sample ${sample} \\
        --species ${params.species} \\
        --out ${sample}.matrix.tsv \\
        2>&1 | tee ${sample}.prepare.log
    """
}
