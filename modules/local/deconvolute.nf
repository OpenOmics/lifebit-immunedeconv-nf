process DECONVOLUTE {
    tag "${sample}:${method}"
    label 'process_medium'
    publishDir "${params.outdir}/per_sample/${method}", mode: 'copy'

    input:
    tuple val(sample), path(tpm), val(method)

    output:
    path "${sample}.${method}.long.tsv",    emit: long_tsv
    path "${sample}.${method}.summary.tsv", emit: summary_tsv
    path "versions.yml",                    emit: versions

    script:
    def tumor_flag  = params.tumor  ? '--tumor'  : ''
    def arrays_flag = params.arrays ? '--arrays' : ''
    def cs_opts     = ''
    if (method in ['cibersort','cibersort_abs']) {
        if (!params.cibersort_binary || !params.cibersort_mat)
            error "cibersort requires --cibersort_binary and --cibersort_mat"
        cs_opts = "--cibersort_binary ${params.cibersort_binary} " +
                  "--cibersort_mat ${params.cibersort_mat}"
    }
    """
    Rscript /opt/scripts/run_deconv.R \\
        --tpm ${tpm} \\
        --method ${method} \\
        --species ${params.species} \\
        --out ${sample}.${method}.long.tsv \\
        --summary ${sample}.${method}.summary.tsv \\
        ${tumor_flag} ${arrays_flag} ${cs_opts}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R:            \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
        immunedeconv: \$(Rscript -e 'cat(as.character(packageVersion("immunedeconv")))')
    END_VERSIONS
    """
}
