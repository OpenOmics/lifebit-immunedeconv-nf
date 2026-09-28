process MERGE_RESULTS {
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path long_tsvs,    stageAs: 'long/*'
    path summary_tsvs, stageAs: 'summary/*'

    output:
    path 'deconv_all_long.tsv'
    path 'cd4_cd8_summary_all.tsv'
    path 'cd4_cd8_wide.tsv'

    script:
    """
    merge_deconv.R \\
        --long_dir long \\
        --summary_dir summary \\
        --out_long deconv_all_long.tsv \\
        --out_summary cd4_cd8_summary_all.tsv \\
        --out_wide cd4_cd8_wide.tsv
    """
}
