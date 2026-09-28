#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================
    lifebit-immunedeconv-nf
    Bulk RNA-seq CD4/CD8 T-cell deconvolution via `immunedeconv` (R)
========================================================================
*/

// ---- Parameters (overridable via CLI or -params-file) ------------------
params.input      = null                          // samplesheet CSV: sample,tpm
params.outdir     = 'results'
params.methods    = 'quantiseq,epic,mcp_counter'  // comma list
params.species    = 'human'                       // 'human' or 'mouse'
params.tumor      = false
params.arrays     = false

// Optional CIBERSORT (user must supply — not in container)
params.cibersort_binary = null
params.cibersort_mat    = null

// ---- Print banner ------------------------------------------------------
log.info """
========================================================================
 lifebit-immunedeconv-nf  v0.1.0
------------------------------------------------------------------------
 input       : ${params.input}
 outdir      : ${params.outdir}
 methods     : ${params.methods}
 species     : ${params.species}
 tumor       : ${params.tumor}
 arrays      : ${params.arrays}
========================================================================
""".stripIndent()

if (!params.input) exit 1, "ERROR: --input samplesheet.csv is required"

// ---- Includes ----------------------------------------------------------
include { PREPARE_MATRIX } from './modules/local/prepare_matrix.nf'
include { DECONVOLUTE    } from './modules/local/deconvolute.nf'
include { MERGE_RESULTS  } from './modules/local/merge_results.nf'

// ---- Workflow ----------------------------------------------------------
workflow {

    // Parse samplesheet:  sample,tpm
    ch_samples = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample || !row.tpm)
                error "Samplesheet must have columns: sample,tpm"
            tuple(row.sample, file(row.tpm, checkIfExists: true))
        }

    // Normalize each raw expression file to canonical `symbol \t <sample>` format
    PREPARE_MATRIX(ch_samples)

    // Cross prepared matrices with methods
    ch_methods = Channel.from(params.methods.tokenize(','))
    ch_jobs    = PREPARE_MATRIX.out.matrix.combine(ch_methods)   // (sample, matrix, method)

    // Run deconvolution
    DECONVOLUTE(ch_jobs)

    // Collect all per-(sample,method) outputs and merge
    MERGE_RESULTS(
        DECONVOLUTE.out.long_tsv.collect(),
        DECONVOLUTE.out.summary_tsv.collect()
    )
}

workflow.onComplete {
    log.info "Pipeline done. Results in: ${params.outdir}"
    log.info "Status: ${workflow.success ? 'OK' : 'FAILED'}"
}
