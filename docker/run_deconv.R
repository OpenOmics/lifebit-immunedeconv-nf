#!/usr/bin/env Rscript
# run_deconv.R
# Run immunedeconv on a bulk RNA-seq expression matrix and emit a tidy TSV.
#
# Input:  TPM/CPM matrix TSV with a gene-symbol column (HGNC) and one column per sample.
# Output: long-format TSV with columns: sample, cell_type, fraction, method
#         Also emits a CD4-vs-CD8 summary TSV.
#
# Supported methods (no license): quantiseq, mcp_counter, epic, xcell, abis,
#                                 consensus_tme, estimate, timer
# CIBERSORT/CIBERSORTx: pass --cibersort_binary and --cibersort_mat (LM22).

suppressPackageStartupMessages({
  library(optparse)
  library(immunedeconv)
  library(readr)
})

opt_list <- list(
  make_option(c("-i","--tpm"),   type="character", help="Path to TPM matrix (TSV, gene symbols in col 1)"),
  make_option(c("-m","--method"),type="character", default="quantiseq",
              help="Deconvolution method [default: %default]"),
  make_option(c("-o","--out"),   type="character", default="deconv_long.tsv",
              help="Output long-format TSV [default: %default]"),
  make_option(c("-s","--summary"), type="character", default="cd4_cd8_summary.tsv",
              help="CD4/CD8 summary TSV [default: %default]"),
  make_option("--tumor", action="store_true", default=FALSE,
              help="Set if samples are tumor (affects quantiseq/epic)"),
  make_option("--arrays", action="store_true", default=FALSE,
              help="Set if input is microarray data (affects quantiseq)"),
  make_option("--cibersort_binary", type="character", default=NA,
              help="Path to CIBERSORT.R (only for --method cibersort/cibersort_abs)"),
  make_option("--cibersort_mat", type="character", default=NA,
              help="Path to LM22.txt signature matrix")
)
opt <- parse_args(OptionParser(option_list=opt_list))

if (is.null(opt$tpm)) stop("--tpm is required")

# ---- Load expression matrix ----------------------------------------------
message("[deconv] reading ", opt$tpm)
expr <- as.data.frame(read_tsv(opt$tpm, show_col_types=FALSE))
rownames(expr) <- expr[[1]]
expr[[1]] <- NULL
expr <- as.matrix(expr)
mode(expr) <- "numeric"
message("[deconv] matrix: ", nrow(expr), " genes x ", ncol(expr), " samples")

# ---- CIBERSORT setup (optional) ------------------------------------------
if (opt$method %in% c("cibersort","cibersort_abs")) {
  if (is.na(opt$cibersort_binary) || is.na(opt$cibersort_mat))
    stop("cibersort methods require --cibersort_binary and --cibersort_mat")
  set_cibersort_binary(opt$cibersort_binary)
  set_cibersort_mat(opt$cibersort_mat)
}

# ---- Run deconvolution ----------------------------------------------------
message("[deconv] method = ", opt$method)
res <- deconvolute(
  gene_expression = expr,
  method          = opt$method,
  tumor           = opt$tumor,
  arrays          = opt$arrays
)

# ---- Long format ----------------------------------------------------------
long <- reshape(as.data.frame(res),
                varying = setdiff(colnames(res), "cell_type"),
                v.names = "fraction",
                timevar = "sample",
                times   = setdiff(colnames(res), "cell_type"),
                direction = "long")
long$id <- NULL
long$method <- opt$method
long <- long[, c("sample","cell_type","fraction","method")]
write_tsv(long, opt$out)
message("[deconv] wrote ", opt$out)

# ---- CD4/CD8 summary ------------------------------------------------------
# Cell-type labels vary by method; grep the common patterns.
cd4_pat <- "T cell CD4|CD4\\+ T|T helper|Th |Tregs?|regulatory T"
cd8_pat <- "T cell CD8|CD8\\+ T|Cytotoxic"

is_cd4 <- grepl(cd4_pat, long$cell_type, ignore.case = TRUE)
is_cd8 <- grepl(cd8_pat, long$cell_type, ignore.case = TRUE)

agg <- do.call(rbind, lapply(split(long, long$sample), function(df) {
  data.frame(
    sample     = unique(df$sample),
    cd4_total  = sum(df$fraction[is_cd4[long$sample == unique(df$sample)]], na.rm=TRUE),
    cd8_total  = sum(df$fraction[is_cd8[long$sample == unique(df$sample)]], na.rm=TRUE),
    method     = opt$method,
    stringsAsFactors = FALSE
  )
}))
agg$cd4_cd8_ratio <- ifelse(agg$cd8_total > 0, agg$cd4_total / agg$cd8_total, NA_real_)
write_tsv(agg, opt$summary)
message("[deconv] wrote ", opt$summary)
