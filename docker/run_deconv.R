#!/usr/bin/env Rscript
# run_deconv.R
# Run immunedeconv on a *pre-prepared* expression matrix.
# Input matrix contract (produced by bin/prepare_matrix.R):
#   - Column 1: HGNC symbol (human) or MGI symbol (mouse)
#   - Columns 2..N: sample TPM/CPM values
#   - No NAs, no all-zero rows, no duplicate symbols

suppressPackageStartupMessages({
  library(optparse)
  library(immunedeconv)
  library(readr)
})

opt_list <- list(
  make_option(c("-i","--tpm"),     type="character", help="Path to prepared matrix (TSV, symbols in col 1)"),
  make_option(c("-m","--method"),  type="character", default="quantiseq"),
  make_option(c("-o","--out"),     type="character", default="deconv_long.tsv"),
  make_option(c("-s","--summary"), type="character", default="cd4_cd8_summary.tsv"),
  make_option(c("--species"),      type="character", default="human"),
  make_option(c("--tumor"),        action="store_true", default=FALSE),
  make_option(c("--arrays"),       action="store_true", default=FALSE),
  make_option(c("--cibersort_binary"), type="character", default=NA),
  make_option(c("--cibersort_mat"),    type="character", default=NA)
)
opt <- parse_args(OptionParser(option_list=opt_list))

if (is.null(opt$tpm)) stop("--tpm is required")
if (!opt$species %in% c("human","mouse")) stop("--species must be 'human' or 'mouse'")

message("[deconv] reading ", opt$tpm)
expr <- as.data.frame(read_tsv(opt$tpm, show_col_types=FALSE))
rownames(expr) <- expr[[1]]
expr[[1]] <- NULL
expr <- as.matrix(expr)
mode(expr) <- "numeric"
message("[deconv] matrix: ", nrow(expr), " genes x ", ncol(expr), " samples")

if (any(!is.finite(expr)))
  stop("[deconv] input matrix contains non-finite values - did prepare_matrix run?")
if (any(duplicated(rownames(expr))))
  stop("[deconv] input matrix has duplicate row names - did prepare_matrix run?")

if (opt$method %in% c("cibersort","cibersort_abs")) {
  if (is.na(opt$cibersort_binary) || is.na(opt$cibersort_mat))
    stop("cibersort methods require --cibersort_binary and --cibersort_mat")
  set_cibersort_binary(opt$cibersort_binary)
  set_cibersort_mat(opt$cibersort_mat)
}

message("[deconv] method = ", opt$method, " (species=", opt$species, ")")
res <- if (opt$species == "mouse") {
  deconvolute_mouse(gene_expression_matrix = expr, method = opt$method)
} else {
  deconvolute(gene_expression = expr, method = opt$method,
              tumor = opt$tumor, arrays = opt$arrays)
}

long <- reshape(as.data.frame(res),
                varying   = setdiff(colnames(res), "cell_type"),
                v.names   = "fraction",
                timevar   = "sample",
                times     = setdiff(colnames(res), "cell_type"),
                direction = "long")
long$id <- NULL
long$method <- opt$method
long <- long[, c("sample","cell_type","fraction","method")]
write_tsv(long, opt$out)
message("[deconv] wrote ", opt$out)

cd4_pat <- "T cell CD4|CD4\\+ T|T helper|Th |Tregs?|regulatory T"
cd8_pat <- "T cell CD8|CD8\\+ T|Cytotoxic"

is_cd4 <- grepl(cd4_pat, long$cell_type, ignore.case = TRUE)
is_cd8 <- grepl(cd8_pat, long$cell_type, ignore.case = TRUE)

agg <- do.call(rbind, lapply(split(long, long$sample), function(df) {
  s <- unique(df$sample)
  data.frame(
    sample        = s,
    cd4_total     = sum(df$fraction[is_cd4[long$sample == s]], na.rm=TRUE),
    cd8_total     = sum(df$fraction[is_cd8[long$sample == s]], na.rm=TRUE),
    method        = opt$method,
    stringsAsFactors = FALSE
  )
}))
agg$cd4_cd8_ratio <- ifelse(agg$cd8_total > 0, agg$cd4_total / agg$cd8_total, NA_real_)
write_tsv(agg, opt$summary)
message("[deconv] wrote ", opt$summary)
