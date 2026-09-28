#!/usr/bin/env Rscript
# run_deconv.R
# Run immunedeconv on a bulk RNA-seq expression matrix and emit a tidy TSV.
#
# Input:  TPM/CPM matrix TSV. Column 1 = gene identifier (HGNC/MGI symbol,
#         Ensembl ID with or without version, or Entrez ID). Remaining
#         columns = one per sample.
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
  library(AnnotationDbi)
})

opt_list <- list(
  make_option(c("-i","--tpm"),   type="character", help="Path to TPM matrix (TSV, gene ID in col 1)"),
  make_option(c("-m","--method"),type="character", default="quantiseq",
              help="Deconvolution method [default: %default]"),
  make_option(c("-o","--out"),   type="character", default="deconv_long.tsv",
              help="Output long-format TSV [default: %default]"),
  make_option(c("-s","--summary"), type="character", default="cd4_cd8_summary.tsv",
              help="CD4/CD8 summary TSV [default: %default]"),
  make_option(c("--species"), type="character", default="human",
              help="human|mouse [default: %default]"),
  make_option(c("--id_type"), type="character", default="auto",
              help="auto|symbol|ensembl|entrez [default: %default]"),
  make_option(c("--tumor"), action="store_true", default=FALSE,
              help="Set if samples are tumor (affects quantiseq/epic)"),
  make_option(c("--arrays"), action="store_true", default=FALSE,
              help="Set if input is microarray data (affects quantiseq)"),
  make_option(c("--cibersort_binary"), type="character", default=NA,
              help="Path to CIBERSORT.R (only for --method cibersort/cibersort_abs)"),
  make_option(c("--cibersort_mat"), type="character", default=NA,
              help="Path to LM22.txt signature matrix")
)
opt <- parse_args(OptionParser(option_list=opt_list))

if (is.null(opt$tpm)) stop("--tpm is required")
if (!opt$species %in% c("human","mouse"))
  stop("--species must be 'human' or 'mouse'")

# ---- Load expression matrix ----------------------------------------------
message("[deconv] reading ", opt$tpm)
expr <- as.data.frame(read_tsv(opt$tpm, show_col_types=FALSE))
gene_col <- expr[[1]]
expr[[1]] <- NULL
expr <- as.matrix(expr)
mode(expr) <- "numeric"
rownames(expr) <- gene_col
message("[deconv] matrix: ", nrow(expr), " genes x ", ncol(expr), " samples")

# ---- Detect gene ID type -------------------------------------------------
detect_id_type <- function(ids) {
  s <- head(ids[!is.na(ids)], 200)
  if (all(grepl("^ENS(MUS)?G[0-9]+", s))) return("ensembl")
  if (all(grepl("^[0-9]+$", s)))          return("entrez")
  "symbol"
}
id_type <- if (opt$id_type == "auto") detect_id_type(rownames(expr)) else opt$id_type
message("[deconv] detected id_type = ", id_type)

# ---- Map to gene symbols (HGNC / MGI) ------------------------------------
if (id_type != "symbol") {
  orgdb_pkg <- if (opt$species == "human") "org.Hs.eg.db" else "org.Mm.eg.db"
  suppressPackageStartupMessages(library(orgdb_pkg, character.only = TRUE))
  orgdb   <- get(orgdb_pkg)
  key_col <- switch(id_type, ensembl = "ENSEMBL", entrez = "ENTREZID")
  sym_col <- if (opt$species == "human") "SYMBOL" else "SYMBOL"  # MGI symbols live in SYMBOL for org.Mm.eg.db

  # Strip Ensembl version suffix (ENSG00000123.4 -> ENSG00000123)
  keys_in <- if (id_type == "ensembl") sub("\\..*$", "", rownames(expr)) else rownames(expr)

  map <- suppressMessages(AnnotationDbi::select(
    orgdb, keys = unique(keys_in), columns = sym_col, keytype = key_col
  ))
  map <- map[!is.na(map[[sym_col]]) & nzchar(map[[sym_col]]), ]
  # Keep first symbol per key to avoid one-to-many blow-up
  map <- map[!duplicated(map[[key_col]]), ]

  sym <- map[[sym_col]][match(keys_in, map[[key_col]])]
  keep <- !is.na(sym)
  message("[deconv] mapped ", sum(keep), "/", length(sym), " ",
          id_type, " -> symbol (", sum(!keep), " dropped)")
  expr <- expr[keep, , drop = FALSE]
  sym  <- sym[keep]

  # Collapse duplicate symbols by mean (common when multiple Ensembl -> same symbol)
  if (any(duplicated(sym))) {
    message("[deconv] collapsing ", sum(duplicated(sym)),
            " duplicate symbols by mean")
    expr <- rowsum(expr, group = sym, reorder = FALSE) /
            as.vector(table(sym)[unique(sym)])
    # rowsum already sets rownames to the group names
  } else {
    rownames(expr) <- sym
  }
  message("[deconv] final matrix: ", nrow(expr), " genes x ", ncol(expr), " samples")
}

# ---- CIBERSORT setup (optional) ------------------------------------------
if (opt$method %in% c("cibersort","cibersort_abs")) {
  if (is.na(opt$cibersort_binary) || is.na(opt$cibersort_mat))
    stop("cibersort methods require --cibersort_binary and --cibersort_mat")
  set_cibersort_binary(opt$cibersort_binary)
  set_cibersort_mat(opt$cibersort_mat)
}

# ---- Run deconvolution ----------------------------------------------------
message("[deconv] method = ", opt$method)
if (opt$species == "mouse") {
  # immunedeconv exposes a separate entry point for mouse (mMCPcounter,
  # seqImmuCC, DCQ, BASE); it handles MGI symbols directly.
  res <- deconvolute_mouse(
    gene_expression_matrix = expr,
    method                 = opt$method
  )
} else {
  res <- deconvolute(
    gene_expression = expr,
    method          = opt$method,
    tumor           = opt$tumor,
    arrays          = opt$arrays
  )
}

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
