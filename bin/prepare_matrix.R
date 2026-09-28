#!/usr/bin/env Rscript
# prepare_matrix.R
# Normalize an input expression file into the canonical form required by
# immunedeconv:  column 1 = HGNC/MGI symbol, remaining columns = samples.
#
# Handles the following common input layouts automatically:
#   1. symbol \\t sample1 [\\t sample2 ...]                       (canonical)
#   2. gene_id \\t GeneName \\t sample1 [\\t sample2 ...]           (RSEM/STAR style)
#   3. ensembl_id \\t sample1 [\\t sample2 ...]                   (no symbol column)
#   4. entrez_id \\t sample1 [\\t sample2 ...]                    (numeric IDs)
#
# When both an Ensembl ID and a GeneName are present, GeneName wins (no
# mapping needed). Otherwise IDs are mapped via org.Hs.eg.db / org.Mm.eg.db.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(AnnotationDbi)
})

opt <- parse_args(OptionParser(option_list=list(
  make_option(c("-i","--input"),   type="character", help="Input TPM/CPM file"),
  make_option(c("-s","--sample"),  type="character", help="Sample name (used if single-sample file has generic header)"),
  make_option(c("--species"),      type="character", default="human", help="human|mouse [default: %default]"),
  make_option(c("-o","--out"),     type="character", default="matrix.tsv", help="Output TSV [default: %default]")
)))

if (is.null(opt$input))  stop("--input is required")
if (is.null(opt$sample)) stop("--sample is required")

# ---- Read raw file -------------------------------------------------------
message("[prep] reading ", opt$input)
raw <- as.data.frame(read_tsv(opt$input, show_col_types=FALSE,
                              guess_max=100000, progress=FALSE))
message("[prep] raw shape: ", nrow(raw), " rows x ", ncol(raw), " cols")
message("[prep] raw columns: ", paste(colnames(raw), collapse=", "))

# ---- Detect layout -------------------------------------------------------
# Which columns are numeric (sample expression) vs. character (ID/metadata)?
col_is_numeric <- sapply(raw, function(x) is.numeric(x) || all(is.na(x)))
id_cols   <- colnames(raw)[!col_is_numeric]
data_cols <- colnames(raw)[col_is_numeric]

if (length(data_cols) < 1)
  stop("[prep] no numeric columns detected — is this a valid TPM file?")
message("[prep] id columns:   ", paste(id_cols, collapse=", "))
message("[prep] data columns: ", paste(data_cols, collapse=", "))

# ---- Pick the symbol column ---------------------------------------------
# Prefer an explicit gene-symbol column if we can find one.
symbol_candidates <- c("GeneName","gene_name","symbol","Symbol",
                       "hgnc_symbol","mgi_symbol","Gene","gene")
sym_col <- intersect(symbol_candidates, id_cols)[1]

detect_id_type <- function(ids) {
  s <- head(ids[!is.na(ids)], 200)
  if (all(grepl("^ENS(MUS)?G[0-9]+", s))) return("ensembl")
  if (all(grepl("^[0-9]+$", s)))          return("entrez")
  "symbol"
}

if (!is.na(sym_col)) {
  message("[prep] using symbol column: ", sym_col)
  symbols <- raw[[sym_col]]
} else {
  # No dedicated symbol column — take the first ID column and map
  if (length(id_cols) < 1)
    stop("[prep] no ID column found")
  first_id <- raw[[id_cols[1]]]
  id_type  <- detect_id_type(first_id)
  message("[prep] no symbol column; mapping ", id_cols[1], " (", id_type, ") -> symbol")

  if (id_type == "symbol") {
    symbols <- first_id
  } else {
    orgdb_pkg <- if (opt$species == "human") "org.Hs.eg.db" else "org.Mm.eg.db"
    suppressPackageStartupMessages(library(orgdb_pkg, character.only=TRUE))
    orgdb   <- get(orgdb_pkg)
    key_col <- switch(id_type, ensembl="ENSEMBL", entrez="ENTREZID")

    # Strip Ensembl version suffix
    keys_in <- if (id_type == "ensembl") sub("\\\\..*\$", "", first_id) else first_id
    map <- suppressMessages(AnnotationDbi::select(
      orgdb, keys=unique(keys_in), columns="SYMBOL", keytype=key_col))
    map <- map[!is.na(map\$SYMBOL) & nzchar(map\$SYMBOL), ]
    map <- map[!duplicated(map[[key_col]]), ]
    symbols <- map\$SYMBOL[match(keys_in, map[[key_col]])]
  }
}

# ---- Extract expression matrix ------------------------------------------
expr <- as.matrix(raw[, data_cols, drop=FALSE])
mode(expr) <- "numeric"

# If the file is single-sample and column header looks generic, rename to sample id
if (ncol(expr) == 1) {
  colnames(expr) <- opt\$sample
  message("[prep] single-sample file; renaming data column to '", opt\$sample, "'")
}

# ---- Clean: drop unmapped, non-finite, collapse duplicates --------------
keep <- !is.na(symbols) & nzchar(symbols)
message("[prep] dropping ", sum(!keep), " rows with no symbol (",
        round(100*sum(!keep)/length(symbols),1), "%)")
symbols <- symbols[keep]
expr    <- expr[keep, , drop=FALSE]

# Replace non-finite with 0
n_nf <- sum(!is.finite(expr))
if (n_nf > 0) {
  message("[prep] replacing ", n_nf, " non-finite values with 0")
  expr[!is.finite(expr)] <- 0
}

# Collapse duplicate symbols by mean (typical for Ensembl -> symbol mapping)
if (any(duplicated(symbols))) {
  n_dup <- sum(duplicated(symbols))
  message("[prep] collapsing ", n_dup, " duplicate symbols by mean")
  expr <- rowsum(expr, group=symbols, reorder=FALSE, na.rm=TRUE) /
          as.vector(table(symbols)[unique(symbols)])
  expr[!is.finite(expr)] <- 0
} else {
  rownames(expr) <- symbols
}

# Drop all-zero rows (they can't inform any signature)
nz <- rowSums(expr) > 0
n_zero <- sum(!nz)
if (n_zero > 0) {
  message("[prep] dropping ", n_zero, " all-zero genes")
  expr <- expr[nz, , drop=FALSE]
}

message("[prep] final matrix: ", nrow(expr), " genes x ", ncol(expr), " samples")

# ---- Write in canonical format -------------------------------------------
out_df <- data.frame(gene_symbol=rownames(expr), expr, check.names=FALSE)
write_tsv(out_df, opt\$out)
message("[prep] wrote ", opt\$out)
