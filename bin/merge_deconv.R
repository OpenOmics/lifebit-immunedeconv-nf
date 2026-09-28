#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(readr)
})
opt <- parse_args(OptionParser(option_list=list(
  make_option("--long_dir"),     make_option("--summary_dir"),
  make_option("--out_long"),     make_option("--out_summary"),
  make_option("--out_wide")
)))

read_all <- function(dir) {
  files <- list.files(dir, pattern="\\.tsv$", full.names=TRUE)
  do.call(rbind, lapply(files, read_tsv, show_col_types=FALSE))
}

long <- read_all(opt$long_dir)
summ <- read_all(opt$summary_dir)

write_tsv(long, opt$out_long)
write_tsv(summ, opt$out_summary)

# Wide: one row per sample, columns = <method>_cd4 / <method>_cd8 / <method>_ratio
wide <- reshape(
  summ,
  idvar     = "sample",
  timevar   = "method",
  direction = "wide"
)
write_tsv(wide, opt$out_wide)
