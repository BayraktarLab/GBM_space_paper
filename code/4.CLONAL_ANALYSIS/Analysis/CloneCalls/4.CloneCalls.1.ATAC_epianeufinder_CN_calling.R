#!/usr/bin/env Rscript

# Run epiAneufinder on 10x Cell Ranger ARC ATAC outputs.
#
# Expected layout:
#   <input_root>/<SECTION>/<SAMPLE>/outs/fragments.tsv.gz
#
# Example:
#   ../../data/raw/AT7-BRA-3-FO-3/cellranger-arc201_count_70d7aab5b03db0b188aa0f7937c7ebc9/outs/fragments.tsv.gz
#
# Output layout:
#   ../../atac_cn_calling/<SECTION>/<SAMPLE>/1mb/epiAneufinder_results/

library(epiAneufinder)

get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- match(flag, args)
  if (is.na(hit)) {
    return(default)
  }
  args[[hit + 1]]
}

input_root <- get_arg("--input-root", "../../data/raw")
out_root <- get_arg("--out-root", "../../atac_cn_calling")
blacklist <- get_arg("--blacklist", "../../data/blacklist.bed")
ncores <- as.integer(get_arg("--ncores", "4"))

fragment_files <- list.files(
  input_root,
  pattern = "fragments.tsv.gz$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(fragment_files) == 0) {
  stop(sprintf("No fragments.tsv.gz files found under %s", input_root))
}

for (fragments in fragment_files) {
  parts <- strsplit(normalizePath(fragments), .Platform$file.sep)[[1]]
  section_pos <- grep("^AT[0-9]+-", parts)

  if (length(section_pos) == 0) {
    next
  }

  section_pos <- section_pos[[1]]
  section <- parts[[section_pos]]
  sample <- parts[[section_pos + 1]]
  patient <- strsplit(section, "-", fixed = TRUE)[[1]][[1]]

  sample_out <- file.path(out_root, section, sample, "1mb")
  dir.create(sample_out, recursive = TRUE, showWarnings = FALSE)

  message(sprintf("Running epiAneufinder: %s / %s", section, sample))

  epiAneufinder(
    input = fragments,
    outdir = sample_out,
    blacklist = blacklist,
    windowSize = 1e6,
    genome = "BSgenome.Hsapiens.UCSC.hg38",
    exclude = c("chrX", "chrY", "chrM"),
    ncores = ncores,
    minFrags = 20000,
    minsizeCNV = 0,
    k = 4,
    plotKaryo = TRUE
  )

}
