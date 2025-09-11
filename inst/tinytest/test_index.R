library(Bioc.gff)
library(BiocIO)
library(Rsamtools)
library(GenomicRanges)

## Setup: create a compressed and indexed GFF file
gff_file <- system.file("extdata", "genes.gff3", package = "Bioc.gff")
gff_data <- import(gff_file)
sorted_gff <- sort(gff_data)

temp_dir <- tempfile()
dir.create(temp_dir)

sorted_gff_path <- file.path(temp_dir, "sorted.gff3")
export(sorted_gff, sorted_gff_path)

bgzipped_path <- bgzip(sorted_gff_path, overwrite = TRUE)
indexed_path <- indexTabix(bgzipped_path, "gff")

tbi_file <- TabixFile(bgzipped_path)
gff_file_obj <- GFFFile(bgzipped_path)
query_range <- GRanges("chr12", IRanges(11200000, 11200200))
mgr <- BiocIO:::manager()

## Test BiocFile method, no index present
unindexed_gff_obj <- GFFFile(gff_file)
res_no_index <- Bioc.gff:::queryForResource(mgr, unindexed_gff_obj)
expect_false(attr(res_no_index, "usedWhich"))
expect_equal(res_no_index[[1L]], gff_file)

res_no_index_which <- Bioc.gff:::queryForResource(
    mgr, unindexed_gff_obj, which = query_range
)
expect_false(attr(res_no_index_which, "usedWhich"))
expect_equal(res_no_index_which[[1L]], gff_file)

## Cleanup
BiocIO:::release(mgr, res_no_index)
unlink(temp_dir, recursive = TRUE)
