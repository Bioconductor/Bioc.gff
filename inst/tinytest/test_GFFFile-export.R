## Regression test for Bioconductor/Bioc.gff#1: export() called score() and
## strand() without importing them from BiocGenerics, so it failed unless
## BiocGenerics was attached. Other test files attach GenomicRanges (and with
## it BiocGenerics), so check the imports environment directly as well.

bioc_gff_imports <- parent.env(asNamespace("Bioc.gff"))
for (fn in c("score", "strand"))
    expect_true(
        exists(fn, envir = bioc_gff_imports, inherits = FALSE),
        info = paste0("'", fn, "' is imported from BiocGenerics")
    )

gr <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(c(1, 20, 40), width = 10),
    strand = c("+", "-", "*"),
    type = "exon", gene_id = "g1", transcript_id = "t1"
)
out <- tempfile(fileext = ".gtf")
expect_silent(BiocIO::export(gr, Bioc.gff::GTFFile(out)))

lines <- readLines(out)
lines <- lines[!startsWith(lines, "#")]
expect_identical(length(lines), 3L)
fields <- strsplit(lines, "\t", fixed = TRUE)
expect_identical(vapply(fields, `[[`, character(1L), 6L), rep(".", 3L))
expect_identical(vapply(fields, `[[`, character(1L), 7L), c("+", "-", "."))

unlink(out)
