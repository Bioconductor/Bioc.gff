## A failed NCBI lookup for the ##species pragma should warn and import the
## file with empty metadata, not fail in `metadata<-`. A Taxonomy Browser URL
## with no 'id' or 'name' makes metadataFromNCBI() error without a request.

species_url <- paste0(Bioc.gff:::.NCBI_TAX_URL, "?mode=Info")

expect_warning(
    res <- Bioc.gff:::.parseSpeciesAsMetadata(
        paste("##species", species_url)
    ),
    pattern = "failed to retrieve organism information from NCBI"
)
expect_identical(res, list())

gff <- tempfile(fileext = ".gff3")
writeLines(
    c(
        "##gff-version 3",
        paste("##species", species_url),
        "chr1\t.\tgene\t1\t10\t.\t+\t.\tID=g1"
    ),
    gff
)
expect_warning(
    gr <- BiocIO::import(gff),
    pattern = "no 'id' or 'name' query parameter"
)
expect_identical(length(gr), 1L)
expect_identical(S4Vectors::metadata(gr), list())

unlink(gff)
