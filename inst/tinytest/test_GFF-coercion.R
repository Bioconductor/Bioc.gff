library(GenomicRanges)
library(tinytest)

## Test asGFF on a simple GRangesList
grl <- GRangesList(
    tx1 = GRanges("chr1", IRanges(c(1, 10), c(5, 15)), "+"),
    tx2 = GRanges("chr1", IRanges(c(20, 30), c(25, 35)), "-")
)

gff <- asGFF(grl)

expect_true(inherits(gff, "GRanges"))
expect_equal(length(gff), 6)
expect_equal(sum(gff$type == "mRNA"), 2)
expect_equal(sum(gff$type == "exon"), 4)
parents <- gff[gff$type == "mRNA"]
children <- gff[gff$type == "exon"]
expect_true(all(children$Parent %in% parents$ID))

## Test error when elements are on different sequences
grl_error_seq <- GRangesList(
    tx1 = GRanges(c("chr1", "chr2"), IRanges(c(1, 10), c(5, 15)), "+")
)
expect_error(
    asGFF(grl_error_seq),
    "Elements in a group must be on same sequence and strand"
)

## Test error when elements are on different strands
grl_error_strand <- GRangesList(
    tx1 = GRanges("chr1", IRanges(c(1, 10), c(5, 15)), c("+", "-"))
)
expect_error(
    asGFF(grl_error_strand),
    "Elements in a group must be on same sequence and strand"
)
