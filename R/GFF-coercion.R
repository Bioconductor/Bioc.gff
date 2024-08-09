# Coercion ----------------------------------------------------------------

#' @name asGFF
#'
#' @title Coerce to GFF structure
#'
#' @description Coerce the structure of an object to one following GFF-like
#'   conventions, i.e., using the `Parent` GFF3 attribute to encode the
#'   hierarchical structure. This object is then suitable for export as GFF3.
#'
#' @param x Generally, a tabular object to structure as GFF(3)
#'
#' @param parentType The value to store in the `type` column for the
#' top-level (e.g., transcript) ranges.
#'
#' @param childType The value to store in the `type` column for the child
#' (e.g., exon) ranges.
#'
#' @param \dots Arguments to pass to methods
#'
#' @return For the `GRangesList` method: A `GRanges`, with the columns: `ID`
#'   (unique identifier), `Name` (from `names(x)`, and the names on each
#'   element of `x`, if any), `type` (as given by `parentType` and `childType`),
#'   and `Parent` (to relate each child range to its parent at the top-level).
#'
#' @author Michael Lawrence
#'
#' @examples
#' \dontrun{
#'   library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#'   library(GenomicFeatures)
#'   exons <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene)
#'   mcols(asGFF(exons))
#' }
#' @export
setGeneric("asGFF", function(x, ...) standardGeneric("asGFF"))

#' @describeIn asGFF Coerce to GFF GRanges structure
#'
#' @exportMethod asGFF
setMethod("asGFF", "GRangesList",
          function(x, parentType = "mRNA", childType = "exon") {
            parent_range <- range(x)
            if (!all(elementNROWS(parent_range) == 1))
              stop("Elements in a group must be on same sequence and strand")
            parents <- unlist(parent_range, use.names = FALSE)
            children <- unlist(x, use.names = FALSE)
            makeId <- function(x, prefix) {
                paste(prefix, seq_len(length(x)), sep = "")
            }
            parentIds <- makeId(parents, parentType)
            values(parents)$type <- parentType
            values(parents)$ID <- parentIds
            values(parents)$Name <- names(x)
            values(children)$type <- childType
            values(children)$ID <- makeId(children, childType)
            values(children)$Name <- names(children)
            values(children)$Parent <- rep.int(parentIds, elementNROWS(x))
            allColumns <- union(colnames(values(parents)),
                                colnames(values(children)))
            values(children) <- rectifyDataFrame(values(children), allColumns)
            values(parents) <- rectifyDataFrame(values(parents), allColumns)
            c(parents, children)
          })

#' @importFrom S4Vectors DataFrame
rectifyDataFrame <- function(x, allColumns) {
    x[setdiff(allColumns, colnames(x))] <- DataFrame(NA)
    x[allColumns]
}

### FIXME: We wrote this but never tested it, and it is not yet
### used. People should use GFF3 instead of this.
###
### KNOWN ISSUES:
### 1) The stop codon should not be included in the CDS (annoying)
### 2) pmapFromTranscripts() does not yet support our usage of it
### 3) Needs to move to GenomicFeatures

setGeneric("asGTF", function(x, ...) standardGeneric("asGTF"))

#' @importFrom BiocGenerics width end
#' @importFrom IRanges PartitioningByEnd
#' @importFrom utils head
frame <- function(x) {
    cs <- cumsum(width(x))
    ucs <- unlist(cs, use.names=FALSE)
    ucs[end(PartitioningByEnd(x))] <- 0L
    ucs <- c(0L, head(ucs, -1L))
    ucs %% 3L
}

#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges mcols mcols<-
setMethod("asGTF", "GRangesList",
          function(x) {
              tx_ids <- names(x)
              if (is.null(tx_ids)) {
                  tx_ids <- seq_along(x)
              }
              processFeatures <- function(f) {
                  ans <- unlist(f, use.names=FALSE)
                  ans$frame <- frame(f)
                  if (is.null(ans$gene_id)) {
                      ans$gene_id <- ""
                  }
                  if (is.null(ans$transcript_id)) {
                      ans$transcript_id <- tx_ids[togroup(f)]
                  }
                  ans
              }
              start_codon_tx <-
                  GenomicFeatures::pmapFromTranscripts(IRanges(1L, 3L), x)
              start_codon <- processFeatures(start_codon_tx)
              mcols(start_codon)$type <- "start_codon"
              stop_ranges <- IRanges(end=sum(width(x)), width=3L)
              stop_codon_tx <-
                  GenomicFeatures::pmapFromTranscripts(stop_ranges, x)
              stop_codon <- processFeatures(stop_codon_tx)
              mcols(stop_codon)$type <- "stop_codon"
              codons <- c(start_codon, stop_codon)
              cds <- processFeatures(x)
              mcols(cds)$type <- "CDS"
              values(codons) <- rectifyDataFrame(values(codons), colnames(cds))
              c(codons, cds)
          })

## setMethod("asGTF", "TxDb",
##           function(x, by) {
##               cds <- cds(x, columns="tx_id")
##               cds <- cds[togroup(cds$tx_id)]
##               cds$transcript_id <- unlist(cds$tx_id)
##               cds$tx_id <- NULL
##               txGene <- transcriptsBy(x)
##               txToGene <- setNames(names(txGene)[togroup(txGene)],
##                                    unlist(txGene, use.names=FALSE)$tx_id)
##               cds$gene_id <- txToGene[cds$transcript_id]
##               processUTRs <- function(utr) {
##                   ans <- unlist(utr, use.names=FALSE)
##                   ans$transcript_id <- names(utr)[togroup(utr)]
##                   ans$gene_id <- txToGene[ans$transcript_id]
##                   ans
##               }
##               three_utr <- processUTRs(threeUTRsByTranscript(x))
##               three_utr$type <- "3UTR"
##               five_utr <- processUTRs(fiveUTRsByTranscript(x))
##               five_utr$type <- "5UTR"
##               c(asGTF(split(cds, cds$transcript_id)), five_utr, three_utr)
##           })
