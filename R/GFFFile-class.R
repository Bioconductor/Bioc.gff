#' @include readGFF.R

# GFF (General Feature Format) support (all three versions, plus G --------

#' @name GFFFile-class
#'
#' @aliases class:GFFFile class:GFF1File class:GFF2File class:GFF3File
#' class:GVFFile class:GTFFile GFFFile-class GFF1File-class GFF2File-class
#' GFF3File-class GVFFile-class GTFFile-class GFFFile GFF1File GFF2File
#' GFF3File GVFFile GTFFile import,GFFFile,ANY,ANY-method import.gff
#' import.gff1 import.gff2 import.gff3 import.gff,ANY-method
#' import.gff1,ANY-method import.gff2,ANY-method import.gff3,ANY-method
#' export,ANY,GFFFile,ANY-method export,GenomicRanges,GFFFile,ANY-method
#' export,GenomicRangesList,GFFFile,ANY-method
#' export,GRangesList,GFFFile,ANY-method export,GRangesList,GTFFile,ANY-method
#' export.gff export.gff,ANY-method export.gff1 export.gff1,ANY-method
#' export.gff2 export.gff2,ANY-method export.gff3 export.gff3,ANY-method
#' genome,GFFFile-method
#'
#' @docType class
#'
#' @title GFFFile objects
#'
#' @description These functions support the import and export of the GFF format,
#'   of which there are three versions and several flavors.
#'
#' @details The Generic Feature Format (GFF) format is a tab-separated table of
#' intervals. There are three different versions of GFF, and they all have the
#' same number of columns. In GFF1, the last column is a grouping factor,
#' whereas in the later versions the last column holds application-specific
#' attributes, with some conventions defined for those commonly used. This
#' attribute support facilitates specifying extensions to the format. These
#' include GTF (Gene Transfer Format, an extension of GFF2) and GVF (Genome
#' Variation Format, an extension of GFF3).  The rtracklayer package recognizes
#' the "gtf" and "gvf" extensions and parses the extra attributes
#' into columns of the result; however, it does not perform any
#' extension-specific processing. Both GFF1 and GFF2 have been proclaimed
#' obsolete; however, the UCSC Genome Browser only supports GFF1 (and GTF), and
#' GFF2 is still in broad use.
#'
#' GFF is distinguished from the simpler BED format by its flexible attribute
#' support and its hierarchical structure, as specified by the `group`
#' column in GFF1 (only one level of grouping) and the `Parent` attribute
#' in GFF3. GFF2 does not specify a convention for representing hierarchies,
#' although its GTF extension provides this for gene structures. The
#' combination of support for hierarchical data and arbitrary descriptive
#' attributes makes GFF(3) the preferred format for representing gene models.
#'
#' Although GFF features a `score` column, large quantitative data belong
#' in a format like \link[=BigWigFile]{BigWig} and alignments from
#' high-throughput experiments belong in \link[Rsamtools:BamFile]{BAM}. For
#' variants, the VCF format (supported by the VariantAnnotation package) seems
#' to be more widely adopted than the GVF extension.
#'
#' A note on the UCSC track line metaformat: track lines are a means for
#' passing hints to visualization tools like the UCSC Genome Browser and the
#' Integrated Genome Browser (IGB), and they allow multiple tracks to be
#' concatenated in the same file. Since GFF is not a UCSC format, it is not
#' common to annotate GFF data with track lines, but rtracklayer still supports
#' it. To export or import GFF data in the track line format, call
#' \code{\link{export.ucsc}} or \code{\link{import.ucsc}}.
#'
#' The following is the mapping of GFF elements to a `GRanges` object.  NA
#' values are allowed only where indicated.  These appear as a "." in
#' the file. GFF requires that all columns are included, so `export`
#' generates defaults for missing columns.
#'
#' \describe{
#' \item{seqid, start, end}{the `ranges` component.}
#' \item{source}{character vector in the `source` column; defaults to
#'   "rtracklayer" on export.}
#' \item{type}{character vector in the `type` column; defaults to
#'   "sequence_feature" in the output, i.e., SO:0000110.}
#' \item{score}{numeric vector (NA's allowed) in the `score` column,
#'   accessible via the `score` accessor; defaults to `NA` upon export.}
#' \item{strand}{strand factor (NA's allowed) in the `strand` column,
#'   accessible via the `strand` accessor; defaults to `NA` upon export.}
#' \item{phase}{integer vector, either 0, 1 or 2 (NA's allowed); defaults to
#'   `NA` upon export.}
#' \item{group}{a factor (GFF1 only); defaults to the `seqid` (e.g.,
#'   chromosome) on export.}
#' }
#'
#' In GFF versions 2 and 3, attributes map to arbitrary columns in the result.
#' In GFF3, some attributes (`Parent`, `Alias`, `Note`, `DBxref` and
#' `Ontology_term`) can have multiple, comma-separated values; these columns are
#' thus always `CharacterList` objects.
#'
#' @param con A path, URL, connection or `GFFFile` object. For the functions
#'   ending in `.gff`, `.gff1`, etc, the file format is indicated by the
#'   function name. For the base `export` and `import` functions, the format
#'   must be indicated another way. If `con` is a path, URL or connection,
#'   either the file extension or the `format` argument needs to be one of
#'   "gff", "gff1" "gff2", "gff3", "gvf", or "gtf". Compressed files ("gz",
#'   "bz2" and "xz") are handled transparently.
#'
#' @param object The object to export, should be a `GRanges` or something
#'   coercible to a `GRanges`. If the object has a method for `asGFF`, it is
#'   called prior to coercion. This makes it possible to export a `GRangesList`
#'   or `TxDb` in a way that preserves the hierarchical structure. For exporting
#'   multiple tracks, in the UCSC track line metaformat, pass a
#'   `GenomicRangesList`, or something coercible to one.
#'
#' @param format If not missing, should be one of "gff", "gff1" "gff2", "gff3",
#'   "gvf", or "gtf".
#'
#' @param version If the format is given as "gff", i.e., it does not specify a
#'   version, then this should indicate the GFF version as one of \dQuote{} (for
#'   import only, from the `gff-version` directive in the file or "1" if none),
#'   "1", "2" or "3".
#'
#' @param text If `con` is missing, a character vector to use as the input.
#'
#' @param genome The identifier of a genome, or a `Seqinfo`, or `NA` if unknown.
#'   Typically, this is a UCSC identifier like "hg19". An attempt will be made
#'   to derive the `Seqinfo` on the return value using either an installed
#'   BSgenome package or UCSC, if network access is available.
#'
#' @param colnames A character vector naming the columns to parse. These should
#'   name either fixed fields, like `source` or `type`, or, for GFF2 and GFF3,
#'   any attribute.
#'
#' @param which A `GRanges` or other range-based object supported by
#'   \code{\link[IRanges]{findOverlaps}}. Only the intervals in the file
#'   overlapping the given ranges are returned. This is much more efficient when
#'   the file is indexed with the tabix utility.
#'
#' @param feature.type `NULL` (the default) or a character vector of valid
#'   feature types. If not `NULL`, then only the features of the specified
#'   type(s) are imported.
#'
#' @param sequenceRegionsAsSeqinfo If `TRUE`, attempt to infer the `Seqinfo`
#'   (`seqlevels` and `seqlengths`) from the \dQuote{##sequence-region}
#'   directives as specified by GFF3.
#'
#' @param source The value for the source column in GFF. This is typically the
#'   name of the package or algorithm that generated the feature.
#'
#' @param index If `TRUE`, automatically compress and index the output file with
#'   bgzf and tabix. Note that tabix indexing will sort the data by chromosome
#'   and start. Tabix supports a single track in a file.
#'
#' @param append If `TRUE`, and `con` points to a file path, the data is
#'   appended to the file. Obviously, if `con` is a connection, the data is
#'   always appended.
#'
#' @param ... Arguments to pass down to methods to other methods. For import,
#'   the flow eventually reaches the `GFFFile` method on `import`. When
#'   `trackLine` is `TRUE` or the target format is BED15, the arguments are
#'   passed through `export.ucsc`, so track line parameters are supported.
#'
#' @param x A `GFFFile` object.
#'
#' @return A `GRanges` with the metadata columns described in the details.
#'
#' @section GFFFile objects: The `GFFFile` class extends
#'   \code{\link[BiocIO:BiocFile-class]{BiocFile}} and is a formal
#'   representation of a resource in the GFF format.  To cast a path, URL or
#'   connection to a `GFFFile`, pass it to the `GFFFile` constructor. The
#'   `GFF1File`, `GFF2File`, `GFF3File`, `GVFFile` and `GTFFile` classes all
#'   extend `GFFFile` and indicate a particular version of the format.
#'
#' @author Michael Lawrence
#'
#' @references
#' * GFF1, GFF2: <http://www.sanger.ac.uk/resources/software/gff/spec.html>
#' * GFF3: <http://www.sequenceontology.org/gff3.shtml>
#' * GVF: <http://www.sequenceontology.org/resources/gvf.html>
#' * GTF: <http://mblab.wustl.edu/GTF22.html>
#'
#' @keywords methods classes
#' @examples
#'
#'   test_path <- system.file("tests", package = "rtracklayer")
#'   test_gff3 <- file.path(test_path, "genes.gff3")
#'
#'   ## basic import
#'   test <- import(test_gff3)
#'   test
#'
#'   ## import.gff functions
#'   import.gff(test_gff3)
#'   import.gff3(test_gff3)
#'
#'   ## GFFFile derivatives
#'   test_gff_file <- GFF3File(test_gff3)
#'   import(test_gff_file)
#'   test_gff_file <- GFFFile(test_gff3)
#'   import(test_gff_file)
#'   test_gff_file <- GFFFile(test_gff3, version = "3")
#'   import(test_gff_file)
#'
#'   ## from connection
#'   test_gff_con <- file(test_gff3)
#'   test <- import(test_gff_con, format = "gff")
#'
#'   ## various arguments
#'   import(test_gff3, genome = "hg19")
#'   import(test_gff3, colnames = character())
#'   import(test_gff3, colnames = c("type", "geneName"))
#'
#'   ## 'which'
#'   which <- GRanges("chr10:90000-93000")
#'   import(test_gff3, which = which)
#'
#' \dontrun{
#'   ## 'append'
#'   test_gff3_out <- file.path(tempdir(), "genes.gff3")
#'
#'   export(test[seqnames(test) == "chr10"], test_gff3_out)
#'   export(test[seqnames(test) == "chr12"], test_gff3_out, append = TRUE)
#'   import(test_gff3_out)
#'
#'   ## 'index'
#'   export(test, test_gff3_out, index = TRUE)
#'   test_bed_gz <- paste(test_gff3_out, ".gz", sep = "")
#'   import(test_bed_gz, which = which)
#' }
#'
NULL

# Classes -----------------------------------------------------------------

#' @rdname GFFFile-class
#'
#' @importClassesFrom BiocIO BiocFile
#' @exportClass GFFFile
setClass("GFFFile", contains = "BiocFile")

#' @rdname GFFFile-class
#'
#' @importFrom methods new
#' @param resource `character(1)` or `connection` A low-level resource typically
#'   a path, URL, or connection.
#'
#' @export
GFFFile <- function(resource, version = c("", "1", "2", "3")) {
  version <- match.arg(version)
  new(gffFileClass(version), resource = resource)
}

#' @exportClass GFF1File
setClass("GFF1File", contains = "GFFFile")

#' @export
GFF1File <- function(resource) {
  GFFFile(resource, "1")
}

#' @exportClass GFF2File
setClass("GFF2File", contains = "GFFFile")

#' @export
GFF2File <- function(resource) {
  GFFFile(resource, "2")
}

#' @exportClass GFF3File
setClass("GFF3File", contains = "GFFFile")

#' @export
GFF3File <- function(resource) {
  GFFFile(resource, "3")
}

#' @exportClass GTFFile
setClass("GTFFile", contains = "GFF2File")

#' @export
GTFFile <- function(resource) {
  new("GTFFile", GFF2File(resource))
}

#' @exportClass GVFFile
setClass("GVFFile", contains = "GFF3File")

#' @export
GVFFile <- function(resource) {
  new("GVFFile", GFF3File(resource))
}

# Export ------------------------------------------------------------------

#' @describeIn GFFFile-class
#'
#' @export
setGeneric("export.gff",
           function(object, con, ...) standardGeneric("export.gff"))

#' @describeIn GFFFile-class
#'
#' @exportMethod export.gff
setMethod("export.gff", "ANY",
          function(object, con, ...)
          {
            export(object, con, ...)
          })

#' @describeIn GFFFile-class
#'
#' @importFrom BiocIO export
#' @importFrom methods is as hasMethod
#' @exportMethod export
setMethod("export", c("ANY", "GFFFile"),
          function(object, con, format, ...)
          {
            if (hasMethod("asGFF", class(object)))
              object <- asGFF(object)
            res <- try(as(object, "GRanges"), silent = TRUE)
            if (is(res, "try-error")) {
              res <- try(as(object, "GenomicRangedDataList"), silent = TRUE)
              if (is(res, "try-error"))
                stop("cannot export object of class '", class(object), "'")
            }
            object <- res
            if (!missing(format))
              checkArgFormat(con, format)
            export(object, con, ...)
          })

#' @describeIn GFFFile-class
#'
#' @importClassesFrom GenomicRanges CompressedGRangesList GRangesList
#'   GenomicRanges
#' @importFrom methods callGeneric
#'
#' @exportMethod export
setMethod("export", c("CompressedGRangesList", "GFFFile"),
          function(object, con, format, ...)
          {
            object <- asGFF(object)
            callGeneric()
          }
          )

#' @exportMethod export
setMethod("export", c("GRangesList", "GTFFile"),
          function(object, con, format, ...) {
              stop("export of GRangesList to GTF is not yet supported")
              ## there is a start on asGTF() later in this file
              ## object <- asGTF(object)
              ## callGeneric()
          }
          )

#' @describeIn GFFFile-class
#'
#' @importFrom BiocIO resource
#' @importFrom GenomeInfoDb genome seqnames
#' @importFrom BiocGenerics start
#' @importFrom utils packageVersion relist write.table
#' @importFrom S4Vectors wmsg
#'
#' @exportMethod export
setMethod("export", c("GenomicRanges", "GFFFile"),
          function(object, con, format, version = c("1", "2", "3"),
                   source = "rtracklayer", append = FALSE, index = FALSE)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            if (!missing(version) || !length(gffFileVersion(con)))
              con <- asGFFVersion(con, match.arg(version))
            version <- gffFileVersion(con)

            file <- con
            con <- resource(con)

            if (!append) {
              cat("", file = con) # clear any existing file
              gffComment(con, "gff-version", version)
              sourceVersion <- try(packageVersion(source), TRUE)
              if (!inherits(sourceVersion, "try-error"))
                gffComment(con, "source-version", source, sourceVersion)
              gffComment(con, "date", base::format(Sys.time(), "%Y-%m-%d"))
              genome <- singleGenome(genome(object))
              if (!is.na(genome))
                gffComment(con, "genome-build", paste(".", genome, sep = "\t"))
            }

            if (index)
              object <- .sortBySeqnameAndStart(object)

            seqname <- seqnames(object)
            if (is.null(mcols(object)$ID))
              mcols(object)$ID <- names(object)
            if (version == "3")
              seqname <- urlEncode(seqname, "a-zA-Z0-9.:^*$@!+_?|-")
            if (!is.null(mcols(object)$source) && missing(source))
              source <- mcols(object)$source
            else source <- rep(source, length(object))
            if (version == "3")
              source <- urlEncode(source, "\t\n\r;=%&,", FALSE)
            feature <- mcols(object)$type
            if (is.null(feature))
              feature <- rep("sequence_feature", length(object))
            score <- score(object)
            if (is.null(score)) {
              score <- rep(NA_real_, length(object))
            } else {
              if (!("score" %in% colnames(mcols(object))))
                ## avoid outputting as attribute
                colnames(mcols(object))[1] <- "score"
            }
            strand <- strand(object)
            if (is.null(strand))
                strand <- rep(strand(NA_character_), length(object))
            strand[strand == "*"] <- NA_integer_
            frame <- mcols(object)$phase
            if (is.null(frame)) {
              frame <- rep(NA_integer_, length(object))
              if ("CDS" %in% feature)
                warning(wmsg("The phase information is missing. ",
                             "The written file will contain CDS with ",
                             "no phase information."))
            } else {
              if (anyNA(frame[feature %in% "CDS"]))
                warning(wmsg("The phase information is missing for some CDS. ",
                             "The written file will contain some CDS with ",
                             "no phase information."))
            }

            table <- data.frame(seqname, source, feature, start(object),
                                end(object), score, strand, frame)

            attrs <- NULL
            if (version == "1") {
              attrs <- mcols(object)$group
              if (is.null(attrs))
                attrs <- as.vector(seqname)
            } else {
              builtin <- c("type", "score", "phase", "source")
              custom <- setdiff(colnames(mcols(object)), builtin)
              if (length(custom)) {
                if (version == "3") tvsep <- "=" else tvsep <- " "
                attrs <- mcols(object)
                attrs <- as.data.frame(sapply(custom, function(name) {
                  x <- attrs[[name]]
                  x_flat <- if (is(x, "List")) unlist(x, use.names=FALSE) else x
                  x_char <- as.character(x_flat)
                  x_char <- sub(" *$", "", sub("^ *", "", as.character(x_char)))
                  if (version == "3")
                    x_char <- urlEncode(x_char, "%\t\n\r;=&,", FALSE)
                  if (is(x, "List")) {
                    x_char[is.na(x_char)] <- "."
                    x_char <- pasteCollapse(relist(x_char, x))
                    x_char[elementNROWS(x) == 0] <- NA
                  }
                  ## FIXME: add option so these become "." instead of removing
                  x_char[is.na(x_char)] <- "\r"
                  if (!is.numeric(x_flat) && version != "3")
                    x_char <- paste0("\"", x_char, "\"")
                  paste(name, x_char, sep = tvsep)
                }, simplify = FALSE))
                if (version == "3") sep <- ";" else sep <- "; "
                attrs <- do.call(paste, c(attrs, sep = sep))
                attrs <- gsub("[^;]*?\r\"?(;|$)", "", attrs)
                attrs[nchar(attrs) == 0] <- NA
              }
            }

            scipen <- getOption("scipen")
            options(scipen = 100) # prevent use of scientific notation
            on.exit(options(scipen = scipen))

            if (!is.null(attrs)) { # write out the rows with attributes first
              write.table(cbind(table, attrs)[!is.na(attrs),], con, sep = "\t",
                          na = ".", quote = FALSE, col.names = FALSE,
                          row.names = FALSE, append = TRUE)
              table <- table[is.na(attrs),]
            }

            write.table(table, con, sep = "\t", na = ".", quote = FALSE,
                        col.names = FALSE, row.names = FALSE, append = TRUE)
            if (index)
              tracklayer:::indexTrack(file)
            invisible(NULL)
          })

setClass("UCSCFile", contains = "BiocFile")

UCSCFile <- function(resource) {
  new("UCSCFile", resource = resource)
}

#' @importFrom BiocIO fileFormat
.export_SimpleGRangesList_BiocFile <- function(object, con, format, ...) {
  export(object, UCSCFile(resource(con)), subformat = fileFormat(con), ...)
}

#' @describeIn GFFFile-class
#'
#' @exportMethod export
setMethod("export", c("SimpleGRangesList", "GFFFile"),
          .export_SimpleGRangesList_BiocFile)

#' @describeIn GFFFile-class
#'
#' @export
setGeneric("export.gff1",
           function(object, con, ...) standardGeneric("export.gff1"))

#' @describeIn GFFFile-class
#'
#' @exportMethod export.gff1
setMethod("export.gff1", "ANY",
          function(object, con, ...) export(object, con, "gff1", ...))

#' @describeIn GFFFile-class
#'
#' @export
setGeneric("export.gff2",
           function(object, con, ...) standardGeneric("export.gff2"))

#' @describeIn GFFFile-class
#'
#' @exportMethod export.gff2
setMethod("export.gff2", "ANY",
          function(object, con, ...) export(object, con, "gff2", ...))

#' @describeIn GFFFile-class
#'
#' @export
setGeneric("export.gff3",
           function(object, con, ...) standardGeneric("export.gff3"))

#' @describeIn GFFFile-class
#'
#' @exportMethod export.gff3
setMethod("export.gff3", "ANY",
          function(object, con, ...) export(object, con, "gff3", ...))

# Import ------------------------------------------------------------------

#' @export
setGeneric("import.gff", function(con, ...) standardGeneric("import.gff"))

#' @exportMethod import.gff
setMethod("import.gff", "ANY",
          function(con, ...)
          {
            import(con, "gff", ...)
          })

connection <- BiocIO:::connection
manager <- BiocIO:::manager
release <- BiocIO:::release

#' @describeIn GFFFile-class
#'
#' @importFrom BiocIO import
#' @importFrom utils download.file
#' @importFrom IRanges subsetByOverlaps
#' @importFrom S4Vectors isTRUEorFALSE
#'
#' @exportMethod import
setMethod("import", "GFFFile",
          function(con, format, text, version = c("", "1", "2", "3"),
                   genome = NA, colnames = NULL,
                   which = NULL, feature.type = NULL,
                   sequenceRegionsAsSeqinfo = FALSE)
          {
            if (!missing(format))
              checkArgFormat(con, format)
            if (!missing(version))
              con <- asGFFVersion(con, match.arg(version))
            stopifnot(isTRUEorFALSE(sequenceRegionsAsSeqinfo))

            ## download the file first if it's remote
            if (is.character(resource(con))) {
                uri <- BiocIO:::.parseURI(resource(con))
                if (uri$scheme %in% c("ftp", "http")) {
                    destfile <- tempfile()
                    download.file(resource(con), destfile)
                    con@resource <- destfile
                }
            }

            m <- manager()
            sniff_con <- connection(m, con, "r")
            on.exit(release(m, sniff_con))
            sniffed <- .sniffGFFVersion(sniff_con)
            version <- gffFileVersion(con)
            if (!length(version)) {
              if (is.null(sniffed))
                sniffed <- "1"
              con <- asGFFVersion(con, sniffed)
            }

            if (length(version) && !is.null(sniffed) &&
                !identical(sniffed, version))
              warning("gff-version directive indicates version is ", sniffed,
                      ", not ", version)

            if (is(genome, "Seqinfo") && length(genome) == 0L) {
                genome <- NA_character_
            }

            if (is.na(genome)) {
              genome <- genome(con)
              if (is.null(genome))
                  genome <- NA
            }

### FIXME: a queryForLines() function would be more efficient

            ## Temporarily disable use of Tabix Index.
            ## TODO: Restore use of Tabix Index!
            #con <- queryForResource(m, con, which)
            resource <- queryForResource(m, con)
            on.exit(release(m, resource), add=TRUE)
            ans <- readGFFAsGRanges(resource,
                                    version=version,
                                    colnames=colnames,
                                    filter=list(type=feature.type),
                                    genome=genome,
                                    sequenceRegionsAsSeqinfo=
                                        sequenceRegionsAsSeqinfo,
                                    speciesAsMetadata=TRUE)
            if (!attr(resource, "usedWhich") && !is.null(which))
                ans <- subsetByOverlaps(ans, which)
            ans
          })

#' @describeIn GFFFile-class
#'
#' @export
setGeneric("import.gff1",
           function(con, ...) standardGeneric("import.gff1"))

#' @describeIn GFFFile-class
#'
#' @exportMethod import.gff1
setMethod("import.gff1", "ANY",
          function(con, ...) import(con, "gff1", ...))

#' @describeIn GFFFile-class
#'
#' @export
setGeneric("import.gff2",
           function(con, ...) standardGeneric("import.gff2"))

#' @describeIn GFFFile-class
#'
#' @exportMethod import.gff2
setMethod("import.gff2", "ANY",
          function(con, ...) import(con, "gff2", ...))

#' @describeIn GFFFile-class
#'
#' @export
setGeneric("import.gff3",
           function(con, ...) standardGeneric("import.gff3"))

#' @describeIn GFFFile-class
#'
#' @exportMethod import.gff3
setMethod("import.gff3", "ANY",
          function(con, ...) import(con, "gff3", ...))
