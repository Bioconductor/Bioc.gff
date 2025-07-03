# Utilities ---------------------------------------------------------------

connection <- BiocIO:::connection
manager <- BiocIO:::manager
release <- BiocIO:::release

scanGFFDirectives <- function(con, tag = NULL) {
    m <- manager()
    con <- connection(m, con, "r")
    on.exit(release(m, con))
    directives <- character()
    lines <- line <- readLines(con, n = 1)
    while (grepl("^#", line)) {
        if (grepl("^##", line)) {
            directives <- c(directives, line)
        }
        line <- readLines(con, n = 1)
        if (length(line) == 0L) break
        lines <- c(lines, line)
    }
    pushBack(lines, con)
    sub(
        "^[^[:space:]]* ",
        "",
        grep(paste0("^##", tag), directives, value = TRUE)
    )
}

gffGenomeBuild <- function(x) {
    genome_build <- scanGFFDirectives(x, "genome-build")
    unlist(strsplit(genome_build, "\t", fixed = TRUE))
}

#' @importFrom GenomeInfoDb provider providerVersion genome
setMethod("provider", "GFFFile", function(x) {
    gffGenomeBuild(x)[1]
})

setMethod("providerVersion", "GFFFile", function(x) {
    gffGenomeBuild(x)[2]
})

#' @describeIn GFFFile-class Gets the genome identifier from the "genome-build"
#'   header directive.
#'
#' @exportMethod genome
setMethod("genome", "GFFFile", function(x) providerVersion(x))

gffComment <- function(con, ...)
    cat("##", paste(...), "\n", sep = "", file = con, append = TRUE)

.sniffGFFVersion <- function(con) {
    version <- NULL
    lines <- line <- readLines(con, n = 1)
    while (grepl("^#", line)) {
        if (grepl("^##gff-version", line)) {
            version <- sub("^##gff-version *", "", line)
            break
        }
        line <- readLines(con, n = 1)
        lines <- c(lines, line)
    }
    pushBack(lines, con)
    version
}

gffFileClass <- function(version) {
    paste("GFF", version, "File", sep = "")
}

gffFileVersion <- function(file) {
    versions <- c("1", "2", "3")
    unlist(Filter(function(v) is(file, gffFileClass(v)), versions))
}

asGFFVersion <- function(con, version) {
    if (!is(con, gffFileClass(version))) {
        if (!is(con, "GFFFile"))
            warning(
                "Treating a '",
                class(con),
                "' as GFF version '",
                version,
                "'"
            )
        con <- GFFFile(resource(con), version)
    }
    con
}

#' @importFrom S4Vectors isSingleString
pasteCollapse <- function(x, collapse = ",") {
    if (!is(x, "CharacterList")) stop("'x' must be a CharacterList")
    if (!isSingleString(collapse) || nchar(collapse) != 1L)
        stop("'collapse' must be a single string, with a single character")
    x <- as.list(x)
    .Call(CharacterList_pasteCollapse, x, collapse)
}

## taken from rtracklayer:::checkArgFormat
checkArgFormat <-
    function(con, format) {
        if (
            toupper(format) !=
                substring(
                    toupper(sub("File$", "", class(con))),
                    1,
                    nchar(format)
                )
        )
            stop("Cannot treat a '", class(con), "' as format '", format, "'")
    }

## adapted from rtracklayer:::singleGenome
#' @importFrom BiocBaseUtils isScalarCharacter
singleGenome <-
    function(x) {
        x1 <- unname(unique(x))
        if (!isScalarCharacter(x1, na.ok = TRUE))
            stop("Multiple genomes encountered; only one supported")
        x1
    }

## taken from rtracklayer:::urlEncode
urlEncode <-
    function(str, chars = "-a-zA-Z0-9$_.+!*'(),", keep = TRUE) {
        str <- as.character(str)
        bad <- chars
        if (keep) bad <- gsub(paste("[", chars, "]", sep = ""), "", str)
        bad <- unique(unlist(strsplit(bad, "")))
        code <- paste("%", charToRaw(paste(bad, collapse = "")), sep = "")
        for (i in seq_along(bad))
            str <- gsub(bad[i], code[i], str, fixed = TRUE)
        str
    }

## getMethod(rtracklayer:::sortBySeqnameAndStart, "GenomicRanges")
.sortBySeqnameAndStart <- function(x) {
    x[order(as.factor(seqnames(x)), start(x)), ]
}

## taken from rtracklayer:::uriIsLocal
uriIsLocal <- function(x) {
    x$scheme == "file" || x$scheme == ""
}

## taken from rtracklayer:::indexTrack
#' @importFrom Rsamtools bgzip indexTabix TabixFile
indexTrack <- function(con, ...) {
    indexed <- NULL
    formats <- eval(formals(indexTabix)$format)
    uri <- path(con)
    parsed_uri <- BiocIO:::.parseURI(uri)
    if (!uriIsLocal(parsed_uri)) stop("'con' must be a path to a local file")
    original_path <- parsed_uri$path
    path <- bgzip(original_path, overwrite = TRUE)
    format <- Find(
        function(f) {
            is(con, paste(toupper(f), "File", sep = ""))
        },
        formats
    )
    if (!is.null(format)) {
        indexTabix(path, format, ...)
    } else {
        indexTabix(path, ...)
    }
    indexed <- TabixFile(path)
    unlink(original_path)
    invisible(indexed)
}
