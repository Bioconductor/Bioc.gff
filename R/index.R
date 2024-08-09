# Separate index support --------------------------------------------------

# queryForResource --------------------------------------------------------

## There is some question as to whether the import routines should
## recurse through import.tabix, or perform the query
## internally. Currently, we are taking the latter route, because it
## is consistent with (gzip) decoding and it allows the parsers more
## flexibility.

setGeneric("queryForResource",
           function(manager, x, which = NULL, ...)
               standardGeneric("queryForResource"),
           signature="x")

## Attaches 'usedWhich' attribute, an optimization hint indicating
## that subsetting by 'which' has been performed and is no longer
## necessary. Probably premature.

setMethod("queryForResource", "BiocFile", function(manager, x, which = NULL, ...)
{
    r <- resource(x)
    ans <- structure(r, usedWhich = FALSE)
    if (!is.null(which) && is.character(r)) {
        x_tbi <- paste(r, "tbi", sep = ".")
        if (file.exists(x_tbi))
            ans <- queryForResource(
                manager, Rsamtools::TabixFile(r), which = which, ...
            )
    }
    ans
})

connectionForResource <- BiocIO:::connectionForResource

#' @importClassesFrom Rsamtools TabixFile
setMethod("queryForResource", "TabixFile",
          function(manager, x, which, header = TRUE, ...) {
              tabixHeader <- headerTabix(x)
              si <- Seqinfo(tabixHeader$seqnames)
              if (is.null(which)) {
                  buffer <- connectionForResource(manager, path(x), "r")
                  if (!header)
                      readLines(buffer, tabixHeader$skip)
              } else {
                  buffer <- manage(manager, file())
                  if (header) {
                      skippedLines <- readLines(path(x), tabixHeader$skip)
                      writeLines(skippedLines, buffer)
                  }
                  lines <- unlist(scanTabix(x, param = which), use.names = FALSE)
                  writeLines(lines, buffer)
                  si <- merge(si, seqinfo(which))
              }
              structure(buffer, usedWhich = TRUE, seqinfo = si)
          })

queryForConnection <- function(manager, x, which = NULL, ...) {
    resource <- queryForResource(manager, x, which = which, ...)
    con <- connectionForResource(manager, resource, open = "r")
    structure(con, usedWhich = attr(resource, "usedWhich"))
}
