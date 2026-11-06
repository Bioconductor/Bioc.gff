.NCBI_TAX_URL <- "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi"

.NCBI_EUTILS_URL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

#' @name metadataFromNCBI
isNCBISpeciesURL <- function(url) {
    grepl("taxonomy/browser", url, ignore.case = TRUE)
}

#' Obtain metadata from NCBI
#'
#' These helper functions obtain both the Taxonomy ID and the Organism name
#' for a NCBI Taxonomy Browser URL. They are a modern re-write of the old
#' functions in `rtracklayer`. The Taxonomy Browser web pages are no longer
#' scrapable (they redirect to NCBI Datasets behind a bot check), so the
#' `id` or `name` query parameter is extracted from the URL and resolved
#' with the NCBI E-utilities API via `httr2`. A `name` is treated as a
#' search term and all matching Taxonomy IDs are returned.
#'
#' @param url A URL to the NCBI Taxonomy Browser, typically obtained from a
#' GFF file with the `## species` line.
#'
#' @param id `character()` One or more NCBI Taxonomy IDs.
#'
#' @param name `character(1)` An organism name, e.g., "drosophila miranda".
#'
#' @returns
#' * metadataFromNCBI: A `list` with three elements: `Taxonomy ID`,
#'   `Organism`, and `Rank`. All are `character` vectors of the same length,
#'   with one element per match; a name search may match more than one taxon,
#'   e.g., a genus and a subgenus with the same name.
#' * isNCBISpeciesURL: A `logical` indicating if the URL is from the NCBI
#'   Taxonomy Browser.
#' * taxonomyIDFromNCBI: A `character` vector of all matching Taxonomy IDs.
#' * taxonFromNCBI: A `list` with two `character` vectors, `Organism`
#'   (lower case) and `Rank` (e.g., "species", "genus"), in the same order as
#'   `id`.
#'
#' @keywords internal
#' @examples
#' \donttest{
#'   isNCBISpeciesURL(.NCBI_TAX_URL)
#'
#'   metadataFromNCBI(
#'       paste0(.NCBI_TAX_URL, "?mode=Info&id=9606")
#'   )
#'   metadataFromNCBI(
#'       paste0(.NCBI_TAX_URL, "?id=3702")
#'   )
#'   metadataFromNCBI(
#'       paste0(.NCBI_TAX_URL, "?name=drosophila+melanogaster")
#'   )
#'   metadataFromNCBI(
#'       paste0(.NCBI_TAX_URL, "?name=drosophila+miranda")
#'   )
#' }
#' @importFrom BiocBaseUtils checkInstalled
metadataFromNCBI <- function(url) {
    BiocBaseUtils::checkInstalled(c("httr2", "jsonlite"))
    query <- httr2::url_parse(url)[["query"]]
    if (!is.null(query[["id"]]))
        Taxonomy_ID <- query[["id"]]
    else if (!is.null(query[["name"]]))
        Taxonomy_ID <- taxonomyIDFromNCBI(
            gsub("+", " ", query[["name"]], fixed = TRUE)
        )
    else
        stop("NCBI Taxonomy URL has no 'id' or 'name' query parameter")
    c(list("Taxonomy ID" = Taxonomy_ID), taxonFromNCBI(Taxonomy_ID))
}

.eutils <- function(endpoint, ...) {
    httr2::request(.NCBI_EUTILS_URL) |>
        httr2::req_url_path_append(endpoint) |>
        httr2::req_url_query(db = "taxonomy", retmode = "json", ...) |>
        httr2::req_retry(max_tries = 3L) |>
        httr2::req_perform() |>
        httr2::resp_body_json()
}

#' @rdname metadataFromNCBI
taxonomyIDFromNCBI <- function(name) {
    res <- .eutils("esearch.fcgi", term = name, retmax = 10000L)
    res <- res[["esearchresult"]]
    ids <- unlist(res[["idlist"]])
    if (!length(ids))
        stop("no NCBI Taxonomy ID found for '", name, "'")
    ids
}

#' @rdname metadataFromNCBI
taxonFromNCBI <- function(id) {
    res <- .eutils("esummary.fcgi", id = paste(id, collapse = ","))
    res <- res[["result"]]
    field <- function(i, what) {
        value <- res[[i]][[what]]
        if (is.null(value) || !nzchar(value))
            stop("no NCBI ", what, " found for Taxonomy ID '", i, "'")
        value
    }
    list(
        "Organism" = tolower(vapply(
            id, field, character(1L), what = "scientificname",
            USE.NAMES = FALSE
        )),
        "Rank" = vapply(
            id, field, character(1L), what = "rank", USE.NAMES = FALSE
        )
    )
}
