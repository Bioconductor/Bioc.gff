.NCBI_TAX_URL <- "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi"

#' @name metadataFromNCBI
isNCBISpeciesURL <- function(url) {
    grepl("taxonomy/browser", url, ignore.case = TRUE)
}

#' Obtain metadata from NCBI
#'
#' These helper functions obtain both the Taxonomy ID and the Organism name from
#' the NCBI Taxonomy Browser. They are a modern re-write of the old functions in
#' `rtracklayer`. They use `httr2` and `rvest` to parse the HTML content.
#'
#' @param url A URL to the NCBI Taxonomy Browser, typically obtained from a
#' GFF file with the `## species` line.
#'
#' @keywords internal
#' @examples
#' isNCBISpeciesURL(.NCBI_TAX_URL)
#'
#' metadataFromNCBI(
#'     paste0(.NCBI_TAX_URL, "?mode=Info&id=9606")
#' )
#' metadataFromNCBI(
#'     paste0(.NCBI_TAX_URL, "?id=3702")
#' )
#' metadataFromNCBI(
#'     paste0(.NCBI_TAX_URL, "?name=drosophila+melanogaster")
#' )
#' metadataFromNCBI(
#'     paste0(.NCBI_TAX_URL, "?name=drosophila+miranda")
#' )
metadataFromNCBI <- function(url) {
    res_html <- httr2::request(url) |>
        httr2::req_perform() |>
        httr2::resp_body_html()
    Organism <- parseOrganismFromNCBI(res_html)
    Taxonomy_ID <- parseTaxonomyIDFromNCBI(res_html, url)
    list(
        "Taxonomy ID" = Taxonomy_ID,
        "Organism" = Organism
    )
}

#' @rdname metadataFromNCBI
parseOrganismFromNCBI <- function(html) {
    html |>
        rvest::html_elements("h2") |>
        rvest::html_text() |>
        tolower()
}

#' @rdname metadataFromNCBI
parseTaxonomyIDFromNCBI <- function(html, url) {
    if (grepl("id=", url, fixed = TRUE)) {
        httr2::url_parse(url)[[c("query", "id")]]
    } else {
        html |>
            rvest::html_elements(xpath = "/html/body/form/table[4]") |>
            rvest::html_text() |>
            gsub(".*\n(Taxonomy ID: [0-9]*)\\s.*", "\\1", x = _) |>
            strsplit(":\\s+") |>
            unlist() |>
            utils::tail(1L)
    }
}
