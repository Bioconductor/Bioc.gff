.NCBI_TAX_URL <- Bioc.gff:::.NCBI_TAX_URL

Bioc.gff:::isNCBISpeciesURL(.NCBI_TAX_URL)

response <- Bioc.gff:::metadataFromNCBI(
    paste0(.NCBI_TAX_URL, "?mode=Info&id=9606")
)
expect_identical(
    response,
    list(
        `Taxonomy ID` = "9606", Organism = "homo sapiens", Rank = "species"
    )
)

response <- Bioc.gff:::metadataFromNCBI(
    paste0(.NCBI_TAX_URL, "?id=3702")
)
expect_identical(
    response,
    list(
        `Taxonomy ID` = "3702", Organism = "arabidopsis thaliana",
        Rank = "species"
    )
)

response <- Bioc.gff:::metadataFromNCBI(
    paste0(.NCBI_TAX_URL, "?name=drosophila")
)
expect_true(length(response[["Taxonomy ID"]]) > 1L)
expect_identical(
    length(response[["Taxonomy ID"]]), length(response[["Organism"]])
)
expect_identical(
    length(response[["Taxonomy ID"]]), length(response[["Rank"]])
)
expect_identical(
    response[["Rank"]][response[["Taxonomy ID"]] == "7215"], "genus"
)
expect_true(all(response[["Organism"]] == "drosophila"))
