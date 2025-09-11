.NCBI_TAX_URL <- Bioc.gff:::.NCBI_TAX_URL

Bioc.gff:::isNCBISpeciesURL(.NCBI_TAX_URL)

response <- Bioc.gff:::metadataFromNCBI(
    paste0(.NCBI_TAX_URL, "?mode=Info&id=9606")
)
expect_identical(
    response,
    list(`Taxonomy ID` = "9606", Organism = "homo sapiens")
)

response <- Bioc.gff:::metadataFromNCBI(
    paste0(.NCBI_TAX_URL, "?id=3702")
)
expect_identical(
    response,
    list(`Taxonomy ID` = "3702", Organism = "arabidopsis thaliana")
)
