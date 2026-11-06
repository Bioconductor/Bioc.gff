# Bioc.gff 1.4.0

## SIGNIFICANT USER-VISIBLE CHANGES

* The `##species` metadata added by `import()` now has a `Rank` element
  next to `Taxonomy ID` and `Organism`. A `##species` URL that gives an
  organism `name` is treated as a search term, and every matching taxon
  is returned.

## BUG FIXES AND MINOR IMPROVEMENTS

* Organism names are found again for the `##species` NCBI Taxonomy
  Browser URL. The Taxonomy Browser now sits behind a bot check, so the
  lookup uses the NCBI E-utilities API instead of scraping the page.
* `export()` to GFF and GTF files no longer fails with
  `could not find function "score"` when BiocGenerics is not attached;
  `score` and `strand` are now imported (#1).
* `import()` no longer fails with "replacement 'metadata' value must be a
  list" when the NCBI lookup for `##species` fails. It now warns, gives
  the reason, and returns empty metadata.
* Replace rvest with jsonlite in Suggests.
* Use roxygen2 8.1.1 and regenerate NAMESPACE.
* Add descriptions to `@describeIn` tags and use `@aliases` tags in the
  `GFFFile-class` documentation.
* Add tests for NCBI name searches and `Rank`, GTF export without
  BiocGenerics attached, and GFF import when the NCBI lookup fails.

# Bioc.gff 0.99.5

* Initial Bioconductor submission.
