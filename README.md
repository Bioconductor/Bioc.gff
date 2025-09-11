
# Bioc.gff

## Introduction

The *[Bioc.gff](https://bioconductor.org/packages/3.22/Bioc.gff)*
package provides support for the General Feature Format (GFF) file
format, which is widely used for representing genomic features and
annotations. This package allows users to read, write, and manipulate
GFF versions 1, 2 (GTF), and 3 in R, making it easier to work with
genomic data.

While the
*[rtracklayer](https://bioconductor.org/packages/3.22/rtracklayer)*
package offers robust GFF support, it is a large package with many
features beyond file import. `Bioc.gff` fills a specific niche in the
Bioconductor ecosystem by providing a lightweight, focused solution with
minimal dependencies. This modularity is beneficial for developers of
other packages who need to parse GFF files without inheriting the
extensive dependency footprint of `rtracklayer`. For the end user, it
offers a simple and direct interface for GFF manipulation.

Note that much of the code in the package is ported and adapted from the
*[rtracklayer](https://bioconductor.org/packages/3.22/rtracklayer)*
Bioconductor package. The intention is that the GFF functionality
residing in
*[rtracklayer](https://bioconductor.org/packages/3.22/rtracklayer)* will
be removed from that package in favor of using `Bioc.gff`, thereby
reducing its size and complexity.

# Installation

You can install the `Bioc.gff` package from Bioconductor using the
following:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Bioc.gff")
```

# Loading packages

``` r
library(Bioc.gff)
library(GenomicRanges)
```

# Usage

The `Bioc.gff` package provides functions to read and write GFF files,
as well as to manipulate GFF data structures. The following sections
provide some examples of how to use the package.

## Reading GFF Files

You can read GFF files using the `readGFF` function. This function
supports GFF versions 1, 2, and 3. The following example demonstrates
how to read a GFF file:

``` r
gff_file <- system.file("extdata", "genes.gff3", package = "Bioc.gff")
```

### Using the `import` function

The `import` function deduces that the input string is a GFF file type
and calls the appropriate method for reading the GFF file.

``` r
import(gff_file)
#> GRanges object with 31 ranges and 10 metadata columns:
#>        seqnames      ranges strand |      source     type     score     phase            ID        Name        geneName           Alias
#>           <Rle>   <IRanges>  <Rle> |    <factor> <factor> <numeric> <integer>   <character> <character>     <character> <CharacterList>
#>    [1]    chr10 92828-95504      - | rtracklayer     gene         5      <NA> GeneID:347688       TUBB8 tubulin, beta 8  FLJ40100,TUBB8
#>    [2]    chr10 92828-95178      - | rtracklayer     mRNA        NA      <NA>           873       TUBB8            <NA>                
#>    [3]    chr10 92828-95504      - | rtracklayer     mRNA        NA      <NA>           872       TUBB8            <NA>                
#>    [4]    chr10 92828-94054      - | rtracklayer     exon        NA      <NA>          <NA>        <NA>            <NA>                
#>    [5]    chr10 92997-94054      - | rtracklayer     CDS         NA      <NA>          <NA>        <NA>            <NA>                
#>    ...      ...         ...    ... .         ...      ...       ...       ...           ...         ...             ...             ...
#>   [27]    chr12 89675-89827      + | rtracklayer     CDS         NA      <NA>          <NA>        <NA>            <NA>                
#>   [28]    chr12 90587-90655      + | rtracklayer     exon        NA      <NA>          <NA>        <NA>            <NA>                
#>   [29]    chr12 90587-90655      + | rtracklayer     CDS         NA      <NA>          <NA>        <NA>            <NA>                
#>   [30]    chr12 90796-91263      + | rtracklayer     exon        NA      <NA>          <NA>        <NA>            <NA>                
#>   [31]    chr12 90796-91263      * | rtracklayer     CDS         NA      <NA>          <NA>        <NA>            <NA>                
#>             genome          Parent
#>        <character> <CharacterList>
#>    [1]        hg19                
#>    [2]        <NA>   GeneID:347688
#>    [3]        <NA>   GeneID:347688
#>    [4]        <NA>         872,873
#>    [5]        <NA>         872,873
#>    ...         ...             ...
#>   [27]        <NA>            4644
#>   [28]        <NA>            4644
#>   [29]        <NA>            4644
#>   [30]        <NA>            4644
#>   [31]        <NA>            4644
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

Here the output object is a `GRanges` class, which contains genomic
coordinates and associated metadata. The `GRanges` object indicates the
ranges with the `seqnames`, `ranges`, and `strand`, columns. The
associated metadata is stored in the `mcols` slot, which (in this
example) includes attributes like `source`, `type`, `phase`, `ID`,
`Name`, `geneName`, `Alias` and `genome`.

#### Selectively importing ranges

You can also selectively import ranges from the GFF file using the
`which` argument. This allows you to filter the data based on specific
criteria, such as chromosome or strand. For example, to import only the
ranges on chromosome 1:

``` r
which <- GRanges("chr10:90000-93000")
import(gff_file, which = which, genome = "hg19")
#> 
#> Attaching package: 'rtracklayer'
#> The following objects are masked from 'package:Bioc.gff':
#> 
#>     asGFF, export.gff, export.gff1, export.gff2, export.gff3, GFF1File, GFF2File, GFF3File, GFFcolnames, GFFFile, GTFFile,
#>     GVFFile, import.gff, import.gff1, import.gff2, import.gff3, readGFF
#> GRanges object with 5 ranges and 10 metadata columns:
#>       seqnames      ranges strand |      source     type     score     phase            ID        Name        geneName           Alias
#>          <Rle>   <IRanges>  <Rle> |    <factor> <factor> <numeric> <integer>   <character> <character>     <character> <CharacterList>
#>   [1]    chr10 92828-95504      - | rtracklayer     gene         5      <NA> GeneID:347688       TUBB8 tubulin, beta 8  FLJ40100,TUBB8
#>   [2]    chr10 92828-95178      - | rtracklayer     mRNA        NA      <NA>           873       TUBB8            <NA>                
#>   [3]    chr10 92828-95504      - | rtracklayer     mRNA        NA      <NA>           872       TUBB8            <NA>                
#>   [4]    chr10 92828-94054      - | rtracklayer     exon        NA      <NA>          <NA>        <NA>            <NA>                
#>   [5]    chr10 92997-94054      - | rtracklayer     CDS         NA      <NA>          <NA>        <NA>            <NA>                
#>            genome          Parent
#>       <character> <CharacterList>
#>   [1]        hg19                
#>   [2]        <NA>   GeneID:347688
#>   [3]        <NA>   GeneID:347688
#>   [4]        <NA>         872,873
#>   [5]        <NA>         872,873
#>   -------
#>   seqinfo: 298 sequences (2 circular) from hg19 genome
```

Note that you can indicate the genome build using the `genome` argument,
which is useful for ensuring that the imported data is compatible with
other genomic data you may be working with.

### Using the `readGFF` function

As an alternative, you can use the `readGFF` function directly, which is
more explicit about the GFF version:

``` r
readGFF(gff_file, version = "3")
#> DataFrame with 31 rows and 14 columns
#>        seqid      source     type     start       end     score      strand     phase            ID        Name        geneName
#>     <factor>    <factor> <factor> <integer> <integer> <numeric> <character> <integer>   <character> <character>     <character>
#> 1      chr10 rtracklayer     gene     92828     95504         5           -        NA GeneID:347688       TUBB8 tubulin, beta 8
#> 2      chr10 rtracklayer     mRNA     92828     95178        NA           -        NA           873       TUBB8              NA
#> 3      chr10 rtracklayer     mRNA     92828     95504        NA           -        NA           872       TUBB8              NA
#> 4      chr10 rtracklayer     exon     92828     94054        NA           -        NA            NA          NA              NA
#> 5      chr10 rtracklayer     CDS      92997     94054        NA           -        NA            NA          NA              NA
#> ...      ...         ...      ...       ...       ...       ...         ...       ...           ...         ...             ...
#>               Alias      genome          Parent
#>     <CharacterList> <character> <CharacterList>
#> 1    FLJ40100,TUBB8        hg19                
#> 2                            NA   GeneID:347688
#> 3                            NA   GeneID:347688
#> 4                            NA         872,873
#> 5                            NA         872,873
#> ...             ...         ...             ...
#>  [ reached 'max' / getOption("max.print") -- omitted 5 rows ]
```

This function returns a `DataFrame` object, which contains the same
information as the `GRanges` object but in a tabular format. Note that
one can use the `makeGRangesFromDataFrame` function to convert the
`DataFrame` returned by `readGFF` into a `GRanges` object, which is
often more convenient for downstream analysis:

``` r
readGFF(gff_file, version = "3") |>
    makeGRangesFromDataFrame(
        keep.extra.columns = TRUE
    )
#> GRanges object with 31 ranges and 10 metadata columns:
#>        seqnames      ranges strand |      source     type     score     phase            ID        Name        geneName           Alias
#>           <Rle>   <IRanges>  <Rle> |    <factor> <factor> <numeric> <integer>   <character> <character>     <character> <CharacterList>
#>    [1]    chr10 92828-95504      - | rtracklayer     gene         5      <NA> GeneID:347688       TUBB8 tubulin, beta 8  FLJ40100,TUBB8
#>    [2]    chr10 92828-95178      - | rtracklayer     mRNA        NA      <NA>           873       TUBB8            <NA>                
#>    [3]    chr10 92828-95504      - | rtracklayer     mRNA        NA      <NA>           872       TUBB8            <NA>                
#>    [4]    chr10 92828-94054      - | rtracklayer     exon        NA      <NA>          <NA>        <NA>            <NA>                
#>    [5]    chr10 92997-94054      - | rtracklayer     CDS         NA      <NA>          <NA>        <NA>            <NA>                
#>    ...      ...         ...    ... .         ...      ...       ...       ...           ...         ...             ...             ...
#>   [27]    chr12 89675-89827      + | rtracklayer     CDS         NA      <NA>          <NA>        <NA>            <NA>                
#>   [28]    chr12 90587-90655      + | rtracklayer     exon        NA      <NA>          <NA>        <NA>            <NA>                
#>   [29]    chr12 90587-90655      + | rtracklayer     CDS         NA      <NA>          <NA>        <NA>            <NA>                
#>   [30]    chr12 90796-91263      + | rtracklayer     exon        NA      <NA>          <NA>        <NA>            <NA>                
#>   [31]    chr12 90796-91263      * | rtracklayer     CDS         NA      <NA>          <NA>        <NA>            <NA>                
#>             genome          Parent
#>        <character> <CharacterList>
#>    [1]        hg19                
#>    [2]        <NA>   GeneID:347688
#>    [3]        <NA>   GeneID:347688
#>    [4]        <NA>         872,873
#>    [5]        <NA>         872,873
#>    ...         ...             ...
#>   [27]        <NA>            4644
#>   [28]        <NA>            4644
#>   [29]        <NA>            4644
#>   [30]        <NA>            4644
#>   [31]        <NA>            4644
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

### Reading remote GFF files

You can also read GFF files from remote URLs. The `import` function can
handle URLs directly, allowing you to read GFF files hosted on remote
servers. On the Bioconductor Build System (BBS), we use the
`BiocFileCache` package to cache remote files, which avoids downloading
the file every time the code gets run.

``` r
library(BiocFileCache)
bfc <- BiocFileCache::BiocFileCache()
url <- "https://www.mirbase.org/download/hsa.gff3"
bquery <- bfcquery(bfc, url, "rname", exact = TRUE)
if (!nrow(bquery))
    bfcadd(x = bfc, rname = url, rtype = "web", download = TRUE)
gff_remote <-
    bfcrpath(bfc, rnames = url, exact = TRUE, download = FALSE, rtype = "web")
import(gff_remote)
#> GRanges object with 4801 ranges and 8 metadata columns:
#>          seqnames            ranges strand |   source                     type     score     phase             ID           Alias
#>             <Rle>         <IRanges>  <Rle> | <factor>                 <factor> <numeric> <integer>    <character> <CharacterList>
#>      [1]     chr1       17369-17436      - |       NA miRNA_primary_transcript        NA      <NA>      MI0022705       MI0022705
#>      [2]     chr1       17409-17431      - |       NA miRNA                           NA      <NA>   MIMAT0027618    MIMAT0027618
#>      [3]     chr1       17369-17391      - |       NA miRNA                           NA      <NA>   MIMAT0027619    MIMAT0027619
#>      [4]     chr1       30366-30503      + |       NA miRNA_primary_transcript        NA      <NA>      MI0006363       MI0006363
#>      [5]     chr1       30438-30458      + |       NA miRNA                           NA      <NA>   MIMAT0005890    MIMAT0005890
#>      ...      ...               ...    ... .      ...                      ...       ...       ...            ...             ...
#>   [4797]     chrY   2609229-2609252      + |       NA miRNA                           NA      <NA> MIMAT0023714_1    MIMAT0023714
#>   [4798]     chrY   4606120-4606228      + |       NA miRNA_primary_transcript        NA      <NA>      MI0032313       MI0032313
#>   [4799]     chrY   4606140-4606158      + |       NA miRNA                           NA      <NA>   MIMAT0039763    MIMAT0039763
#>   [4800]     chrY 13479177-13479266      + |       NA miRNA_primary_transcript        NA      <NA>      MI0039722       MI0039722
#>   [4801]     chrY 13479189-13479214      + |       NA miRNA                           NA      <NA>   MIMAT0049014    MIMAT0049014
#>                     Name Derives_from
#>              <character>  <character>
#>      [1]  hsa-mir-6859-1         <NA>
#>      [2] hsa-miR-6859-5p    MI0022705
#>      [3] hsa-miR-6859-3p    MI0022705
#>      [4]  hsa-mir-1302-2         <NA>
#>      [5]    hsa-miR-1302    MI0006363
#>      ...             ...          ...
#>   [4797]    hsa-miR-6089    MI0023563
#>   [4798]    hsa-mir-9985         <NA>
#>   [4799]    hsa-miR-9985    MI0032313
#>   [4800]   hsa-mir-12120         <NA>
#>   [4801]   hsa-miR-12120    MI0039722
#>   -------
#>   seqinfo: 24 sequences from an unspecified genome; no seqlengths
```

Otherwise, you can read the remote GFF file directly using the `readGFF`
function:

``` r
import(url, version = "3")
```

## Conversion to GFF

While `TxDb` objects are highly efficient for querying transcript
annotations within R, it is often necessary to export this data to a
standard, portable format for use with external tools or for sharing
with collaborators. The GFF3 format is a widely accepted standard for
this purpose. Converting a `TxDb` or a derivative object (like a
`GRangesList` of exons) into a GFF3-compatible `GRanges` object allows
for easy export. This is particularly useful for visualizing annotations
in genome browsers like IGV or UCSC, or for input into downstream
analysis pipelines that expect GFF3 files.

``` r
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
exonsBy(txdb, by = "tx") |>
    asGFF()
#> GRanges object with 825453 ranges and 7 metadata columns:
#>                  seqnames        ranges strand |        type          ID        Name   exon_id   exon_name exon_rank      Parent
#>                     <Rle>     <IRanges>  <Rle> | <character> <character> <character> <integer> <character> <integer> <character>
#>        [1]           chr1   11874-14409      + |        mRNA       mRNA1           1      <NA>        <NA>      <NA>        <NA>
#>        [2]           chr1   11874-14409      + |        mRNA       mRNA2           2      <NA>        <NA>      <NA>        <NA>
#>        [3]           chr1   11874-14409      + |        mRNA       mRNA3           3      <NA>        <NA>      <NA>        <NA>
#>        [4]           chr1   69091-70008      + |        mRNA       mRNA4           4      <NA>        <NA>      <NA>        <NA>
#>        [5]           chr1 321084-321115      + |        mRNA       mRNA5           5      <NA>        <NA>      <NA>        <NA>
#>        ...            ...           ...    ... .         ...         ...         ...       ...         ...       ...         ...
#>   [825449] chrUn_gl000241   22732-22846      - |        exon  exon742489        <NA>    289961        <NA>         6   mRNA82957
#>   [825450] chrUn_gl000241   20433-20481      - |        exon  exon742490        <NA>    289960        <NA>         7   mRNA82957
#>   [825451] chrUn_gl000243   11501-11530      + |        exon  exon742491        <NA>    289967        <NA>         1   mRNA82958
#>   [825452] chrUn_gl000243   13608-13637      + |        exon  exon742492        <NA>    289968        <NA>         1   mRNA82959
#>   [825453] chrUn_gl000247     5787-5816      - |        exon  exon742493        <NA>    289969        <NA>         1   mRNA82960
#>   -------
#>   seqinfo: 93 sequences (1 circular) from hg19 genome
```

## GFF to TxDb

You can convert a `GFF` object to a `TxDb` object using the
`makeTxDbFromGFF` in the
*[txdbmaker](https://bioconductor.org/packages/3.22/txdbmaker)* package.
This is useful for creating a transcript database from GFF annotations,
which can then be used for various genomic analyses.

``` r
library(txdbmaker)
txdb <- makeTxDbFromGFF(
    file = gff_remote,
    format = "gff3",
    dataSource = "https://www.mirbase.org/download/hsa.gff3",
    organism = "Homo sapiens",
    taxonomyId = 9606
)
#> Import genomic features from the file as a GRanges object ... OK
#> Prepare the 'metadata' data frame ... OK
#> Make the TxDb object ...
#> Warning in .extract_transcripts_from_GRanges(tx_IDX, gr, mcols0$type, mcols0$ID, : the transcript names ("tx_name" column in the TxDb
#> object) imported from the "Name" attribute are not unique
#> Warning in .makeTxDb_normarg_chrominfo(chrominfo): genome version information is not available for this TxDb object
#> OK
genome <- grepv("genome-build-id", readLines(gff_remote)) |>
    strsplit("# genome-build-id:\\s+") |>
    unlist() |>
    tail(1L)
genome(txdb) <- genome
txdb
#> TxDb object:
#> # Db type: TxDb
#> # Supporting package: GenomicFeatures
#> # Data source: https://www.mirbase.org/download/hsa.gff3
#> # Organism: Homo sapiens
#> # Taxonomy ID: 9606
#> # Genome: NA
#> # Nb of transcripts: 4801
#> # Db created by: txdbmaker package from Bioconductor
#> # Creation time: 2025-09-11 17:01:07 -0400 (Thu, 11 Sep 2025)
#> # txdbmaker version at creation time: 1.5.6
#> # RSQLite version at creation time: 2.4.3
#> # DBSCHEMAVERSION: 1.2
```

# SessionInfo

``` r
sessionInfo()
#> R version 4.5.1 Patched (2025-06-14 r88325)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
#>  [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: America/New_York
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#>  [1] BSgenome.Hsapiens.UCSC.hg19_1.4.3       BSgenome_1.77.1                         rtracklayer_1.69.1                     
#>  [4] txdbmaker_1.5.6                         TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2 GenomicFeatures_1.61.6                 
#>  [7] AnnotationDbi_1.71.1                    Biobase_2.69.0                          BiocStyle_2.37.1                       
#> [10] BiocFileCache_2.99.6                    dbplyr_2.5.1                            Rsamtools_2.25.3                       
#> [13] Biostrings_2.77.2                       XVector_0.49.0                          BiocIO_1.19.0                          
#> [16] GenomicRanges_1.61.2                    Seqinfo_0.99.2                          IRanges_2.43.0                         
#> [19] S4Vectors_0.47.0                        BiocGenerics_0.55.1                     generics_0.1.4                         
#> [22] Bioc.gff_0.99.13                        colorout_1.3-2                         
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyselect_1.2.1            dplyr_1.1.4                 blob_1.2.4                  filelock_1.0.3             
#>  [5] bitops_1.0-9                fastmap_1.2.0               lazyeval_0.2.2              RCurl_1.98-1.17            
#>  [9] GenomicAlignments_1.45.2    promises_1.3.3              rex_1.2.1                   XML_3.99-0.19              
#> [13] digest_0.6.37               lifecycle_1.0.4             processx_3.8.6              KEGGREST_1.49.1            
#> [17] RSQLite_2.4.3               magrittr_2.0.3              compiler_4.5.1              progress_1.2.3             
#> [21] rlang_1.1.6                 tools_4.5.1                 yaml_2.3.10                 knitr_1.50                 
#> [25] prettyunits_1.2.0           S4Arrays_1.9.1              bit_4.6.0                   curl_7.0.0                 
#> [29] DelayedArray_0.35.2         xml2_1.4.0                  abind_1.4-8                 BiocParallel_1.43.4        
#> [33] rsconnect_1.5.1             websocket_1.4.4             covr_3.6.4                  withr_3.0.2                
#> [37] purrr_1.1.0                 grid_4.5.1                  biomaRt_2.65.7              SummarizedExperiment_1.39.1
#> [41] cli_3.6.5                   rmarkdown_2.29              crayon_1.5.3                tinytest_1.4.1             
#> [45] rstudioapi_0.17.1           httr_1.4.7                  rjson_0.2.23                BiocBaseUtils_1.11.2       
#> [49] DBI_1.2.3                   cachem_1.1.0                chromote_0.5.1              stringr_1.5.2              
#> [53] rvest_1.0.5                 parallel_4.5.1              selectr_0.4-2               BiocManager_1.30.26        
#> [57] restfulr_0.0.16             matrixStats_1.5.0           vctrs_0.6.5                 Matrix_1.7-4               
#> [61] jsonlite_2.0.0              hms_1.1.3                   bit64_4.6.0-1               glue_1.8.0                 
#> [65] codetools_0.2-20            ps_1.9.1                    stringi_1.8.7               later_1.4.4                
#> [69] GenomeInfoDb_1.45.10        UCSC.utils_1.5.0            tibble_3.3.0                pillar_1.11.0              
#> [73] rappdirs_0.3.3              htmltools_0.5.8.1           GenomeInfoDbData_1.2.14     R6_2.6.1                   
#> [77] httr2_1.2.1                 lattice_0.22-7              evaluate_1.0.5              png_0.1-8                  
#> [81] memoise_2.0.1               Rcpp_1.1.0                  SparseArray_1.9.1           xfun_0.53                  
#> [85] MatrixGenerics_1.21.0       pkgconfig_2.0.3
```
