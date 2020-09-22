library(biomaRt)
cache <- file.path(tempdir(), "biomart_cache_test")
Sys.setenv(BIOMART_CACHE = cache)

context("Result caching")

go <- c("GO:0051330","GO:0000080","GO:0000114","GO:0000082")
chrom <- c(17,20,"Y")
attributes <- "hgnc_symbol"
filters <- c("go","chromosome_name")
values <- list(go, chrom)
ensembl <- useEnsembl("ensembl", "hsapiens_gene_ensembl", mirror = "www")

test_that("Hashing is order insensitive", {
   expect_identical(
       biomaRt:::.createHash(ensembl, attributes, filters, values),
       biomaRt:::.createHash(ensembl, rev(attributes), rev(filters), rev(values))
   )
})

test_that("Environment variable for cache location is used", {
    expect_message(biomartCacheInfo(),
                   regexp = "biomart_cache_test")
})

test_that("We find cache for previous query", {
    
    mart <- useEnsembl(biomart = "ensembl", 
                       mirror = "useast", 
                       dataset ="mmusculus_gene_ensembl")
    
    expect_message(res <- getBM(filter = "ensembl_gene_id",
                 values = "ENSMUSG00000028798",
                 attributes = c("ensembl_transcript_id", "neugenii_homolog_canonical_transcript_protein"),
                 mart = mart,
                 verbose = TRUE),  
                 regexp = "Cache found")
})

test_that("Cache details are printed", {
    expect_message( biomartCacheInfo(),
                   regexp = "biomaRt cache")  
})

test_that("Cache can be cleared", {
    cache_file <- biomartCacheInfo()
    expect_true( file.exists( cache_file ) )
    expect_silent( biomaRt:::biomartCacheClear() )
    expect_false( file.exists( cache_file) )
})