library(biomaRt)

########################
context('getBM()')
########################

test_that("Fail with no arguments", {
    expect_error(getBM(), "You must provide a valid Mart object")
})

test_that("Fail when no dataset is specified", {
     ensembl=useMart("ensembl")
     expect_error(getBM(mart = ensembl), "No dataset selected, please select a dataset first")
})

 
