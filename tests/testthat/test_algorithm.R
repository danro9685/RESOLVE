context("SparseSignaturesPlus")

data("ssm560_reduced")
test_that("SparseSignaturesPlus input is correct", {
    expect_equal(colnames(ssm560_reduced),c("sample","chrom","start","end","ref","alt"))
})
