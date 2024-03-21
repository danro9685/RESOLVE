context("RESOLVE")

library("data.table")
library("BSgenome.Hsapiens.1000genomes.hs37d5")
data(ssm560_reduced)
res_sbs = getSBSCounts(data = ssm560_reduced, reference = BSgenome.Hsapiens.1000genomes.hs37d5)
test_that("RESOLVE can import data for SBS signatures analysis", {
    expect_equal(ncol(res_sbs),96)
})

res_mnv = getMNVCounts(data = ssm560_reduced)
test_that("RESOLVE can import data for MNV signatures analysis", {
    expect_equal(ncol(res_mnv),78)
})
