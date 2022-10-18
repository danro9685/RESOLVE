context("RESOLVE")

data(background)
data(patients)
set.seed(12345)
res_denovo = signaturesDecomposition(x = patients[1:3,], 
                                     K = 3:4, 
                                     background_signature = background, 
                                     nmf_runs = 2, 
                                     sparsify = FALSE, 
                                     num_processes = 1)
test_that("RESOLVE can perform de novo signatures decomposition", {
    expect_equal(names(res_denovo),c("alpha","beta","measures"))
})

set.seed(12345)
res_assignment = signaturesAssignment(x = patients[1:3,], beta = res_denovo$beta[[1]], sparsify = FALSE)
test_that("RESOLVE can perform signatures assignment", {
    expect_equal(names(res_assignment),c("alpha","beta"))
})

set.seed(12345)
res_cv = signaturesCV(x = patients[1:3,], 
                      beta = res_denovo$beta, 
                      cross_validation_iterations = 2, 
                      cross_validation_repetitions = 2, 
                      num_processes = 1)
test_that("RESOLVE can perform optimanl number of signatures estimation", {
    expect_equal(names(res_cv),c("estimates","summary"))
})

set.seed(12345)
res_sig = signaturesSignificance(x = patients[1:3,], 
                                 beta = res_denovo$beta[[1]], 
                                 cosine_thr = 0.95, 
                                 min_contribution = 0.05, 
                                 pvalue_thr = 0.05, 
                                 sparsify = FALSE, 
                                 nboot = 2, 
                                 num_processes = 1)
test_that("RESOLVE can perform the signatures significance estimation", {
    expect_equal(names(res_sig),c("alpha","beta","goodness_fit","bootstrap_estimates"))
})
