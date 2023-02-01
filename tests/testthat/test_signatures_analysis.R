context("RESOLVE")

data(background)
data(patients)
set.seed(12345)
res_denovo = signaturesDecomposition(x = patients[seq_len(3),seq_len(2)], 
                                     K = 3:4, 
                                     background_signature = background[seq_len(2)], 
                                     nmf_runs = 2, 
                                     num_processes = 1)
test_that("RESOLVE can perform de novo signatures decomposition", {
    expect_equal(names(res_denovo),c("alpha","beta","unexplained_mutations","cosine_similarity","measures"))
})

set.seed(12345)
res_assignment = signaturesAssignment(x = patients[seq_len(3),seq_len(2)], beta = res_denovo$beta[[1]])
test_that("RESOLVE can perform signatures assignment", {
    expect_equal(names(res_assignment),c("alpha","beta","unexplained_mutations"))
})

set.seed(12345)
res_cv = signaturesCV(x = patients[seq_len(3),seq_len(2)], 
                      beta = res_denovo$beta, 
                      cross_validation_iterations = 2, 
                      cross_validation_repetitions = 2, 
                      num_processes = 1)
test_that("RESOLVE can perform optimanl number of signatures estimation", {
    expect_equal(names(res_cv),c("estimates","summary"))
})

set.seed(12345)
res_sig = signaturesSignificance(x = patients[seq_len(3),seq_len(2)], 
                                 beta = res_denovo$beta[[1]], 
                                 cosine_thr = 0.95, 
                                 min_contribution = 0.05, 
                                 pvalue_thr = 0.05, 
                                 nboot = 5, 
                                 num_processes = 1)
test_that("RESOLVE can perform the signatures significance estimation", {
    expect_equal(names(res_sig),c("alpha","beta","unexplained_mutations","goodness_fit","bootstrap_estimates"))
})
