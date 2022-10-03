Sys.setenv("R_TESTS" = "")

library("testthat")
library("SparseSignaturesPlus")

test_check("SparseSignaturesPlus")
