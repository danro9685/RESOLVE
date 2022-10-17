Sys.setenv("R_TESTS" = "")

library("testthat")
library("RESOLVE")

test_check("RESOLVE")
