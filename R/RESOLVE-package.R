#' An R package for the efficient analysis of mutational signatures from cancer genomes.
#'
#' This package provides a framework that allows for the efficient extraction and assignment of mutational signatures.
#' It implements a novel algorithm that enables: (i) the efficient extraction, (ii) exposure estimation, and
#' (iii) confidence assessment during the computational inference of mutational signatures.
#' RESOLVE performs de novo signatures extraction by regularized Non-Negative Matrix Factorization. The method incorporates
#' a background signature during the inference step and adopts elastic net regression with LASSO penalty to reduce the 
#' impact of overfitting. The estimation of the optimal number of signatures is performed by bi-cross-validation.
#'
#' @name RESOLVE-package
#' @aliases RESOLVE
#' @docType package
#' @keywords package
#' @import parallel
#' @import IRanges
#' @import GenomeInfoDb
#' @import BSgenome.Hsapiens.1000genomes.hs37d5
#' @importFrom glmnet cv.glmnet
#' @importFrom lsa cosine
#' @importFrom stats coef runif wilcox.test
#' @importFrom data.table data.table dcast .N
#' @importFrom Biostrings DNAStringSet complement reverseComplement subseq
#' @importFrom BSgenome getSeq
#' @importFrom GenomicRanges GRanges seqnames
#' @importFrom methods is
#'
NULL
