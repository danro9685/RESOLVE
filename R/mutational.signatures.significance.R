#' Perform a robust estimation by bootstrap of alpha coefficients until a given level of cosine similarity is reached
#' given a set of observed counts x and discovered signatures beta.
#'
#' @examples
#' data(background)
#' data(patients)
#' set.seed(12345)
#' beta <- signaturesDecomposition(x = patients[seq_len(3),],
#'                                K = 3:4,
#'                                background_signature = background,
#'                                nmf_runs = 2,
#'                                sparsify = FALSE,
#'                                num_processes = 1)
#' set.seed(12345)
#' res <- signaturesSignificance(x = patients[seq_len(3),],
#'                              beta = beta$beta[[1]],
#'                              cosine_thr = 0.95,
#'                              min_contribution = 0.05,
#'                              pvalue_thr = 0.05,
#'                              sparsify = FALSE,
#'                              nboot = 2,
#'                              num_processes = 1)
#'
#' @title signaturesSignificance
#' @param x counts matrix for a set of n patients and m categories. These can be, e.g., trinucleotides counts for n patients and 96 trinucleotides.
#' @param beta discovered signatures to be used for the fit of alpha.
#' @param cosine_thr Level of cosine similarity to be reached for the fit of alpha.
#' @param min_contribution Minimum contribution of a signature to be considered significant.
#' @param pvalue_thr Pvalue level to be used to assess significance.
#' @param sparsify boolean; Shall I perform regularization using LASSO?
#' @param nboot Number of bootstrap iterations to be performed.
#' @param num_processes Number of processes to be used during parallel execution. To execute in single process mode,
#' this parameter needs to be set to either NA or NULL.
#' @param verbose boolean; Shall I print information messages?
#' @return A list with the bootstrap estimates. It includes 4 elements:
#'              alpha: matrix of the discovered exposure values considering significant signatures as estimated by bootstrap.
#'              beta: matrix of the discovered signatures.
#'              goodness_fit: vector reporting cosine similarities between predictions and observations.
#'              bootstrap_estimates: list of matrices reporting results by bootstrap estimates.
#' @export signaturesSignificance
#' @import nnls
#' @import parallel
#' @importFrom glmnet cv.glmnet
#' @importFrom lsa cosine
#'
signaturesSignificance <- function(x, beta, cosine_thr = 0.95, min_contribution = 0.05,
    pvalue_thr = 0.05, sparsify = TRUE, nboot = 100, num_processes = Inf, verbose = TRUE) {
    # check input parameters
    x <- as.matrix(x)
    beta <- as.matrix(beta)
    if (nboot < 5) {
        warning("The minimum value of nboot can be 5...")
        nboot <- 5
    }

    # setting up parallel execution
    parallel <- NULL
    close_parallel <- FALSE
    if (is.na(num_processes) || is.null(num_processes)) {
        num_processes <- 1
    } else if (num_processes == Inf) {
        cores <- as.integer((detectCores() - 1))
        num_processes <- min(cores, nboot)
    } else {
        num_processes <- min(num_processes, nboot)
    }

    if (verbose) {
        message("Estimating the contribution of each signature to the fit with a total of ",
            nboot, " bootstrap iterations...", "\n")
        if (num_processes > 1) {
            message("Executing ", num_processes, " processes in parallel...", "\n")
        }
    }

    # performing inference
    if (num_processes == 1) {
        # performing inference sequentially
        alpha <- list()
        for (boot_iteration in seq_len(nboot)) {
            if (verbose) {
                message("Performing iteration ", boot_iteration, " out of ",
                  nboot, "...", "\n")
            }

            curr_alpha <- lapply(X = seq_len(nrow(x)), FUN = function(counts) {
                curr_patient <- x[counts, , drop = FALSE]
                curr_num_counts <- sum(curr_patient)
                curr_counts_distribution <- as.numeric(curr_patient/rowSums(curr_patient))
                curr_boot_sampling <- table(sample(x = colnames(curr_patient), size = curr_num_counts,
                  replace = TRUE, prob = curr_counts_distribution))
                curr_patient[1, ] <- 0
                curr_patient[1, names(curr_boot_sampling)] <- as.numeric(curr_boot_sampling)
                return(.sigs_significance(x = curr_patient, beta = beta, cosine_thr = cosine_thr,
                  sparsify = sparsify))
            })
            alpha[[boot_iteration]] <- Reduce("rbind", curr_alpha)

        }

    } else {
        # performing inference in parallel
        parallel <- makeCluster(num_processes, outfile = "")
        close_parallel <- TRUE
        res_clusterEvalQ <- clusterEvalQ(parallel, library("nnls", warn.conflicts = FALSE,
            quietly = TRUE, verbose = FALSE))
        res_clusterEvalQ <- clusterEvalQ(parallel, library("glmnet", warn.conflicts = FALSE,
            quietly = TRUE, verbose = FALSE))
        res_clusterEvalQ <- clusterEvalQ(parallel, library("lsa", warn.conflicts = FALSE,
            quietly = TRUE, verbose = FALSE))
        clusterExport(parallel, varlist = c(".sigs_significance"), envir = environment())
        clusterExport(parallel, varlist = c("verbose", "nboot"), envir = environment())
        clusterExport(parallel, varlist = c("x", "beta", "cosine_thr", "sparsify"),
            envir = environment())
        clusterSetRNGStream(parallel, iseed = round(runif(1) * 1e+05))
        rm(res_clusterEvalQ)
        gc(verbose = FALSE)
        alpha <- parLapply(parallel, seq_len(nboot), function(boot_iteration) {
            if (verbose) {
                message("Performing iteration ", boot_iteration, " out of ",
                  nboot, "...", "\n")
            }

            curr_alpha <- lapply(X = seq_len(nrow(x)), FUN = function(counts) {
                curr_patient <- x[counts, , drop = FALSE]
                curr_num_counts <- sum(curr_patient)
                curr_counts_distribution <- as.numeric(curr_patient/rowSums(curr_patient))
                curr_boot_sampling <- table(sample(x = colnames(curr_patient), size = curr_num_counts,
                  replace = TRUE, prob = curr_counts_distribution))
                curr_patient[1, ] <- 0
                curr_patient[1, names(curr_boot_sampling)] <- as.numeric(curr_boot_sampling)
                return(.sigs_significance(x = curr_patient, beta = beta, cosine_thr = cosine_thr,
                  sparsify = sparsify))
            })
            curr_alpha <- Reduce("rbind", curr_alpha)

            return(curr_alpha)

        })

    }

    # process results
    results <- alpha[[1]]
    if (length(alpha) > 1) {
        for (i in 2:length(alpha)) {
            results <- results + alpha[[i]]
        }
        results <- results/length(alpha)
    }
    alpha_distibution <- alpha
    alpha <- results

    if (verbose) {
        message("Estimating level of significance for each signature...", "\n")
    }
    pvalues <- alpha
    pvalues[which(pvalues == 0)] <- NA
    for (i in seq_len(nrow(pvalues))) {
        for (j in seq_len(ncol(pvalues))) {
            if (!is.na(pvalues[i, j])) {
                pvalues[i, j] <- NA
                curr_values <- NULL
                for (k in seq_len(length(alpha_distibution))) {
                  curr_values <- c(curr_values, (alpha_distibution[[k]][i, ]/sum(alpha_distibution[[k]][i,
                    ]))[j])
                }
                pvalues[i, j] <- wilcox.test(as.numeric(curr_values), alternative = "greater",
                  mu = min_contribution)$p.value
            }
        }
    }
    goodness_fit <- NULL
    for (i in seq_len(nrow(x))) {
        # estimate goodness of fit
        goodness_fit <- c(goodness_fit, as.numeric(cosine(as.numeric(x[i, ]), as.numeric((alpha[i,
            ] %*% beta)))))
    }
    names(goodness_fit) <- rownames(x)
    bootstrap <- list(estimate = alpha, pvalues = pvalues, goodness_fit = goodness_fit)

    if (verbose) {
        message("Performing fit of alpha considering only signatures with significant contribution...",
            "\n")
    }

    # perform final fit of alpha using only significant signatures
    alpha <- array(0, c(nrow(x), nrow(beta)))
    rownames(alpha) <- rownames(x)
    colnames(alpha) <- rownames(beta)
    goodness_fit <- NULL
    for (i in seq_len(nrow(x))) {
        # perform fit
        curr_sigs <- names(which(bootstrap$pvalues[i, ] < pvalue_thr))
        curr_alpha <- array(NA, c(1, length(curr_sigs)))
        curr_beta <- beta[curr_sigs, , drop = FALSE]
        if (sparsify == TRUE && nrow(curr_beta) > 1) {
            reg_out <- tryCatch({
                res <- cv.glmnet(t(curr_beta), as.vector(x[i, ]), type.measure = "mse",
                    nfolds = 10, nlambda = 10, family = "gaussian", lower.limits = 0)
                curr_alpha[1, ] <- as.numeric(res$glmnet.fit$beta[, ncol(res$glmnet.fit$beta)])
            }, error = function( e ) {
                curr_alpha[1, ] <- 0
            }, finally = {
                gc(verbose = FALSE)
            })
        } else {
            curr_alpha[1, ] <- nnls(t(curr_beta), as.vector(x[i, ]))$x
        }
        alpha[i, rownames(curr_beta)] <- as.numeric(curr_alpha)
        # estimate goodness of fit
        goodness_fit <- c(goodness_fit, as.numeric(cosine(as.numeric(x[i, ]), as.numeric((alpha[i,
            ] %*% beta)))))
    }
    names(goodness_fit) <- rownames(x)

    # close parallel
    if (close_parallel) {
        stopCluster(parallel)
    }
    rm(close_parallel)
    rm(parallel)
    gc(verbose = FALSE)

    # save results
    results <- list(alpha = alpha, beta = beta, goodness_fit = goodness_fit, bootstrap_estimates = bootstrap)

    # return the estimeted signatures contributions
    return(results)

}

# iteratively estimate alpha coefficients until a given level of cosine
# similarity is reached
.sigs_significance <- function(x, beta, cosine_thr = 0.95,
    sparsify = TRUE) {
    # initialization
    alpha <- array(NA, c(1, nrow(beta)))
    rownames(alpha) <- rownames(x)
    colnames(alpha) <- rownames(beta)
    if (sparsify == TRUE && nrow(beta) > 1) {
        reg_out <- tryCatch({
            res <- cv.glmnet(t(beta), as.vector(x[1, ]),
                type.measure = "mse", nfolds = 10, nlambda = 10,
                family = "gaussian", lower.limits = 0)
            alpha[1, ] <- as.numeric(res$glmnet.fit$beta[,
                ncol(res$glmnet.fit$beta)])
        }, error = function( e ) {
            alpha[1, ] <- 0
        }, finally = {
            gc(verbose = FALSE)
        })
    } else {
        alpha[1, ] <- nnls(t(beta), as.vector(x[1, ]))$x
    }

    # iteratively include signatures into the fit until a given level of cosine
    # similarity is reached
    sigs <- colnames(alpha[, which(alpha[1, ] > 0),
        drop = FALSE])[sort.int(alpha[, which(alpha[1,
        ] > 0)], decreasing = TRUE, index.return = TRUE)$ix]
    alpha <- array(0, c(1, nrow(beta)))
    rownames(alpha) <- rownames(x)
    colnames(alpha) <- rownames(beta)
    for (i in seq_len(length(sigs))) {
        # consider current set of signatures and perform fit of alpha
        curr_alpha <- array(NA, c(1, i))
        curr_beta <- beta[sigs[seq_len(i)], , drop = FALSE]
        if (sparsify == TRUE && nrow(curr_beta) > 1) {
            reg_out <- tryCatch({
                res <- cv.glmnet(t(curr_beta), as.vector(x[1,
                    ]), type.measure = "mse", nfolds = 10,
                    nlambda = 10, family = "gaussian", lower.limits = 0)
                curr_alpha[1, ] <- as.numeric(res$glmnet.fit$beta[,
                    ncol(res$glmnet.fit$beta)])
            }, error = function( e ) {
                curr_alpha[1, ] <- 0
            }, finally = {
                gc(verbose = FALSE)
            })
        } else {
            curr_alpha[1, ] <- nnls(t(curr_beta), as.vector(x[1,
                ]))$x
        }
        # estimate goodness of fit
        curr_predicted_counts <- curr_alpha %*% curr_beta
        curr_goodness_fit <- as.numeric(cosine(as.numeric(x),
            as.numeric(curr_predicted_counts)))
        if (curr_goodness_fit > cosine_thr) {
            break

        }
    }
    alpha[, rownames(curr_beta)] <- curr_alpha

    # return the estimated alpha coefficients
    return(alpha)

}
