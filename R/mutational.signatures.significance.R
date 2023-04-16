#' Perform a robust estimation of alpha coefficients by bootstrap to reach a certain level of cosine similarity 
#' given a set of observed counts x and discovered signatures beta.
#'
#' @examples
#' data(background)
#' data(patients)
#' set.seed(12345)
#' beta <- signaturesDecomposition(x = patients[seq_len(3),seq_len(2)],
#'                                 K = 3:4,
#'                                 background_signature = background[seq_len(2)],
#'                                 nmf_runs = 2,
#'                                 num_processes = 1)
#' set.seed(12345)
#' res <- signaturesSignificance(x = patients[seq_len(3),seq_len(2)],
#'                               beta = beta$beta[[1]],
#'                               cosine_thr = 0.95,
#'                               min_contribution = 0.05,
#'                               pvalue_thr = 0.05,
#'                               nboot = 5,
#'                               num_processes = 1)
#'
#' @title signaturesSignificance
#' @param x Counts matrix for a set of n patients and m categories. These can be, e.g., SBS, MNV, CN or CN counts;
#' in the case of SBS it would be an n patients x 96 trinucleotides matrix.
#' @param beta Discovered signatures to be used for the fit of alpha.
#' @param cosine_thr Level of cosine similarity to be reached for the fit of alpha.
#' @param min_contribution Minimum contribution of a signature to be considered significant.
#' @param pvalue_thr Pvalue level to be used to assess significance.
#' @param nboot Number of bootstrap iterations to be performed.
#' @param num_processes Number of processes to be used during parallel execution. To execute in single process mode,
#' this parameter needs to be set to either NA or NULL.
#' @param verbose Boolean. Shall I print information messages?
#' @return A list with the bootstrap estimates. It includes 5 elements:
#'              alpha: matrix of the discovered exposure values considering significant signatures as estimated by bootstrap.
#'              beta: matrix of the discovered signatures.
#'              unexplained_mutations: number of unexplained mutations per sample.
#'              goodness_fit: vector reporting cosine similarities between predictions and observations.
#'              bootstrap_estimates: list of matrices reporting results by bootstrap estimates.
#' @export signaturesSignificance
#'
signaturesSignificance <- function( x, beta, cosine_thr = 0.95, min_contribution = 0.05,
    pvalue_thr = 0.05, nboot = 100, num_processes = Inf, verbose = TRUE ) {

    # check input parameters
    x <- as.matrix(x)
    beta <- as.matrix(beta)
    if (nboot < 5) {
        warning("The minimum value of nboot must be 5...")
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

    # performing the inference
    if (num_processes == 1) { # performing the inference sequentially
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
                return(.signatures_significance(x = curr_patient, beta = beta, cosine_thr = cosine_thr))
            })
            alpha[[boot_iteration]] <- Reduce("rbind", curr_alpha)
        }
    } else { # performing the inference in parallel
        parallel <- makeCluster(num_processes, outfile = "")
        close_parallel <- TRUE
        res_clusterEvalQ <- clusterEvalQ(parallel, library("glmnet", warn.conflicts = FALSE,
            quietly = TRUE, verbose = FALSE))
        res_clusterEvalQ <- clusterEvalQ(parallel, library("lsa", warn.conflicts = FALSE,
            quietly = TRUE, verbose = FALSE))
        clusterExport(parallel, varlist = c(".signatures_significance", "verbose"), envir = environment())
        clusterExport(parallel, varlist = c("nboot", "x", "beta", "cosine_thr"), envir = environment())
        clusterSetRNGStream(parallel, iseed = round(runif(1) * 1e+05))
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
                return(.signatures_significance(x = curr_patient, beta = beta, cosine_thr = cosine_thr))
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
    pvalues[which(alpha == 0)] <- NA
    for (i in seq_len(nrow(pvalues))) {
        for (j in seq_len(ncol(pvalues))) {
            if (!is.na(pvalues[i, j])) {
                pvalues[i, j] <- NA
                curr_values <- NULL
                for (k in seq_len(length(alpha_distibution))) {
                    curr_values <- c(curr_values, ((alpha_distibution[[k]][i, j])/sum((alpha_distibution[[k]][i, ]))))
                }
                pvalues[i, j] <- wilcox.test(as.numeric(curr_values), alternative = "greater",
                    mu = min_contribution)$p.value
            }
        }
    }

    bootstrap <- list(estimate = alpha, pvalues = pvalues)

    if (verbose) {
        message("Performing fit of alpha considering only signatures with significant contribution...", "\n")
    }

    # perform final fit of alpha using only significant signatures
    alpha <- matrix(0, nrow = nrow(x), ncol = nrow(beta))
    rownames(alpha) <- rownames(x)
    colnames(alpha) <- rownames(beta)
    goodness_fit <- NULL
    unexplained_mutations <- rep(NA, nrow(x))
    names(unexplained_mutations) <- rownames(x)
    for (i in seq_len(nrow(x))) {
        # perform fit of alpha
        curr_sigs <- names(which(bootstrap$pvalues[i, ] < pvalue_thr))
        curr_alpha <- matrix(NA, nrow=1, ncol=length(curr_sigs))
        curr_beta <- beta[curr_sigs, , drop = FALSE]
        if (nrow(curr_beta) > 1) {
            curr_alpha[1, ] <- tryCatch({
                res <- cv.glmnet(x = t(curr_beta), y = as.vector(x[i, ]), 
                        type.measure = "mse", nfolds = 10, nlambda = 100, 
                        family = "gaussian", alpha = 1, lower.limits = 0, 
                        maxit = 1e+05)
                res <- as.numeric(coef(res,s=res$lambda.min))
                unexplained_mutations[i] <- res[1]
                res <- res[-1]
                res
            }, error = function( e ) {
                res <- (sum(x[i,])/ncol(curr_alpha))
                return(res)
            }, finally = {
                gc(verbose = FALSE)
            })
        } else {
            fit_inputs <- cbind(t(curr_beta), rep(0, ncol(curr_beta)))
            curr_alpha[1, ] <- tryCatch({
                res <- cv.glmnet(x = fit_inputs, y = as.vector(x[i, ]), 
                        type.measure = "mse", nfolds = 10, nlambda = 100, 
                        family = "gaussian", alpha = 1, lower.limits = 0, 
                        maxit = 1e+05)
                res <- as.numeric(coef(res,s=res$lambda.min))
                res <- res[-length(res)]
                unexplained_mutations[i] <- res[1]
                res <- res[-1]
                res
            }, error = function( e ) {
                res <- (sum(x[i,])/ncol(curr_alpha))
                return(res)
            }, finally = {
                gc(verbose = FALSE)
            })
        }
        alpha[i, rownames(curr_beta)] <- as.numeric(curr_alpha)
        # estimate goodness of fit
        pred_counts <- ((alpha[i, ] %*% beta) + unexplained_mutations[i])
        goodness_fit <- c(goodness_fit, as.numeric(cosine(as.numeric(x[i, ]), as.numeric(pred_counts))))
    }
    names(goodness_fit) <- rownames(x)

    # close parallel
    if (close_parallel) {
        stopCluster(parallel)
    }
    rm(close_parallel)
    rm(parallel)
    gc(verbose = FALSE)

    # save the results
    results <- list(alpha = alpha, beta = beta, unexplained_mutations = unexplained_mutations, goodness_fit = goodness_fit, bootstrap_estimates = bootstrap)

    # return the estimeted signatures contributions
    return(results)

}
