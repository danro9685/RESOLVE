#' Perform the assignment of K somatic mutational signatures provided as input to samples given a set of observed counts x.
#' This function can be used to estimate different types of mutational signatures such as: SBS (single base substitutions) and 
#' MNV (multi-nucleotide variant) (see Degasperi, Andrea, et al. 'Substitution mutational signatures in whole-genome–sequenced cancers 
#' in the UK population.' Science 376.6591 (2022): abl9283), CX (chromosomal instability) (see Drews, Ruben M., et al. 'A pan-cancer 
#' compendium of chromosomal instability.' Nature 606.7916 (2022): 976-983) and CN (copy number) signatures (see Steele, 
#' Christopher D., et al. 'Signatures of copy number alterations in human cancer.' Nature 606.7916 (2022): 984-991).
#'
#' @examples
#' data(background)
#' data(patients)
#' set.seed(12345)
#' beta <- signaturesDecomposition(x = patients[seq_len(10),],
#'                                 K = 3,
#'                                 background_signature = background,
#'                                 nmf_runs = 2,
#'                                 sparsify = FALSE,
#'                                 num_processes = 1)
#' set.seed(12345)
#' res <- signaturesAssignment(x = patients[seq_len(10),], beta = beta$beta[[1]], sparsify = FALSE)
#'
#' @title signaturesAssignment
#' @param x Counts matrix for a set of n patients and m categories. These can be, e.g., SBS, MNV, CN or CN counts;
#' in the case of SBS it would be an n patients x 96 trinucleotides matrix.
#' @param beta Matrix of the discovered signatures to be used for the assignment.
#' @param normalize_counts If true, the input counts matrix x is normalized such that the patients have the same number of mutation.
#' @param sparsify Boolean. Shall I perform regularization?
#' @param verbose Boolean. Shall I print information messages?
#' @return A list with the discovered signatures. It includes 2 elements:
#'              alpha: matrix of the discovered exposure values.
#'              beta: matrix of the discovered signatures.
#' @export signaturesAssignment
#' @import glmnet
#'
signaturesAssignment <- function( x, beta, normalize_counts = FALSE,
    sparsify = TRUE, verbose = TRUE ) {

    # check input parameters
    x <- as.matrix(x)
    if (normalize_counts) {
        x_not_normalized <- x
        x <- (x/rowSums(x)) * 2500
    }

    if (verbose) {
        message("Performing signatures assignment...", "\n")
    }

    # initialize alpha with an empty matrix
    alpha <- matrix(NA, nrow=nrow(x), ncol=nrow(beta))
    rownames(alpha) <- rownames(x)
    colnames(alpha) <- rownames(beta)

    # perform signatures assignment
    if (sparsify) {
        for (j in seq_len(nrow(alpha))) {
            if (nrow(beta) > 1) {
                alpha[j, ] <- tryCatch({
                    res <- cv.glmnet(x = t(beta), y = as.vector(x[j, ]), 
                        type.measure = "mse", nfolds = 10, nlambda = 100, 
                        family = "gaussian", lower.limits = 0)
                    res <- as.numeric(coef(res,s=res$lambda.min))
                    res <- as.numeric(res[1]+res[-1])
                    if(any(res<0)) {
                        res[which(res<0)] <- 0
                    }
                    res
                }, error = function( e ) {
                    res <- (sum(as.vector(x[j,]))/ncol(alpha))
                    return(res)
                }, finally = {
                    gc(verbose = FALSE)
                })
            } else {
                res_inputs <- cbind(t(beta), rep(0, ncol(beta)))
                alpha[j, ] <- tryCatch({
                    res <- cv.glmnet(x = res_inputs, y = as.vector(x[j, ]), 
                        type.measure = "mse", nfolds = 10, nlambda = 100, 
                        family = "gaussian", lower.limits = 0)
                    res <- as.numeric(coef(res,s=res$lambda.min))
                    res <- as.numeric(res[1]+res[-1])[-length(res)]
                    if(any(res<0)) {
                        res[which(res<0)] <- 0
                    }
                    res
                }, error = function( e ) {
                    res <- (sum(as.vector(x[j,]))/ncol(alpha))
                    return(res)
                }, finally = {
                    gc(verbose = FALSE)
                })
            }
            gc(verbose = FALSE)
        }
    } else {
        for (j in seq_len(nrow(alpha))) {
            if (nrow(beta) > 1) {
                alpha[j, ] <- tryCatch({
                    res <- glmnet(x = t(beta), y = as.vector(x[j, ]), 
                        lambda = 0, family = "gaussian", lower.limits = 0)
                    res <- as.numeric(coef(res,s=res$lambda))
                    res <- as.numeric(res[1]+res[-1])
                    if(any(res<0)) {
                        res[which(res<0)] <- 0
                    }
                    res
                }, error = function( e ) {
                    res <- (sum(as.vector(x[j,]))/ncol(alpha))
                    return(res)
                }, finally = {
                    gc(verbose = FALSE)
                })
            } else {
                res_inputs <- cbind(t(beta), rep(0, ncol(beta)))
                alpha[j, ] <- tryCatch({
                    res <- glmnet(x = res_inputs, y = as.vector(x[j, ]), 
                        lambda = 0, family = "gaussian", lower.limits = 0)
                    res <- as.numeric(coef(res,s=res$lambda))
                    res <- as.numeric(res[1]+res[-1])[-length(res)]
                    if(any(res<0)) {
                        res[which(res<0)] <- 0
                    }
                    res
                }, error = function( e ) {
                    res <- (sum(as.vector(x[j,]))/ncol(alpha))
                    return(res)
                }, finally = {
                    gc(verbose = FALSE)
                })
            }
            gc(verbose = FALSE)
        }
    }

    # rescale alpha to the original magnitude
    if (normalize_counts) {
        alpha <- alpha * (rowSums(x_not_normalized)/2500)
    }

    # save results
    results <- list(alpha = alpha, beta = beta)
    rm(alpha)
    rm(beta)
    gc(verbose = FALSE)

    # return the assigned signatures
    return(results)

}

#' Perform signatures discovery and rank estimation for a range of K somatic mutational signatures given a set of observed counts x.
#' This function can be used to estimate different types of mutational signatures such as: SBS (single base substitutions) and 
#' MNV (multi-nucleotide variant) (see Degasperi, Andrea, et al. 'Substitution mutational signatures in whole-genome–sequenced cancers 
#' in the UK population.' Science 376.6591 (2022): abl9283), CX (chromosomal instability) (see Drews, Ruben M., et al. 'A pan-cancer 
#' compendium of chromosomal instability.' Nature 606.7916 (2022): 976-983) and CN (copy number) signatures (see Steele, 
#' Christopher D., et al. 'Signatures of copy number alterations in human cancer.' Nature 606.7916 (2022): 984-991).
#'
#' @examples
#' data(background)
#' data(patients)
#' set.seed(12345)
#' res <- signaturesDecomposition(x = patients[seq_len(10),],
#'                                K = 3:4,
#'                                background_signature = background,
#'                                nmf_runs = 2,
#'                                sparsify = FALSE,
#'                                num_processes = 1)
#'
#' @title signaturesDecomposition
#' @param x Counts matrix for a set of n patients and m categories. These can be, e.g., SBS, MNV, CN or CN counts;
#' in the case of SBS it would be an n patients x 96 trinucleotides matrix.
#' @param K Either one value or a range of numeric values (each of them greater than 0) indicating the number of signatures 
#' to be considered.
#' @param background_signature Background signature to be used.
#' @param normalize_counts If true, the input counts matrix x is normalized such that the patients have the same number of mutation.
#' @param nmf_runs Number of iteration (minimum 1) of NMF to be performed for a robust estimation of beta.
#' @param sparsify Boolean. Shall I perform regularization?
#' @param num_processes Number of processes to be used during parallel execution. To execute in single process mode,
#' this parameter needs to be set to either NA or NULL.
#' @param verbose Boolean. Shall I print information messages?
#' @return A list with the discovered signatures and related rank measures. It includes 3 elements:
#'              alpha: list of matrices of the discovered exposure values for each possible rank in the range K.
#'              beta: list of matrices of the discovered signatures for each possible rank in the range K.
#'              measures: a data.frame containing the quality measures for each possible rank in the range K.
#' @export signaturesDecomposition
#' @import glmnet
#' @import NMF
#' @import parallel
#' @importFrom cluster silhouette
#' @importFrom lsa cosine
#'
signaturesDecomposition <- function( x, K, background_signature = NULL,
    normalize_counts = FALSE, nmf_runs = 100, sparsify = TRUE,
    num_processes = Inf, verbose = TRUE ) {

    # check input parameters
    x <- as.matrix(x)
    K <- sort(unique(K))
    if (length(K) == 0) {
        K <- 1
    }
    if (min(K) < 1) {
        warning("The minimum value of K must be 1, removing invalid values...")
        K <- K[K >= 1]
        if (length(K) == 0) {
            K <- 1
            warning("No valid value of K is provided: setting K equals to 1...")
        }
    }
    if (normalize_counts) {
        x_not_normalized <- x
        x <- (x/rowSums(x)) * 2500
    }
    if (nmf_runs < 1) {
        warning("The minimum value of nmf_runs must be 1...")
        nmf_runs <- 1
    }

    # setting up parallel execution
    pbackend <- NULL
    close_parallel <- FALSE
    if (is.na(num_processes) || is.null(num_processes)) {
        num_processes <- 1
    } else if (num_processes == Inf) {
        cores <- as.integer((detectCores() - 1))
        num_processes <- min(cores, nmf_runs)
    } else {
        num_processes <- min(num_processes, nmf_runs)
    }
    if (num_processes == 1) {
        pbackend <- "seq"
    } else {
        pbackend <- makeCluster(num_processes)
        res_clusterEvalQ <- clusterEvalQ(pbackend, library("lsa",
            warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
        res_clusterEvalQ <- clusterEvalQ(pbackend, library("NMF",
            warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
        res_clusterEvalQ <- clusterEvalQ(pbackend, library("glmnet",
            warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
        clusterExport(pbackend, varlist = c(".nmf_seed",
            ".nmf_ls", ".nmf_reg"), envir = environment())
        clusterSetRNGStream(pbackend, iseed = round(runif(1) *
            1e+05))
        close_parallel <- TRUE
        rm(res_clusterEvalQ)
        gc(verbose = FALSE)
    }

    if (verbose) {
        message("Performing signatures discovery and rank estimation...", "\n")
        if (num_processes > 1) {
            message("Executing ", num_processes, " processes in parallel...", "\n")
        }
    }

    # handle special case when solution with only background_signature is also considered
    rank0 <- FALSE
    if (K[1] == 1 && !is.null(background_signature)) {
        if (verbose) {
            message("Performing inference for K=1...", "\n")
        }
        rank0 <- TRUE
        K <- K[-1]
        rank0_K <- 1
        rank0_beta <- matrix(background_signature, nrow = 1)
        rownames(rank0_beta) <- "Background"
        colnames(rank0_beta) <- colnames(x)
        rank0_alpha <- signaturesAssignment(x = x, beta = rank0_beta,
            normalize_counts = FALSE, sparsify = sparsify,
            verbose = FALSE)$alpha
        pred_rank0 <- rank0_alpha %*% rank0_beta
        rank0_rss <- rss(pred_rank0, x)
        rank0_evar <- evar(pred_rank0, x)
        rank0_sparseness_alpha <- sparseness(as.vector(rank0_alpha))
        rank0_sparseness_beta <- sparseness(as.vector(rank0_beta))
        if (normalize_counts) {
            rank0_alpha <- rank0_alpha * (rowSums(x_not_normalized)/2500)
        }
        if (nmf_runs > 1) {
            rank0_measures <- matrix(c(rank0_K, rank0_rss,
                rank0_evar, NA, NA, rank0_sparseness_alpha,
                rank0_sparseness_beta, 1, 1, 1), nrow = 1)
            colnames(rank0_measures) <- c("Rank", "Residual_Sum_Squares",
                "Explained_Variance", "Silhouette_Alpha",
                "Silhouette_Beta", "Sparseness_Alpha", "Sparseness_Beta",
                "Cophenetic_Coefficient", "Dispersion_Coefficient",
                "Silhouette_Consensus")
        } else {
            rank0_measures <- matrix(c(rank0_K, rank0_rss,
                rank0_evar, NA, NA, rank0_sparseness_alpha,
                rank0_sparseness_beta), nrow = 1)
            colnames(rank0_measures) <- c("Rank", "Residual_Sum_Squares",
                "Explained_Variance", "Silhouette_Alpha",
                "Silhouette_Beta", "Sparseness_Alpha", "Sparseness_Beta")
        }
    }

    # perform signatures discovery and rank estimation
    if (length(K) > 0) {

        alpha <- list()
        beta <- list()
        measures <- NULL
        if (rank0) {
            alpha[["1_signatures"]] <- rank0_alpha
            beta[["1_signatures"]] <- rank0_beta
            measures <- rank0_measures
        }

        for (i in seq_len(length(K))) {

            if (verbose) {
                message("Performing inference for K=", K[i], "...", "\n")
            }

            # perform the inference for the current K
            results <- NULL
            gc(verbose = FALSE)
            if (sparsify) {
                results <- nmf(x = x, rank = K[i], method = .nmf_reg, 
                        background = background_signature, objective = .nmf_objective, 
                        mixed = FALSE, seed = .nmf_seed, rng = round(runif(1) * 1e+05), 
                        nrun = nmf_runs, .pbackend = pbackend)
            } else {
                results <- nmf(x = x, rank = K[i], method = .nmf_ls, 
                        background = background_signature, objective = .nmf_objective, 
                        mixed = FALSE, seed = .nmf_seed, rng = round(runif(1) * 1e+05), 
                        nrun = nmf_runs, .pbackend = pbackend)
            }
            gc(verbose = FALSE)

            alpha[[paste0(K[i], "_signatures")]] <- basis(results)
            beta[[paste0(K[i], "_signatures")]] <- coef(results)

            # rescale alpha to the original magnitude
            if (normalize_counts) {
                alpha[[paste0(K[i], "_signatures")]] <- alpha[[paste0(K[i], "_signatures")]] * (rowSums(x_not_normalized)/2500)
            }

            # compute and save quality measures
            curr_rss <- rss(results, x)
            curr_evar <- evar(results, x)
            curr_silhouette_alpha <- tryCatch({
                mean(silhouette(results, what = "features")[, "sil_width"])
            }, error = function(e) {
                NA
            })
            curr_silhouette_beta <- tryCatch({
                mean(silhouette(results, what = "samples")[, "sil_width"])
            }, error = function(e) {
                NA
            })
            curr_sparseness_alpha <- sparseness(as.vector(alpha[[paste0(K[i], "_signatures")]]))
            curr_sparseness_beta <- sparseness(as.vector(beta[[paste0(K[i], "_signatures")]]))
            if (nmf_runs > 1) {
                curr_silhouette_consensus <- tryCatch({
                  mean(silhouette(results, what = "chc")[, "sil_width"])
                }, error = function(e) {
                  NA
                })
                curr_measures <- matrix(c(K[i], curr_rss,
                  curr_evar, curr_silhouette_alpha, curr_silhouette_beta,
                  curr_sparseness_alpha, curr_sparseness_beta,
                  cophcor(results), dispersion(results),
                  curr_silhouette_consensus), nrow = 1)
                colnames(curr_measures) <- c("Rank", "Residual_Sum_Squares",
                  "Explained_Variance", "Silhouette_Alpha",
                  "Silhouette_Beta", "Sparseness_Alpha",
                  "Sparseness_Beta", "Cophenetic_Coefficient",
                  "Dispersion_Coefficient", "Silhouette_Consensus")
            } else {
                curr_measures <- matrix(c(K[i], curr_rss,
                  curr_evar, curr_silhouette_alpha, curr_silhouette_beta,
                  curr_sparseness_alpha, curr_sparseness_beta),
                  nrow = 1)
                colnames(curr_measures) <- c("Rank", "Residual_Sum_Squares",
                  "Explained_Variance", "Silhouette_Alpha",
                  "Silhouette_Beta", "Sparseness_Alpha",
                  "Sparseness_Beta")
            }
            measures <- rbind(measures, curr_measures)
            rm(results)
            rm(curr_measures)
            gc(verbose = FALSE)

        }
        rownames(measures) <- seq_len(nrow(measures))

    } else {

        alpha <- list()
        beta <- list()
        alpha[["1_signatures"]] <- rank0_alpha
        beta[["1_signatures"]] <- rank0_beta
        measures <- rank0_measures
        rownames(measures) <- seq_len(nrow(measures))

    }

    # save results
    results <- list(alpha = alpha, beta = beta, measures = measures)
    rm(alpha)
    rm(beta)
    rm(measures)
    gc(verbose = FALSE)

    # close pbackend
    if (close_parallel) {
        stopCluster(pbackend)
    }
    rm(close_parallel)
    rm(pbackend)
    gc(verbose = FALSE)

    # return the discovered signatures
    return(results)

}

#' Perform the assessment of different signaturesDecomposition solutions by cross-validation for K (beta, as estimated by 
#' signaturesDecomposition) somatic mutational signatures given a set of observations x and discovered signatures beta.
#'
#' @examples
#' data(background)
#' data(patients)
#' set.seed(12345)
#' sigs <- signaturesDecomposition(x = patients[seq_len(10),],
#'                                 K = 3:4,
#'                                 background_signature = background,
#'                                 nmf_runs = 2,
#'                                 sparsify = FALSE,
#'                                 num_processes = 1)
#' set.seed(12345)
#' res <- signaturesCV(x = patients[seq_len(10),],
#'                     beta = sigs$beta,
#'                     cross_validation_iterations = 2,
#'                     cross_validation_repetitions = 2,
#'                     num_processes = 1)
#'
#' @title signaturesCV
#' @param x Counts matrix for a set of n patients and m categories. These can be, e.g., SBS, MNV, CN or CN counts;
#' in the case of SBS it would be an n patients x 96 trinucleotides matrix.
#' @param beta A set of inferred signatures as returned by signaturesDecomposition function.
#' @param normalize_counts If true, the input counts matrix x is normalized such that the patients have the same number of mutation.
#' @param cross_validation_entries Percentage of cells in the counts matrix to be replaced by 0s during cross-validation.
#' @param cross_validation_iterations For each configuration, the first time the signatures are fitted form a matrix with a
#' percentage of values replaced by 0s. This may result in poor fit/results. Then, we perform predictions of these entries and replace them with
#' such predicted values. This parameter is the number of restarts to be performed to improve this estimate and obtain more stable solutions.
#' @param cross_validation_repetitions Number of time cross-validation should be repeated. Higher values result in better estimate, but are 
#' computationally more expensive.
#' @param num_processes Number of processes to be used during parallel execution. To execute in single process mode,
#' this parameter needs to be set to either NA or NULL.
#' @param verbose Boolean. Shall I print information messages?
#' @return A list of 2 elements: estimates and summary. Here, cv_estimates reports the mean squared error for each configuration of performed
#' cross-validation; rank_estimates reports mean and median values for each value of K.
#' @export signaturesCV
#' @import glmnet
#' @import parallel
#'
signaturesCV <- function( x, beta, normalize_counts = FALSE, cross_validation_entries = 0.01,
    cross_validation_iterations = 5, cross_validation_repetitions = 100, num_processes = Inf,
    verbose = TRUE ) {

    # check input parameters
    x <- as.matrix(x)
    if (normalize_counts) {
        x <- (x/rowSums(x)) * 2500
    }

    # setting up parallel execution
    parallel <- NULL
    close_parallel <- FALSE
    if (is.na(num_processes) || is.null(num_processes)) {
        num_processes <- 1
    } else if (num_processes == Inf) {
        cores <- as.integer((detectCores() - 1))
        num_processes <- min(cores, cross_validation_repetitions)
    } else {
        num_processes <- min(num_processes, cross_validation_repetitions)
    }

    if (verbose) {
        message("Estimating the optimal number of signatures with a total of ", cross_validation_repetitions,
            " cross-validation repetitions...", "\n")
        if (num_processes > 1) {
            message("Executing ", num_processes, " processes in parallel...", "\n")
        }
    }

    # structure to save the results
    cv_estimates <- matrix(NA, nrow=cross_validation_repetitions, ncol=length(beta))
    rownames(cv_estimates) <- paste0("Repetition_", seq_len(cross_validation_repetitions))
    colnames(cv_estimates) <- seq_len(ncol(cv_estimates))

    # perform a total of cross_validation_repetitions repetitions of cross-validation
    valid_entries <- which(x > 0, arr.ind = TRUE)

    # performing the inference
    if (num_processes == 1) { # performing the inference sequentially

        for (cv_repetitions in seq_len(cross_validation_repetitions)) {

            if (verbose) {
                message("Performing repetition ", cv_repetitions, " out of ",
                    cross_validation_repetitions, "...", "\n")
            }

            # randomly set the cross-validation entries for the current iteration
            cv_entries <- valid_entries[sample(seq_len(nrow(valid_entries)), 
                size = round(nrow(valid_entries) * cross_validation_entries), replace = FALSE), ]

            # consider all the possible values of K
            cont <- 0
            for (num_signs in seq_len(length(beta))) {

                k <- nrow(beta[[num_signs]])
                if (cv_repetitions == 1) {
                    colnames(cv_estimates)[num_signs] <- paste0(k, "_signatures")
                }

                if (verbose) {
                    message("Performing estimation for K=", k, "...", "\n")
                }

                # repeat the estimation for a number of cross_validation_iterations
                x_cv <- x
                for (cv_iteration in seq_len(cross_validation_iterations)) {

                    if (verbose) {
                        message("Performing cross-validation iteration ", cv_iteration,
                            " out of ", cross_validation_iterations, "...", "\n")
                    }

                    # set a percentage of cross_validation_entries entries to 0
                    # in order to perform cross-validation
                    if (cv_iteration == 1) {
                        x_cv[cv_entries] <- 0
                    } else {
                        predicted_counts <- curr_results[["alpha"]] %*% curr_results[["beta"]]
                        x_cv[cv_entries] <- predicted_counts[cv_entries]
                    }

                    # perform the inference
                    curr_results <- .nmf_fit(x = x_cv, beta = beta[[num_signs]])
                    gc(verbose = FALSE)

                }

                # save an estimate of the current solution
                curr_predicted_counts <- curr_results[["alpha"]] %*% curr_results[["beta"]]
                rm(curr_results)
                curr_true_considered_counts <- as.vector(x[cv_entries])
                curr_predicted_considered_counts <- as.vector(curr_predicted_counts[cv_entries])
                error <- median((curr_true_considered_counts - curr_predicted_considered_counts)^2)
                rm(curr_true_considered_counts)
                rm(curr_predicted_considered_counts)
                cv_estimates[cv_repetitions, paste0(k, "_signatures")] <- error
                rm(error)
                gc(verbose = FALSE)

                if (verbose) {
                    cont <- cont + 1
                    message("Progress ", round((cont/length(beta)) * 100, digits = 3), 
                        "%...", "\n")
                }

            }

        }

        # compute mean and median values of estimated cross-validation error
        cv_mean <- NULL
        cv_median <- NULL
        for (i in seq_len(ncol(cv_estimates))) {
            cv_mean <- c(cv_mean, mean(cv_estimates[, i]))
            cv_median <- c(cv_median, median(cv_estimates[, i]))
        }
        cv_summary <- cbind(cv_mean, cv_median)
        rownames(cv_summary) <- colnames(cv_estimates)
        colnames(cv_summary) <- c("Mean cross-validation error", "Median cross-validation error")

    } else { # performing the inference in parallel

        parallel <- makeCluster(num_processes, outfile = "")
        close_parallel <- TRUE
        res_clusterEvalQ <- clusterEvalQ(parallel, library("glmnet", warn.conflicts = FALSE,
            quietly = TRUE, verbose = FALSE))
        clusterExport(parallel, varlist = c(".nmf_fit", "verbose", "cross_validation_repetitions",
            "cross_validation_entries"), envir = environment())
        clusterExport(parallel, varlist = c("cross_validation_iterations", "valid_entries",
            "beta", "x"), envir = environment())
        clusterSetRNGStream(parallel, iseed = round(runif(1) * 1e+05))
        rm(res_clusterEvalQ)
        gc(verbose = FALSE)
        curr_results <- parLapply(parallel, seq_len(cross_validation_repetitions), function(cv_repetitions) {

            if (verbose) {
                message("Performing repetition ", cv_repetitions, " out of ",
                    cross_validation_repetitions, "...", "\n")
            }

            # randomly set the cross-validation entries for the current iteration
            cv_entries <- valid_entries[sample(seq_len(nrow(valid_entries)), 
                size = round(nrow(valid_entries) * cross_validation_entries), replace = FALSE), ]

            # consider all the possible values of K
            cv_errors <- rep(NA, length(beta))
            cv_names <- NULL
            for (num_signs in seq_len(length(beta))) {

                if (cv_repetitions == 1) {
                    cv_names <- c(cv_names, paste0(nrow(beta[[num_signs]]), "_signatures"))
                }

                # repeat the estimation for a number of cross_validation_iterations
                x_cv <- x
                for (cv_iteration in seq_len(cross_validation_iterations)) {

                    # set a percentage of cross_validation_entries entries to 0
                    # in order to perform cross-validation
                    if (cv_iteration == 1) {
                        x_cv[cv_entries] <- 0
                    } else {
                        x_cv[cv_entries] <- predicted_counts[cv_entries]
                    }

                    # perform the inference
                    curr_results <- .nmf_fit(x = x_cv, beta = beta[[num_signs]])
                    gc(verbose = FALSE)
                }

                # save an estimate of the current solution
                curr_predicted_counts <- curr_results[["alpha"]] %*% curr_results[["beta"]]
                rm(curr_results)
                curr_true_considered_counts <- as.vector(x[cv_entries])
                curr_predicted_considered_counts <- as.vector(curr_predicted_counts[cv_entries])
                error <- median((curr_true_considered_counts - curr_predicted_considered_counts)^2)
                rm(curr_true_considered_counts)
                rm(curr_predicted_considered_counts)
                cv_errors[num_signs] <- error
                rm(error)
                gc(verbose = FALSE)

            }
            names(cv_errors) <- cv_names

            return(cv_errors)

        })

        # save the results from parallel computation
        colnames(cv_estimates) <- names(curr_results[[1]])
        for (par_res in seq_len(length(curr_results))) {
            cv_estimates[par_res, ] <- curr_results[[par_res]]
        }
        rm(curr_results)
        gc(verbose = FALSE)

        # compute mean and median values of estimated cross-validation error
        cv_mean <- NULL
        cv_median <- NULL
        for (i in seq_len(ncol(cv_estimates))) {
            cv_mean <- c(cv_mean, mean(cv_estimates[, i]))
            cv_median <- c(cv_median, median(cv_estimates[, i]))
        }
        cv_summary <- cbind(cv_mean, cv_median)
        rownames(cv_summary) <- colnames(cv_estimates)
        colnames(cv_summary) <- c("Mean cross-validation error", "Median cross-validation error")

    }

    # save results
    results <- list(estimates = cv_estimates, summary = cv_summary)
    rm(cv_estimates)
    rm(cv_summary)
    gc(verbose = FALSE)

    # close parallel
    if (close_parallel) {
        stopCluster(parallel)
    }
    rm(close_parallel)
    rm(parallel)
    gc(verbose = FALSE)

    # return results of cross-validation
    return(results)

}

# initialize alpha and beta for the .nmf_ls and .nmf_reg methods
.nmf_seed <- function( model, target ) {

    # initialize alpha with an empty matrix
    alpha <- matrix(NA, nrow=nrow(target), ncol=nbasis(model))
    rownames(alpha) <- rownames(target)
    colnames(alpha) <- paste0("S", seq_len(ncol(alpha)))

    # randomly initialize beta
    beta <- matrix(runif(n = nbasis(model) * ncol(target), min = 0, max = 1),
        nrow = nbasis(model), ncol = ncol(target))
    beta <- (beta/rowSums(beta))  # beta rows (i.e., signatures) must sum to 1
    if (any(is.nan(rowSums(beta)))) {
        beta[which(is.nan(rowSums(beta))), ] <- (1/ncol(beta))
    }
    rownames(beta) <- colnames(alpha)
    colnames(beta) <- colnames(target)

    # update results
    basis(model) <- alpha
    coef(model) <- beta
    rm(alpha)
    rm(beta)
    gc(verbose = FALSE)

    # return the updated model
    return(model)

}

# objective function for the .nmf_ls and .nmf_reg methods
.nmf_objective <- function( x, y, background ) {

    # initialization
    predictions <- basis(x) %*% coef(x)
    rm(x)
    rm(background)
    gc(verbose=FALSE)

    # compute mean cosine similarity comparing observations and predictions
    obj_value <- rep(NA, nrow(y))
    for(i in 1:nrow(y)) {
        obj_value[i] <- as.numeric(cosine(as.numeric(y[i,]),as.numeric(predictions[i,])))
    }
    obj_value <- (1 - mean(obj_value,na.rm=TRUE))
    rm(y)
    rm(predictions)
    gc(verbose=FALSE)

    # return the estimated objective value
    return(obj_value)

}

# perform NMF by Non-negative least squares
.nmf_ls <- function( x, seed, background ) {

    # initialization
    alpha <- basis(seed)  # exposures matrix
    beta <- coef(seed)  # signatures matrix
    if (!is.null(background)) {
        # if specified, the first row of beta is the background signature
        beta[1, ] <- background
    }
    n <- nrow(x)  # n is the number of observations in x, i.e., the patients
    J <- ncol(x)  # J is the number of categories, e.g., 96 trinucleotides
    K <- nrow(beta)  # K is the number of signatures to be fitted

    # iteratively fit alpha and beta by Non-negative least squares
    for (i in seq_len(20)) {
        # update alpha, beta is kept fixed
        for (j in seq_len(n)) {
            curr_beta_values <- beta
            if (nrow(curr_beta_values) > 1) {
                alpha[j, ] <- tryCatch({
                    res <- glmnet(x = t(curr_beta_values), y = as.vector(x[j, ]), 
                        lambda = 0, family = "gaussian", lower.limits = 0)
                    res <- as.numeric(coef(res,s=res$lambda))
                    res <- as.numeric(res[1]+res[-1])
                    if(any(res<0)) {
                        res[which(res<0)] <- 0
                    }
                    res
                }, error = function( e ) {
                    return((sum(as.vector(x[j,]))/ncol(alpha)))
                }, finally = {
                    gc(verbose = FALSE)
                })
            } else {
                res_inputs <- cbind(t(curr_beta_values), rep(0, ncol(curr_beta_values)))
                alpha[j, ] <- tryCatch({
                    res <- glmnet(x = res_inputs, y = as.vector(x[j, ]), 
                        lambda = 0, family = "gaussian", lower.limits = 0)
                    res <- as.numeric(coef(res,s=res$lambda))
                    res <- as.numeric(res[1]+res[-1])[-length(res)]
                    if(any(res<0)) {
                        res[which(res<0)] <- 0
                    }
                    res
                }, error = function( e ) {
                    return((sum(as.vector(x[j,]))/ncol(alpha)))
                }, finally = {
                    gc(verbose = FALSE)
                })
            }
        }
        # update beta, alpha is kept fixed
        for (k in seq_len(J)) {
            if (!is.null(background)) {
                # the first signature represents the background model, thus it
                # is not changed during the fit
                curr_alpha_values <- alpha[, 2:K, drop = FALSE]
                if (ncol(curr_alpha_values) > 1) {
                    beta[2:K, k] <- tryCatch({
                        res <- glmnet(x = curr_alpha_values, y = as.vector(x[, k]), 
                          lambda = 0, family = "gaussian", lower.limits = 0)
                        res <- as.numeric(coef(res,s=res$lambda))
                        res <- as.numeric(res[-1])
                        res
                    }, error = function( e ) {
                        return(0)
                    }, finally = {
                        gc(verbose = FALSE)
                    })
                } else {
                    res_inputs <- cbind(curr_alpha_values, rep(0, nrow(curr_alpha_values)))
                    beta[2:K, k] <- tryCatch({
                        res <- glmnet(x = res_inputs, y = as.vector(x[, k]), 
                          lambda = 0, family = "gaussian", lower.limits = 0)
                        res <- as.numeric(coef(res,s=res$lambda))
                        res <- as.numeric(res[-1])
                        res
                    }, error = function( e ) {
                        return(0)
                    }, finally = {
                        gc(verbose = FALSE)
                    })
                }
            } else {
                curr_alpha_values <- alpha
                if (ncol(curr_alpha_values) > 1) {
                    beta[, k] <- tryCatch({
                        res <- glmnet(x = curr_alpha_values, y = as.vector(x[, k]), 
                          lambda = 0, family = "gaussian", lower.limits = 0)
                        res <- as.numeric(coef(res,s=res$lambda))
                        res <- as.numeric(res[-1])
                        res
                    }, error = function( e ) {
                        return(0)
                    }, finally = {
                        gc(verbose = FALSE)
                    })
                } else {
                    res_inputs <- cbind(curr_alpha_values, rep(0, nrow(curr_alpha_values)))
                    beta[, k] <- tryCatch({
                        res <- glmnet(x = res_inputs, y = as.vector(x[, k]), 
                          lambda = 0, family = "gaussian", lower.limits = 0)
                        res <- as.numeric(coef(res,s=res$lambda))
                        res <- as.numeric(res[-1])
                        res
                    }, error = function( e ) {
                        return(0)
                    }, finally = {
                        gc(verbose = FALSE)
                    })
                }
            }
        }
    }

    # normalize beta and perform final fit of alpha
    beta <- (beta/rowSums(beta))
    if (any(is.nan(rowSums(beta)))) {
        beta[which(is.nan(rowSums(beta))), ] <- (1/ncol(beta))
    }
    for (j in seq_len(n)) {
        curr_beta_values <- beta
        if (nrow(curr_beta_values) > 1) {
            alpha[j, ] <- tryCatch({
                res <- glmnet(x = t(curr_beta_values), y = as.vector(x[j, ]), 
                    lambda = 0, family = "gaussian", lower.limits = 0)
                res <- as.numeric(coef(res,s=res$lambda))
                res <- as.numeric(res[1]+res[-1])
                if(any(res<0)) {
                    res[which(res<0)] <- 0
                }
                res
            }, error = function( e ) {
                return((sum(as.vector(x[j,]))/ncol(alpha)))
            }, finally = {
                gc(verbose = FALSE)
            })
        } else {
            res_inputs <- cbind(t(curr_beta_values), rep(0, ncol(curr_beta_values)))
            alpha[j, ] <- tryCatch({
                res <- glmnet(x = res_inputs, y = as.vector(x[j, ]), 
                    lambda = 0, family = "gaussian", lower.limits = 0)
                res <- as.numeric(coef(res,s=res$lambda))
                res <- as.numeric(res[1]+res[-1])[-length(res)]
                if(any(res<0)) {
                    res[which(res<0)] <- 0
                }
                res
            }, error = function( e ) {
                return((sum(as.vector(x[j,]))/ncol(alpha)))
            }, finally = {
                gc(verbose = FALSE)
            })
        }
    }
    if (!is.null(background)) {
        colnames(alpha) <- c("Background", paste0("S", seq_len((ncol(alpha) - 1))))
    } else {
        colnames(alpha) <- paste0("S", seq_len(ncol(alpha)))
    }
    rownames(beta) <- colnames(alpha)

    # update the results
    basis(seed) <- alpha
    coef(seed) <- beta
    rm(x)
    rm(alpha)
    rm(beta)
    gc(verbose = FALSE)

    # return the updated seed
    return(seed)

}

# perform NMF by Non-negative fit regularized by elastic net with LASSO penalty
.nmf_reg <- function( x, seed, background ) {

    # initialization
    alpha <- basis(seed)  # exposures matrix
    beta <- coef(seed)  # signatures matrix
    if (!is.null(background)) {
        # if specified, the first row of beta is the background signature
        beta[1, ] <- background
    }
    n <- nrow(x)  # n is the number of observations in x, i.e., the patients
    J <- ncol(x)  # J is the number of categories, e.g., 96 trinucleotides
    K <- nrow(beta)  # K is the number of signatures to be fitted

    # STEP 1: we get close to a good solution by Non-negative least squares

    # iteratively fit alpha and beta by Non-negative least squares
    for (i in seq_len(20)) {
        # update alpha, beta is kept fixed
        for (j in seq_len(n)) {
            curr_beta_values <- beta
            if (nrow(curr_beta_values) > 1) {
                alpha[j, ] <- tryCatch({
                    res <- glmnet(x = t(curr_beta_values), y = as.vector(x[j, ]), 
                        lambda = 0, family = "gaussian", lower.limits = 0)
                    res <- as.numeric(coef(res,s=res$lambda))
                    res <- as.numeric(res[1]+res[-1])
                    if(any(res<0)) {
                        res[which(res<0)] <- 0
                    }
                    res
                }, error = function( e ) {
                    return((sum(as.vector(x[j,]))/ncol(alpha)))
                }, finally = {
                    gc(verbose = FALSE)
                })
            } else {
                res_inputs <- cbind(t(curr_beta_values), rep(0, ncol(curr_beta_values)))
                alpha[j, ] <- tryCatch({
                    res <- glmnet(x = res_inputs, y = as.vector(x[j, ]), 
                        lambda = 0, family = "gaussian", lower.limits = 0)
                    res <- as.numeric(coef(res,s=res$lambda))
                    res <- as.numeric(res[1]+res[-1])[-length(res)]
                    if(any(res<0)) {
                        res[which(res<0)] <- 0
                    }
                    res
                }, error = function( e ) {
                    return((sum(as.vector(x[j,]))/ncol(alpha)))
                }, finally = {
                    gc(verbose = FALSE)
                })
            }
        }
        # update beta, alpha is kept fixed
        for (k in seq_len(J)) {
            if (!is.null(background)) {
                # the first signature represents the background model, thus it
                # is not changed during the fit
                curr_alpha_values <- alpha[, 2:K, drop = FALSE]
                if (ncol(curr_alpha_values) > 1) {
                    beta[2:K, k] <- tryCatch({
                        res <- glmnet(x = curr_alpha_values, y = as.vector(x[, k]), 
                          lambda = 0, family = "gaussian", lower.limits = 0)
                        res <- as.numeric(coef(res,s=res$lambda))
                        res <- as.numeric(res[-1])
                        res
                    }, error = function( e ) {
                        return(0)
                    }, finally = {
                        gc(verbose = FALSE)
                    })
                } else {
                    res_inputs <- cbind(curr_alpha_values, rep(0, nrow(curr_alpha_values)))
                    beta[2:K, k] <- tryCatch({
                        res <- glmnet(x = res_inputs, y = as.vector(x[, k]), 
                          lambda = 0, family = "gaussian", lower.limits = 0)
                        res <- as.numeric(coef(res,s=res$lambda))
                        res <- as.numeric(res[-1])
                        res
                    }, error = function( e ) {
                        return(0)
                    }, finally = {
                        gc(verbose = FALSE)
                    })
                }
            } else {
                curr_alpha_values <- alpha
                if (ncol(curr_alpha_values) > 1) {
                    beta[, k] <- tryCatch({
                        res <- glmnet(x = curr_alpha_values, y = as.vector(x[, k]), 
                          lambda = 0, family = "gaussian", lower.limits = 0)
                        res <- as.numeric(coef(res,s=res$lambda))
                        res <- as.numeric(res[-1])
                        res
                    }, error = function( e ) {
                        return(0)
                    }, finally = {
                        gc(verbose = FALSE)
                    })
                } else {
                    res_inputs <- cbind(curr_alpha_values, rep(0, nrow(curr_alpha_values)))
                    beta[, k] <- tryCatch({
                        res <- glmnet(x = res_inputs, y = as.vector(x[, k]), 
                          lambda = 0, family = "gaussian", lower.limits = 0)
                        res <- as.numeric(coef(res,s=res$lambda))
                        res <- as.numeric(res[-1])
                        res
                    }, error = function( e ) {
                        return(0)
                    }, finally = {
                        gc(verbose = FALSE)
                    })
                }
            }
        }
    }

    # normalize beta and perform final fit of alpha
    beta <- (beta/rowSums(beta))
    if (any(is.nan(rowSums(beta)))) {
        beta[which(is.nan(rowSums(beta))), ] <- (1/ncol(beta))
    }
    for (j in seq_len(n)) {
        curr_beta_values <- beta
        if (nrow(curr_beta_values) > 1) {
            alpha[j, ] <- tryCatch({
                res <- glmnet(x = t(curr_beta_values), y = as.vector(x[j, ]), 
                    lambda = 0, family = "gaussian", lower.limits = 0)
                res <- as.numeric(coef(res,s=res$lambda))
                res <- as.numeric(res[1]+res[-1])
                if(any(res<0)) {
                    res[which(res<0)] <- 0
                }
                res
            }, error = function( e ) {
                return((sum(as.vector(x[j,]))/ncol(alpha)))
            }, finally = {
                gc(verbose = FALSE)
            })
        } else {
            res_inputs <- cbind(t(curr_beta_values), rep(0, ncol(curr_beta_values)))
            alpha[j, ] <- tryCatch({
                res <- glmnet(x = res_inputs, y = as.vector(x[j, ]), 
                    lambda = 0, family = "gaussian", lower.limits = 0)
                res <- as.numeric(coef(res,s=res$lambda))
                res <- as.numeric(res[1]+res[-1])[-length(res)]
                if(any(res<0)) {
                    res[which(res<0)] <- 0
                }
                res
            }, error = function( e ) {
                return((sum(as.vector(x[j,]))/ncol(alpha)))
            }, finally = {
                gc(verbose = FALSE)
            })
        }
    }
    if (!is.null(background)) {
        colnames(alpha) <- c("Background", paste0("S", seq_len((ncol(alpha) - 1))))
    } else {
        colnames(alpha) <- paste0("S", seq_len(ncol(alpha)))
    }
    rownames(beta) <- colnames(alpha)

    # STEP 2: we finalize the inference by regularizing with elastic net using LASSO penalty

    # now regularize solutions by Non-negative fit with elastic net with LASSO penalty
    for (i in seq_len(10)) {
        # update alpha, beta is kept fixed
        for (j in seq_len(n)) {
            curr_beta_values <- beta
            if (nrow(curr_beta_values) > 1) {
                alpha[j, ] <- tryCatch({
                    res <- cv.glmnet(x = t(curr_beta_values), y = as.vector(x[j, ]), 
                        type.measure = "mse", nfolds = 10, nlambda = 100, 
                        family = "gaussian", lower.limits = 0)
                    res <- as.numeric(coef(res,s=res$lambda.min))
                    res <- as.numeric(res[1]+res[-1])
                    if(any(res<0)) {
                        res[which(res<0)] <- 0
                    }
                    res
                }, error = function( e ) {
                    return((sum(as.vector(x[j,]))/ncol(alpha)))
                }, finally = {
                    gc(verbose = FALSE)
                })
            } else {
                res_inputs <- cbind(t(curr_beta_values), rep(0, ncol(curr_beta_values)))
                alpha[j, ] <- tryCatch({
                    res <- cv.glmnet(x = res_inputs, y = as.vector(x[j, ]), 
                        type.measure = "mse", nfolds = 10, nlambda = 100, 
                        family = "gaussian", lower.limits = 0)
                    res <- as.numeric(coef(res,s=res$lambda.min))
                    res <- as.numeric(res[1]+res[-1])[-length(res)]
                    if(any(res<0)) {
                        res[which(res<0)] <- 0
                    }
                    res
                }, error = function( e ) {
                    return((sum(as.vector(x[j,]))/ncol(alpha)))
                }, finally = {
                    gc(verbose = FALSE)
                })
            }
        }
        # update beta, alpha is kept fixed
        for (k in seq_len(J)) {
            if (!is.null(background)) {
                # the first signature represents the background model, thus it
                # is not changed during the fit
                curr_alpha_values <- alpha[, 2:K, drop = FALSE]
                if (ncol(curr_alpha_values) > 1) {
                    beta[2:K, k] <- tryCatch({
                        res <- cv.glmnet(x = curr_alpha_values, y = as.vector(x[, k]), 
                            type.measure = "mse", nfolds = 10, nlambda = 100, 
                            family = "gaussian", lower.limits = 0)
                        res <- as.numeric(coef(res,s=res$lambda.min))
                        res <- as.numeric(res[-1])
                        res
                    }, error = function( e ) {
                        return(0)
                    }, finally = {
                        gc(verbose = FALSE)
                    })
                } else {
                    res_inputs <- cbind(curr_alpha_values, rep(0, nrow(curr_alpha_values)))
                    beta[2:K, k] <- tryCatch({
                        res <- cv.glmnet(x = res_inputs, y = as.vector(x[, k]), 
                            type.measure = "mse", nfolds = 10, nlambda = 100, 
                            family = "gaussian", lower.limits = 0)
                        res <- as.numeric(coef(res,s=res$lambda.min))
                        res <- as.numeric(res[-1])
                        res
                    }, error = function( e ) {
                        return(0)
                    }, finally = {
                        gc(verbose = FALSE)
                    })
                }
            } else {
                curr_alpha_values <- alpha
                if (ncol(curr_alpha_values) > 1) {
                    beta[, k] <- tryCatch({
                        res <- cv.glmnet(x = curr_alpha_values, y = as.vector(x[, k]), 
                            type.measure = "mse", nfolds = 10, nlambda = 100, 
                            family = "gaussian", lower.limits = 0)
                        res <- as.numeric(coef(res,s=res$lambda.min))
                        res <- as.numeric(res[-1])
                        res
                    }, error = function( e ) {
                        return(0)
                    }, finally = {
                        gc(verbose = FALSE)
                    })
                } else {
                    res_inputs <- cbind(curr_alpha_values, rep(0, nrow(curr_alpha_values)))
                    beta[, k] <- tryCatch({
                        res <- cv.glmnet(x = res_inputs, y = as.vector(x[, k]), 
                            type.measure = "mse", nfolds = 10, nlambda = 100, 
                            family = "gaussian", lower.limits = 0)
                        res <- as.numeric(coef(res,s=res$lambda.min))
                        res <- as.numeric(res[-1])
                        res
                    }, error = function( e ) {
                        return(0)
                    }, finally = {
                        gc(verbose = FALSE)
                    })
                }
            }
        }
    }

    # normalize beta and perform final fit of alpha
    beta <- (beta/rowSums(beta))
    if (any(is.nan(rowSums(beta)))) {
        beta[which(is.nan(rowSums(beta))), ] <- (1/ncol(beta))
    }
    for (j in seq_len(n)) {
        curr_beta_values <- beta
        if (nrow(curr_beta_values) > 1) {
            alpha[j, ] <- tryCatch({
                res <- cv.glmnet(x = t(curr_beta_values), y = as.vector(x[j, ]), 
                    type.measure = "mse", nfolds = 10, nlambda = 100, 
                    family = "gaussian", lower.limits = 0)
                res <- as.numeric(coef(res,s=res$lambda.min))
                res <- as.numeric(res[1]+res[-1])
                if(any(res<0)) {
                    res[which(res<0)] <- 0
                }
                res
            }, error = function( e ) {
                return((sum(as.vector(x[j,]))/ncol(alpha)))
            }, finally = {
                gc(verbose = FALSE)
            })
        } else {
            res_inputs <- cbind(t(curr_beta_values), rep(0, ncol(curr_beta_values)))
            alpha[j, ] <- tryCatch({
                res <- cv.glmnet(x = res_inputs, y = as.vector(x[j, ]), 
                    type.measure = "mse", nfolds = 10, nlambda = 100, 
                    family = "gaussian", lower.limits = 0)
                res <- as.numeric(coef(res,s=res$lambda.min))
                res <- as.numeric(res[1]+res[-1])[-length(res)]
                if(any(res<0)) {
                    res[which(res<0)] <- 0
                }
                res
            }, error = function( e ) {
                return((sum(as.vector(x[j,]))/ncol(alpha)))
            }, finally = {
                gc(verbose = FALSE)
            })
        }
    }
    if (!is.null(background)) {
        colnames(alpha) <- c("Background", paste0("S", seq_len((ncol(alpha) - 1))))
    } else {
        colnames(alpha) <- paste0("S", seq_len(ncol(alpha)))
    }
    rownames(beta) <- colnames(alpha)

    # update the results
    basis(seed) <- alpha
    coef(seed) <- beta
    rm(x)
    rm(alpha)
    rm(beta)
    gc(verbose = FALSE)

    # return the updated seed
    return(seed)

}

# perform fit of a given NMF solution by Non-negative least squares
.nmf_fit <- function( x, beta ) {

    # initialization
    alpha <- matrix(NA, nrow=nrow(x),ncol=nrow(beta))
    rownames(alpha) <- rownames(x)
    colnames(alpha) <- rownames(beta)
    beta <- (beta/rowSums(beta))
    if (any(is.nan(rowSums(beta)))) {
        beta[which(is.nan(rowSums(beta))), ] <- (1/ncol(beta))
    }
    n <- nrow(x)  # n is the number of observations in x, i.e., the patients
    J <- ncol(x)  # J is the number of categories, e.g., 96 trinucleotides

    # iteratively fit alpha and beta by Non-negative least squares
    for (i in seq_len(20)) {
        # update alpha, beta is kept fixed
        for (j in seq_len(n)) {
            curr_beta_values <- beta
            if (nrow(curr_beta_values) > 1) {
                alpha[j, ] <- tryCatch({
                    res <- glmnet(x = t(curr_beta_values), y = as.vector(x[j, ]), 
                        lambda = 0, family = "gaussian", lower.limits = 0)
                    res <- as.numeric(coef(res,s=res$lambda))
                    res <- as.numeric(res[1]+res[-1])
                    if(any(res<0)) {
                        res[which(res<0)] <- 0
                    }
                    res
                }, error = function( e ) {
                    return((sum(as.vector(x[j,]))/ncol(alpha)))
                }, finally = {
                    gc(verbose = FALSE)
                })
            } else {
                res_inputs <- cbind(t(curr_beta_values), rep(0, ncol(curr_beta_values)))
                alpha[j, ] <- tryCatch({
                    res <- glmnet(x = res_inputs, y = as.vector(x[j, ]), 
                        lambda = 0, family = "gaussian", lower.limits = 0)
                    res <- as.numeric(coef(res,s=res$lambda))
                    res <- as.numeric(res[1]+res[-1])[-length(res)]
                    if(any(res<0)) {
                        res[which(res<0)] <- 0
                    }
                    res
                }, error = function( e ) {
                    return((sum(as.vector(x[j,]))/ncol(alpha)))
                }, finally = {
                    gc(verbose = FALSE)
                })
            }
        }
        # update beta, alpha is kept fixed
        for (k in seq_len(J)) {
            curr_alpha_values <- alpha
            if (ncol(curr_alpha_values) > 1) {
                beta[, k] <- tryCatch({
                    res <- glmnet(x = curr_alpha_values, y = as.vector(x[, k]), 
                      lambda = 0, family = "gaussian", lower.limits = 0)
                    res <- as.numeric(coef(res,s=res$lambda))
                    res <- as.numeric(res[-1])
                    res
                }, error = function( e ) {
                    return(0)
                }, finally = {
                    gc(verbose = FALSE)
                })
            } else {
                res_inputs <- cbind(curr_alpha_values, rep(0, nrow(curr_alpha_values)))
                beta[, k] <- tryCatch({
                    res <- glmnet(x = res_inputs, y = as.vector(x[, k]), 
                      lambda = 0, family = "gaussian", lower.limits = 0)
                    res <- as.numeric(coef(res,s=res$lambda))
                    res <- as.numeric(res[-1])
                    res
                }, error = function( e ) {
                    return(0)
                }, finally = {
                    gc(verbose = FALSE)
                })
            }
        }
    }

    # normalize beta and perform final fit of alpha
    beta <- (beta/rowSums(beta))
    if (any(is.nan(rowSums(beta)))) {
        beta[which(is.nan(rowSums(beta))), ] <- (1/ncol(beta))
    }
    for (j in seq_len(n)) {
        curr_beta_values <- beta
        if (nrow(curr_beta_values) > 1) {
            alpha[j, ] <- tryCatch({
                res <- glmnet(x = t(curr_beta_values), y = as.vector(x[j, ]), 
                    lambda = 0, family = "gaussian", lower.limits = 0)
                res <- as.numeric(coef(res,s=res$lambda))
                res <- as.numeric(res[1]+res[-1])
                if(any(res<0)) {
                    res[which(res<0)] <- 0
                }
                res
            }, error = function( e ) {
                return((sum(as.vector(x[j,]))/ncol(alpha)))
            }, finally = {
                gc(verbose = FALSE)
            })
        } else {
            res_inputs <- cbind(t(curr_beta_values), rep(0, ncol(curr_beta_values)))
            alpha[j, ] <- tryCatch({
                res <- glmnet(x = res_inputs, y = as.vector(x[j, ]), 
                    lambda = 0, family = "gaussian", lower.limits = 0)
                res <- as.numeric(coef(res,s=res$lambda))
                res <- as.numeric(res[1]+res[-1])[-length(res)]
                if(any(res<0)) {
                    res[which(res<0)] <- 0
                }
                res
            }, error = function( e ) {
                return((sum(as.vector(x[j,]))/ncol(alpha)))
            }, finally = {
                gc(verbose = FALSE)
            })
        }
    }

    # update the results
    results <- list(alpha = alpha, beta = beta)
    rm(x)
    rm(alpha)
    rm(beta)
    gc(verbose = FALSE)

    # return the results
    return(results)

}
