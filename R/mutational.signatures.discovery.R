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
#' beta <- signaturesDecomposition(x = patients[seq_len(3),seq_len(2)],
#'                                 K = 3,
#'                                 background_signature = background[seq_len(2)],
#'                                 nmf_runs = 2,
#'                                 num_processes = 1)
#' set.seed(12345)
#' res <- signaturesAssignment(x = patients[seq_len(3),seq_len(2)], beta = beta$beta[[1]])
#'
#' @title signaturesAssignment
#' @param x Counts matrix for a set of n patients and m categories. These can be, e.g., SBS, MNV, CN or CN counts;
#' in the case of SBS it would be an n patients x 96 trinucleotides matrix.
#' @param beta Matrix of the discovered signatures to be used for the assignment.
#' @return A list with the discovered signatures. It includes 3 elements:
#'              alpha: matrix of the discovered exposure values.
#'              beta: matrix of the discovered signatures.
#'              unexplained_mutations: number of unexplained mutations per sample.
#' @export signaturesAssignment
#' @import glmnet
#'
signaturesAssignment <- function( x, beta ) {

    # check input parameters
    x <- as.matrix(x)

    # initialize alpha with an empty matrix
    alpha <- matrix(NA, nrow=nrow(x), ncol=nrow(beta))
    rownames(alpha) <- rownames(x)
    colnames(alpha) <- rownames(beta)

    # perform signatures assignment
    unexplained_mutations <- rep(NA, nrow(x))
    names(unexplained_mutations) <- rownames(x)
    for (j in seq_len(nrow(x))) {
        if (nrow(beta) > 1) {
            alpha[j, ] <- tryCatch({
                res <- cv.glmnet(x = t(beta), y = as.vector(x[j, ]), 
                        type.measure = "mse", nfolds = 10, nlambda = 100, 
                        family = "gaussian", alpha = 1, lower.limits = 0, 
                        maxit = 1e+05)
                res <- as.numeric(coef(res,s=res$lambda.min))
                unexplained_mutations[j] <- res[1]
                res <- res[-1]
                res
            }, error = function( e ) {
                res <- (sum(x[j,])/ncol(alpha))
                return(res)
            }, finally = {
                gc(verbose = FALSE)
            })
        } else {
            fit_inputs <- cbind(t(beta), rep(0, ncol(beta)))
            alpha[j, ] <- tryCatch({
                res <- cv.glmnet(x = fit_inputs, y = as.vector(x[j, ]), 
                        type.measure = "mse", nfolds = 10, nlambda = 100, 
                        family = "gaussian", alpha = 1, lower.limits = 0, 
                        maxit = 1e+05)
                res <- as.numeric(coef(res,s=res$lambda.min))
                unexplained_mutations[j] <- res[1]
                res <- res[-length(res)]
                res <- res[-1]
                res
            }, error = function( e ) {
                res <- (sum(x[j,])/ncol(alpha))
                return(res)
            }, finally = {
                gc(verbose = FALSE)
            })
        }
    }

    # save results
    results <- list(alpha = alpha, beta = beta, unexplained_mutations = unexplained_mutations)
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
#' res <- signaturesDecomposition(x = patients[seq_len(3),seq_len(2)],
#'                                K = 3:4,
#'                                background_signature = background[seq_len(2)],
#'                                nmf_runs = 2,
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
#' @param num_processes Number of processes to be used during parallel execution. To execute in single process mode,
#' this parameter needs to be set to either NA or NULL.
#' @param verbose Boolean. Shall I print information messages?
#' @return A list with the discovered signatures and related rank measures. It includes 5 elements:
#'              alpha: list of matrices of the discovered exposure values for each possible rank in the range K.
#'              beta: list of matrices of the discovered signatures for each possible rank in the range K.
#'              unexplained_mutations: number of unexplained mutations per sample.
#'              cosine_similarity: cosine similarity comparing input data x and predictions for each rank in the range K.
#'              measures: a data.frame containing the quality measures for each possible rank in the range K.
#' @export signaturesDecomposition
#' @import glmnet
#' @import lsa
#' @import nnls
#' @import RhpcBLASctl
#' @import parallel
#'
signaturesDecomposition <- function( x, K, background_signature = NULL,
    normalize_counts = FALSE, nmf_runs = 100, num_processes = Inf, verbose = TRUE ) {

    # check input parameters
    x <- as.matrix(x)
    if (length(K) == 0) {
        warning("No valid value of K is provided: setting K equals to 1...")
        K <- 1
    }
    K <- sort(unique(K))
    if (min(K) < 1) {
        warning("The minimum value of K must be 1, removing invalid values...")
        K <- K[K >= 1]
        if (length(K) == 0) {
            K <- 1
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
    parallel <- NULL
    close_parallel <- FALSE
    if (is.na(num_processes) || is.null(num_processes)) {
        num_processes <- 1
    } else if (num_processes == Inf) {
        cores <- as.integer((detectCores() - 1))
        num_processes <- min(cores, nmf_runs)
    } else {
        num_processes <- min(num_processes, nmf_runs)
    }
    if(num_processes > 1) {
        parallel <- makeCluster(num_processes, outfile = "")
        close_parallel <- TRUE
        res_clusterEvalQ <- clusterEvalQ(parallel, library("glmnet", warn.conflicts = FALSE,
            quietly = TRUE, verbose = FALSE))
        res_clusterEvalQ <- clusterEvalQ(parallel, library("lsa", warn.conflicts = FALSE,
            quietly = TRUE, verbose = FALSE))
        res_clusterEvalQ <- clusterEvalQ(parallel, library("nnls", warn.conflicts = FALSE,
            quietly = TRUE, verbose = FALSE))
        clusterExport(parallel, varlist = c(".fit_nmf", ".fit_seed", ".fit_regularized",
            ".fit_objective"), envir = environment())
        clusterExport(parallel, varlist = c("x", "background_signature",
            "verbose", "nmf_runs"), envir = environment())
        clusterSetRNGStream(parallel, iseed = round(runif(1) * 1e+05))
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
        rank0_beta <- matrix(background_signature, nrow = 1, ncol = ncol(x))
        rownames(rank0_beta) <- "S1"
        colnames(rank0_beta) <- colnames(x)
        rank0_res <- signaturesAssignment(x = x, beta = rank0_beta)
        rank0_alpha <- rank0_res$alpha
        rank0_unexplained_mutations <- rank0_res$unexplained_mutations
        rank0_predictions <- ((rank0_alpha %*% rank0_beta) + rank0_unexplained_mutations)
        rank0_mse <- mean(((x - rank0_predictions)^2), na.rm = TRUE)
        rank0_obj <- .fit_objective(x = x, model = list(alpha = rank0_alpha, beta = rank0_beta, 
            unexplained_mutations = rank0_unexplained_mutations))
        rank0_cosine_similarity <- rank0_obj$cosine_similarities
        rank0_measures <- matrix(NA, nrow = 1, ncol = 4)
        rownames(rank0_measures) <- paste0("K=",1)
        colnames(rank0_measures) <- c("Stability", "Mean Squared Error", "Explained Variance", "Predictions")
        rank0_measures[1, "Stability"] <- 1
        rank0_measures[1, "Mean Squared Error"] <- rank0_mse
        pred_ss <- sum(((x - rank0_predictions)^2), na.rm = TRUE)
        total_ss <- sum((x^2), na.rm = TRUE)
        rank0_measures[1, "Explained Variance"] <- (1 - (pred_ss / total_ss))
        rank0_measures[1, "Predictions"] <- rank0_obj$goodness_fit
        if (normalize_counts) {
            rank0_alpha <- rank0_alpha * (rowSums(x_not_normalized)/2500)
            rank0_unexplained_mutations <- rank0_unexplained_mutations * (rowSums(x_not_normalized)/2500)
        }
    }

    # perform signatures discovery and rank estimation
    if (length(K) > 0) {

        alpha <- list()
        beta <- list()
        unexplained_mutations <- list()
        cosine_similarity <- list()
        measures <- NULL
        if (rank0) {
            alpha[["1_signatures"]] <- rank0_alpha
            beta[["1_signatures"]] <- rank0_beta
            unexplained_mutations[["1_signatures"]] <- rank0_unexplained_mutations
            cosine_similarity[["1_signatures"]] <- rank0_cosine_similarity
            measures <- rbind(measures,rank0_measures)
        }

        for (i in seq_len(length(K))) {

            rank = K[i]
            if (verbose) {
                message("Performing inference for K=", rank, "...", "\n")
            }

            # perform the inference for the current K
            if (num_processes == 1) { # performing the inference sequentially
                results <- list()
                for (iteration in seq_len(nmf_runs)) {
                    if (verbose) {
                        message("Performing NMF run ", iteration, " out of ",
                            nmf_runs, "...", "\n")
                    }
                    results[[iteration]] <- .fit_nmf(x = x, K = rank, background = background_signature)
                }
            } else { # performing the inference in parallel
                clusterExport(parallel, varlist = c("rank"), envir = environment())
                results <- parLapply(parallel, seq_len(nmf_runs), function(iteration) {
                    if (verbose) {
                        message("Performing NMF run ", iteration, " out of ",
                            nmf_runs, "...", "\n")
                    }
                    return(.fit_nmf(x = x, K = rank, background = background_signature))
                })
            }

            # select the best model based on mean cosine similarity
            models_objective <- rep(NA, length(results))
            for(j in 1:length(results)) {
                models_objective[j] <- results[[j]]$objective$goodness_fit
            }
            best_model <- which.max(models_objective)[1]
            alpha[[paste0(rank, "_signatures")]] <- results[[best_model]]$model$alpha
            beta[[paste0(rank, "_signatures")]] <- results[[best_model]]$model$beta
            unexplained_mutations[[paste0(rank, "_signatures")]] <- results[[best_model]]$model$unexplained_mutations
            cosine_similarity[[paste0(rank, "_signatures")]] <- results[[best_model]]$objective$cosine_similarities

            # compute stability among the inferred models
            if(nmf_runs == 1) {
                model_stability <- 1
            }
            else {
                model_stability <- rep(NA, (length(results)-1))
                cont <- 0
                for(j in 1:length(results)) {
                    if(j != best_model) {
                        cont <- cont + 1
                        model_stability[cont] <- .fit_stability(model1 = beta[[paste0(K[i], "_signatures")]], 
                            model2 = results[[j]]$model$beta)
                    }
                }
                model_stability <- mean(model_stability, na.rm = TRUE)
            }

            # save quality measures for the given solution
            curr_measures <- matrix(NA, nrow = 1, ncol = 4)
            rownames(curr_measures) <- paste0("K=",rank)
            colnames(curr_measures) <- c("Stability", "Mean Squared Error", "Explained Variance", "Predictions")
            curr_measures[1, "Stability"] <- model_stability
            curr_rank_predictions <- ((alpha[[paste0(rank, "_signatures")]] %*% beta[[paste0(rank, "_signatures")]]) 
                + unexplained_mutations[[paste0(rank, "_signatures")]])
            curr_rank_mse <- mean(((x - curr_rank_predictions)^2), na.rm = TRUE)
            curr_measures[1, "Mean Squared Error"] <- curr_rank_mse
            curr_pred_ss <- sum(((x - curr_rank_predictions)^2), na.rm = TRUE)
            curr_total_ss <- sum((x^2), na.rm = TRUE)
            curr_measures[1, "Explained Variance"] <- (1 - (curr_pred_ss / curr_total_ss))
            curr_measures[1, "Predictions"] <- results[[best_model]]$objective$goodness_fit
            measures <- rbind(measures, curr_measures)

            # rescale alpha to the original magnitude
            if (normalize_counts) {
                alpha[[paste0(K[i], "_signatures")]] <- alpha[[paste0(K[i], "_signatures")]] * (rowSums(x_not_normalized)/2500)
                unexplained_mutations[[paste0(K[i], "_signatures")]] <- unexplained_mutations[[paste0(K[i], "_signatures")]] * (rowSums(x_not_normalized)/2500)
            }
            rm(results)
            gc(verbose = FALSE)

        }

    } else {

        alpha <- list()
        alpha[["1_signatures"]] <- rank0_alpha
        beta <- list()
        beta[["1_signatures"]] <- rank0_beta
        unexplained_mutations <- list()
        unexplained_mutations[["1_signatures"]] <- rank0_unexplained_mutations
        cosine_similarity <- list()
        cosine_similarity[["1_signatures"]] <- rank0_cosine_similarity
        measures <- rank0_measures

    }

    # save results
    results <- list(alpha = alpha, beta = beta, unexplained_mutations = unexplained_mutations, cosine_similarity = cosine_similarity, measures = measures)
    rm(alpha)
    rm(beta)
    rm(unexplained_mutations)
    rm(cosine_similarity)
    rm(measures)

    # close parallel
    if (close_parallel) {
        stopCluster(parallel)
    }
    rm(close_parallel)
    rm(parallel)
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
#' sigs <- signaturesDecomposition(x = patients[seq_len(3),seq_len(2)],
#'                                 K = 3:4,
#'                                 background_signature = background[seq_len(2)],
#'                                 nmf_runs = 2,
#'                                 num_processes = 1)
#' set.seed(12345)
#' res <- signaturesCV(x = patients[seq_len(3),seq_len(2)],
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
#' @import lsa
#' @import nnls
#' @import RhpcBLASctl
#' @import parallel
#'
signaturesCV <- function( x, beta, normalize_counts = FALSE, cross_validation_entries = 0.01,
    cross_validation_iterations = 5, cross_validation_repetitions = 100, num_processes = Inf, verbose = TRUE ) {

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
    cv_estimates <- matrix(NA, nrow = cross_validation_repetitions, ncol = length(beta))
    rownames(cv_estimates) <- paste0("Repetition_", seq_len(cross_validation_repetitions))
    colnames(cv_estimates) <- seq_len(ncol(cv_estimates))

    # perform a total of cross_validation_repetitions repetitions of cross-validation
    valid_entries <- which(x > 0, arr.ind = TRUE)
    num_entries <- round((nrow(valid_entries) * cross_validation_entries))

    # performing the inference
    if (num_processes == 1) { # performing the inference sequentially

        for (cv_repetitions in seq_len(cross_validation_repetitions)) {

            if (verbose) {
                message("Performing repetition ", cv_repetitions, " out of ",
                    cross_validation_repetitions, "...", "\n")
            }

            # randomly set the cross-validation entries for the current iteration
            cv_entries <- valid_entries[sample(x = seq_len(nrow(valid_entries)), 
                size = num_entries, replace = FALSE), , drop = FALSE]

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

                    # set a percentage of cross_validation_entries entries to 0 to perform cross-validation
                    if (cv_iteration == 1) {
                        x_cv[cv_entries] <- 0
                    } else {
                        predicted_counts <- curr_results[["alpha"]] %*% curr_results[["beta"]]
                        x_cv[cv_entries] <- predicted_counts[cv_entries]
                    }

                    # perform the inference
                    curr_results <- .fit_model(x = x_cv, beta = beta[[num_signs]])
                    gc(verbose = FALSE)

                }

                # perform final fit of alpha
                curr_results <- signaturesAssignment(x = x_cv, beta = curr_results[["beta"]])

                # save an estimate of the current solution
                curr_predicted_counts <- ((curr_results[["alpha"]] %*% curr_results[["beta"]]) 
                    + curr_results[["unexplained_mutations"]])
                rm(curr_results)
                curr_true_considered_counts <- as.vector(x[cv_entries])
                curr_predicted_considered_counts <- as.vector(curr_predicted_counts[cv_entries])
                error <- (1 - as.numeric(cosine(curr_true_considered_counts,curr_predicted_considered_counts)))
                rm(curr_true_considered_counts)
                rm(curr_predicted_considered_counts)
                cv_estimates[cv_repetitions, paste0(k, "_signatures")] <- error

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
            cv_mean <- c(cv_mean, mean(cv_estimates[, i], na.rm = TRUE))
            cv_median <- c(cv_median, median(cv_estimates[, i], na.rm = TRUE))
        }
        cv_summary <- cbind(cv_mean, cv_median)
        rownames(cv_summary) <- colnames(cv_estimates)
        colnames(cv_summary) <- c("Mean cross-validation error", "Median cross-validation error")

    } else { # performing the inference in parallel

        parallel <- makeCluster(num_processes, outfile = "")
        close_parallel <- TRUE
        res_clusterEvalQ <- clusterEvalQ(parallel, library("glmnet", warn.conflicts = FALSE,
            quietly = TRUE, verbose = FALSE))
        res_clusterEvalQ <- clusterEvalQ(parallel, library("lsa", warn.conflicts = FALSE,
            quietly = TRUE, verbose = FALSE))
        res_clusterEvalQ <- clusterEvalQ(parallel, library("nnls", warn.conflicts = FALSE,
            quietly = TRUE, verbose = FALSE))
        clusterExport(parallel, varlist = c(".fit_model", "signaturesAssignment", "verbose", 
            "cross_validation_repetitions", "cross_validation_entries"), envir = environment())
        clusterExport(parallel, varlist = c("cross_validation_iterations", "valid_entries", "num_entries",
            "beta", "x"), envir = environment())
        clusterSetRNGStream(parallel, iseed = round(runif(1) * 1e+05))
        curr_results <- parLapply(parallel, seq_len(cross_validation_repetitions), function(cv_repetitions) {

            if (verbose) {
                message("Performing repetition ", cv_repetitions, " out of ",
                    cross_validation_repetitions, "...", "\n")
            }

            # randomly set the cross-validation entries for the current iteration
            cv_entries <- valid_entries[sample(x = seq_len(nrow(valid_entries)), 
                size = num_entries, replace = FALSE), , drop = FALSE]

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

                    # set a percentage of cross_validation_entries entries to 0 to perform cross-validation
                    if (cv_iteration == 1) {
                        x_cv[cv_entries] <- 0
                    } else {
                        predicted_counts <- curr_results[["alpha"]] %*% curr_results[["beta"]]
                        x_cv[cv_entries] <- predicted_counts[cv_entries]
                    }

                    # perform the inference
                    curr_results <- .fit_model(x = x_cv, beta = beta[[num_signs]])
                    gc(verbose = FALSE)
                }

                # perform final fit of alpha
                curr_results <- signaturesAssignment(x = x_cv, beta = curr_results[["beta"]])

                # save an estimate of the current solution
                curr_predicted_counts <- ((curr_results[["alpha"]] %*% curr_results[["beta"]]) 
                    + curr_results[["unexplained_mutations"]])
                rm(curr_results)
                curr_true_considered_counts <- as.vector(x[cv_entries])
                curr_predicted_considered_counts <- as.vector(curr_predicted_counts[cv_entries])
                error <- (1 - as.numeric(cosine(curr_true_considered_counts,curr_predicted_considered_counts)))
                rm(curr_true_considered_counts)
                rm(curr_predicted_considered_counts)
                cv_errors[num_signs] <- error

            }
            names(cv_errors) <- cv_names

            return(cv_errors)

        })

        # save the results from parallel computation
        colnames(cv_estimates) <- names(curr_results[[1]])
        for (par_res in seq_len(length(curr_results))) {
            cv_estimates[par_res, ] <- curr_results[[par_res]]
        }

        # compute mean and median values of estimated cross-validation error
        cv_mean <- NULL
        cv_median <- NULL
        for (i in seq_len(ncol(cv_estimates))) {
            cv_mean <- c(cv_mean, mean(cv_estimates[, i], na.rm = TRUE))
            cv_median <- c(cv_median, median(cv_estimates[, i], na.rm = TRUE))
        }
        cv_summary <- cbind(cv_mean, cv_median)
        rownames(cv_summary) <- colnames(cv_estimates)
        colnames(cv_summary) <- c("Mean cross-validation error", "Median cross-validation error")

    }

    # save results
    results <- list(estimates = cv_estimates, summary = cv_summary)
    rm(cv_estimates)
    rm(cv_summary)

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
