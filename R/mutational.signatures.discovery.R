#' Perform the assignment of K somatic mutational signatures provided as input to samples given a set of observed counts x.
#' This function can be used to estimate different types of mutational signatures such as: SBS (single base substitutions) and MNV (multi-nucleotide variant) 
#' (see Degasperi, Andrea, et al. "Substitution mutational signatures in whole-genome–sequenced cancers in the UK population." Science 376.6591 (2022): abl9283), 
#' CX (chromosomal instability) (see Drews, Ruben M., et al. "A pan-cancer compendium of chromosomal instability." Nature 606.7916 (2022): 976-983) and CN (copy number) 
#' signatures (see Steele, Christopher D., et al. "Signatures of copy number alterations in human cancer." Nature 606.7916 (2022): 984-991).
#'
#' @examples
#' data(background)
#' data(patients)
#' set.seed(12345)
#' beta <- signaturesDecomposition(x = patients[seq_len(3),], 
#'                                K = 3, 
#'                                background_signature = background, 
#'                                nmf_runs = 2, 
#'                                sparsify = FALSE, 
#'                                num_processes = 1)
#' set.seed(12345)
#' res <- signaturesAssignment(x = patients[seq_len(3),], beta = beta$beta[[1]], sparsify = FALSE)
#'
#' @title signaturesAssignment
#' @param x counts matrix for a set of n patients and m categories. These can be, e.g., SBS, MNV, CN or CN counts;
#' in the case of SBS it would be an n patients x 96 trinucleotides matrix.
#' @param beta matrix of the discovered signatures to be used for the assignment.
#' @param normalize_counts if true, the input counts matrix x is normalized such that the patients have the same number of mutation.
#' @param sparsify boolean; Shall I perform regularization using LASSO?
#' @param verbose boolean; Shall I print information messages?
#' @return A list with the discovered signatures. It includes 2 elements:
#'              alpha: matrix of the discovered exposure values.
#'              beta: matrix of the discovered signatures.
#' @export signaturesAssignment
#' @import nnls
#' @importFrom glmnet cv.glmnet
#'
signaturesAssignment <- function( x, beta, normalize_counts = FALSE, sparsify = TRUE, verbose = TRUE ) {
    
    # check input parameters
    x <- as.matrix(x)
    if(normalize_counts) {
        x_not_normalized <- x
        x <- (x/rowSums(x))*2500
    }

    if(verbose) {
        cat("Performing signatures assignment...","\n")
    }

    # initialize alpha with an empty matrix
    alpha <- array(NA,c(nrow(x),nrow(beta)))
    rownames(alpha) <- rownames(x)
    colnames(alpha) <- rownames(beta)
    
    # perform signatures assignment
    if(sparsify) {
        for(j in seq_len(nrow(alpha))) {
            curr_beta_values <- beta
            if(nrow(curr_beta_values)>1) {
                res <- cv.glmnet(t(curr_beta_values),as.vector(x[j,]),type.measure="mse",nfolds=10,nlambda=10,family="gaussian",lower.limits=0.00)
                alpha[j,] <- as.numeric(res$glmnet.fit$beta[,ncol(res$glmnet.fit$beta)])
            }
            else {
                res_inputs <- cbind(t(curr_beta_values),matrix(rep(0,ncol(curr_beta_values)),ncol=1))
                res <- cv.glmnet(res_inputs,as.vector(x[j,]),type.measure="mse",nfolds=10,nlambda=10,family="gaussian",lower.limits=0.00)
                res <- as.numeric(res$glmnet.fit$beta[,ncol(res$glmnet.fit$beta)])
                alpha[j,] <- res[-length(res)]
            }
        }
    }
    else {
        for(j in seq_len(nrow(alpha))) {
            alpha[j,] <- nnls(t(beta),as.vector(x[j,]))$x
        }
    }

    # rescale alpha to the original magnitude
    if(normalize_counts) {
        alpha <- alpha * (rowSums(x_not_normalized)/2500)
    }

    # return the assigned signatures
    results <- list(alpha=alpha,beta=beta)
    return(results)
    
}

#' Perform signatures discovery and rank estimation for a range of K somatic mutational signatures given a set of observed counts x.
#' This function can be used to estimate different types of mutational signatures such as: SBS (single base substitutions) and MNV (multi-nucleotide variant) 
#' (see Degasperi, Andrea, et al. "Substitution mutational signatures in whole-genome–sequenced cancers in the UK population." Science 376.6591 (2022): abl9283), 
#' CX (chromosomal instability) (see Drews, Ruben M., et al. "A pan-cancer compendium of chromosomal instability." Nature 606.7916 (2022): 976-983) and CN (copy number) 
#' signatures (see Steele, Christopher D., et al. "Signatures of copy number alterations in human cancer." Nature 606.7916 (2022): 984-991).
#'
#' @examples
#' data(background)
#' data(patients)
#' set.seed(12345)
#' res <- signaturesDecomposition(x = patients[seq_len(3),], 
#'                               K = 3:4, 
#'                               background_signature = background, 
#'                               nmf_runs = 2, 
#'                               sparsify = FALSE, 
#'                               num_processes = 1)
#'
#' @title signaturesDecomposition
#' @param x counts matrix for a set of n patients and m categories. These can be, e.g., SBS, MNV, CN or CN counts;
#' in the case of SBS it would be an n patients x 96 trinucleotides matrix.
#' @param K either one value or a range of numeric values (each of them greater than 0) indicating the number of signatures to be considered.
#' @param background_signature background signature to be used.
#' @param normalize_counts if true, the input counts matrix x is normalized such that the patients have the same number of mutation.
#' @param nmf_runs number of iteration (minimum 1) of NMF to be performed for a robust estimation of beta.
#' @param sparsify boolean; Shall I perform regularization using LASSO?
#' @param num_processes Number of processes to be used during parallel execution. To execute in single process mode, 
#' this parameter needs to be set to either NA or NULL.
#' @param verbose boolean; Shall I print information messages?
#' @return A list with the discovered signatures and related rank measures. It includes 3 elements:
#'              alpha: list of matrices of the discovered exposure values for each possible rank in the range K.
#'              beta: list of matrices of the discovered signatures for each possible rank in the range K.
#'              measures: a data.frame containing the quality measures for each possible rank in the range K.
#' @export signaturesDecomposition
#' @import NMF
#' @import nnls
#' @import parallel
#' @importFrom cluster silhouette
#' @importFrom glmnet cv.glmnet
#'
signaturesDecomposition <- function( x, K, background_signature = NULL, normalize_counts = FALSE, nmf_runs = 50, sparsify = TRUE, num_processes = Inf, verbose = TRUE ) {
    
    # check input parameters
    x <- as.matrix(x)
    if(any(colSums(x)==0)) {
        invalid_cols <- as.numeric(which(colSums(x)==0))
        for(inv_cols in invalid_cols) {
            x[sample(seq_len(length(x[,inv_cols])),size=1),inv_cols] <- 1e-05
        }
    }
    K <- sort(unique(K))
    if(min(K)<1) {
        warning("The minimum value of K can be 1, removing invalid values...")
        K <- K[K>=1]
        if(length(K)==0) {
            K <- 1
            warning("No valid value of K is provided: setting K equals to 1...")
        }
    }
    if(normalize_counts) {
        x_not_normalized <- x
        x <- (x/rowSums(x))*2500
    }
    if(nmf_runs<1) {
        warning("The minimum value of nmf_runs can be 1...")
        nmf_runs <- 1
    }

    # setting up parallel execution
    pbackend <- NULL
    close_parallel <- FALSE
    if(is.na(num_processes) || is.null(num_processes)) {
        num_processes <- 1
    }
    else if(num_processes==Inf) {
        cores <- as.integer((detectCores()-1))
        num_processes <- min(cores,nmf_runs)
    }
    else {
        num_processes <- min(num_processes,nmf_runs)
    }
    if(num_processes==1) {
        pbackend <- "seq"
    }
    else {
        pbackend <- makeCluster(num_processes)
        res_clusterEvalQ <- clusterEvalQ(pbackend,library("NMF",warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE))
        res_clusterEvalQ <- clusterEvalQ(pbackend,library("nnls",warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE))
        res_clusterEvalQ <- clusterEvalQ(pbackend,library("glmnet",warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE))
        clusterExport(pbackend,varlist=c(".nmf_seed",".nmf_nnls",".nmf_lasso"),envir=environment())
        clusterSetRNGStream(pbackend,iseed=round(runif(1)*100000))
        close_parallel <- TRUE
        rm(res_clusterEvalQ)
        gc(verbose=FALSE)
    }

    if(verbose) {
        cat("Performing signatures discovery and rank estimation...","\n")
        if(num_processes>1) {
            cat("Executing",num_processes,"processes in parallel...","\n")
        }
    }

    # handle special case when solution with only background_signature is also considered
    rank0 <- FALSE
    if(K[1]==1&&!is.null(background_signature)) {
        if(verbose) {
            cat("Performing inference for K=1...","\n")
        }
        rank0 <- TRUE
        K <- K[-1]
        rank0_K <- 1
        rank0_beta <- matrix(background_signature,nrow=1)
        rownames(rank0_beta) <- "Background"
        colnames(rank0_beta) <- colnames(x)
        rank0_alpha <- signaturesAssignment(x=x,beta=rank0_beta,normalize_counts=FALSE,sparsify=sparsify,verbose=FALSE)$alpha
        rank0_rss <- rss(rank0_alpha%*%rank0_beta,x)
        rank0_evar <- evar(rank0_alpha%*%rank0_beta,x)
        rank0_sparseness_alpha <- sparseness(as.vector(rank0_alpha))
        rank0_sparseness_beta <- sparseness(as.vector(rank0_beta))
        if(normalize_counts) {
            rank0_alpha <- rank0_alpha * (rowSums(x_not_normalized)/2500)
        }
        if(nmf_runs>1) {
            rank0_measures <- matrix(c(rank0_K,rank0_rss,rank0_evar,NA,NA,rank0_sparseness_alpha,rank0_sparseness_beta,1,1,1),nrow=1)
            colnames(rank0_measures) <- c("Rank","Residual_Sum_Squares","Explained_Variance","Silhouette_Alpha","Silhouette_Beta","Sparseness_Alpha","Sparseness_Beta","Cophenetic_Coefficient","Dispersion_Coefficient","Silhouette_Consensus")
        }
        else {
            rank0_measures <- matrix(c(rank0_K,rank0_rss,rank0_evar,NA,NA,rank0_sparseness_alpha,rank0_sparseness_beta),nrow=1)
            colnames(rank0_measures) <- c("Rank","Residual_Sum_Squares","Explained_Variance","Silhouette_Alpha","Silhouette_Beta","Sparseness_Alpha","Sparseness_Beta")
        }
    }

    # perform signatures discovery and rank estimation
    if(length(K)>0) {

        alpha <- list()
        beta <- list()
        measures <- NULL
        if(rank0) {
            alpha[["1_signatures"]] <- rank0_alpha
            beta[["1_signatures"]] <- rank0_beta
            measures <- rank0_measures
        }
        
        for(i in seq_len(length(K))) {

            if(verbose) {
                cat(paste0("Performing inference for K=",K[i],"..."),"\n")
            }

            # perform the inference for current K
            results <- NULL
            if(sparsify) {
                while(is.null(results)) {
                    results <- tryCatch({
                        results <- nmf(x=x,rank=K[i],method=.nmf_lasso,background=background_signature,seed=.nmf_seed,rng=round(runif(1)*10000),nrun=nmf_runs,.pbackend=pbackend)
                        gc(verbose=FALSE)
                        results
                    }, error = function(e) {
                        cat(paste0("An error has occurred: ",e$message),"\n")
                        gc(verbose=FALSE)
                        NULL
                    }, finally = {
                        gc(verbose=FALSE)
                    })
                }
            }
            else {
                while(is.null(results)) {
                    results <- tryCatch({
                        results <- nmf(x=x,rank=K[i],method=.nmf_nnls,background=background_signature,seed=.nmf_seed,rng=round(runif(1)*10000),nrun=nmf_runs,.pbackend=pbackend)
                        gc(verbose=FALSE)
                        results
                    }, error = function(e) {
                        cat(paste0("An error has occurred: ",e$message),"\n")
                        gc(verbose=FALSE)
                        NULL
                    }, finally = {
                        gc(verbose=FALSE)
                    })
                }
            }

            alpha[[paste0(K[i],"_signatures")]] <- basis(results)
            beta[[paste0(K[i],"_signatures")]] <- coef(results)

            # rescale alpha to the original magnitude
            if(normalize_counts) {
                alpha[[paste0(K[i],"_signatures")]] <- alpha[[paste0(K[i],"_signatures")]] * (rowSums(x_not_normalized)/2500)
            }

            # compute and save quality measures
            curr_rss <- rss(results,x)
            curr_evar <- evar(results,x)
            curr_silhouette_alpha <- tryCatch({
                mean(silhouette(results,what="features")[,"sil_width"])
            }, error = function(e) {
                NA
            })
            curr_silhouette_beta <- tryCatch({
                mean(silhouette(results,what="samples")[,"sil_width"])
            }, error = function(e) {
                NA
            })
            curr_sparseness_alpha <- sparseness(as.vector(alpha[[paste0(K[i],"_signatures")]]))
            curr_sparseness_beta <- sparseness(as.vector(beta[[paste0(K[i],"_signatures")]]))
            if(nmf_runs>1) {
                curr_silhouette_consensus <- tryCatch({
                    mean(silhouette(results,what="chc")[,"sil_width"])
                }, error = function(e) {
                    NA
                })
                curr_measures <- matrix(c(K[i],curr_rss,curr_evar,curr_silhouette_alpha,curr_silhouette_beta,curr_sparseness_alpha,curr_silhouette_beta,cophcor(results),dispersion(results),curr_silhouette_consensus),nrow=1)
                colnames(curr_measures) <- c("Rank","Residual_Sum_Squares","Explained_Variance","Silhouette_Alpha","Silhouette_Beta","Sparseness_Alpha","Sparseness_Beta","Cophenetic_Coefficient","Dispersion_Coefficient","Silhouette_Consensus")
            }
            else {
                curr_measures <- matrix(c(K[i],curr_rss,curr_evar,curr_silhouette_alpha,curr_silhouette_beta,curr_sparseness_alpha,curr_silhouette_beta),nrow=1)
                colnames(curr_measures) <- c("Rank","Residual_Sum_Squares","Explained_Variance","Silhouette_Alpha","Silhouette_Beta","Sparseness_Alpha","Sparseness_Beta")
            }
            measures <- rbind(measures,curr_measures)
            rm(results)
            rm(curr_measures)
            gc(verbose=FALSE)

        }
        rownames(measures) <- seq_len(nrow(measures))

    }
    else {

        alpha <- list()
        beta <- list()
        alpha[["1_signatures"]] <- rank0_alpha
        beta[["1_signatures"]] <- rank0_beta
        measures <- rank0_measures
        rownames(measures) <- seq_len(nrow(measures))

    }

    # save results
    results <- list(alpha=alpha,beta=beta,measures=measures)
    rm(alpha)
    rm(beta)
    rm(measures)
    gc(verbose=FALSE)
    
    # close pbackend
    if(close_parallel) {
        stopCluster(pbackend)
    }
    rm(close_parallel)
    rm(pbackend)
    gc(verbose=FALSE)

    # return the discovered signatures
    return(results)
    
}

#' Perform the assessment of different signaturesDecomposition solutions by cross validation for K (beta, as estimated by signaturesDecomposition) somatic 
#' mutational signatures given a set of observations x and discovered signatures beta.
#'
#' @examples
#' data(background)
#' data(patients)
#' set.seed(12345)
#' sigs <- signaturesDecomposition(x = patients[seq_len(3),], 
#'                                K = 3:4, 
#'                                background_signature = background, 
#'                                nmf_runs = 2, 
#'                                sparsify = FALSE, 
#'                                num_processes = 1)
#' set.seed(12345)
#' res <- signaturesCV(x = patients[seq_len(3),], 
#'                    beta = sigs$beta, 
#'                    cross_validation_iterations = 2, 
#'                    cross_validation_repetitions = 2, 
#'                    num_processes = 1)
#'
#' @title signaturesCV
#' @param x counts matrix for a set of n patients and m categories. These can be, e.g., trinucleotides counts for n patients and 96 trinucleotides.
#' @param beta a set of inferred signatures as returned by signaturesDecomposition function.
#' @param normalize_counts if true, the input counts matrix x is normalized such that the patients have the same number of mutation.
#' @param cross_validation_entries Percentage of cells in the counts matrix to be replaced by 0s during cross validation.
#' @param cross_validation_iterations For each configuration, the first time the signatures are fitted form a matrix with a 
#' percentage of values replaced by 0s. This may result in poor fit/results. Then, we perform predictions of these entries and replace them with 
#' such predicted values. This parameter is the number of restarts to be performed to improve this estimate and obtain more stable solutions.
#' @param cross_validation_repetitions Number of time cross-validation should be repeated. Higher values result in better estimate, but are computationally 
#' more expensive.
#' @param num_processes Number of processes to be used during parallel execution. To execute in single process mode, 
#' this parameter needs to be set to either NA or NULL.
#' @param verbose boolean; Shall I print information messages?
#' @return A list of 2 elements: estimates and summary. Here, cv_estimates reports the mean squared error for each configuration of performed 
#' cross validation; rank_estimates reports mean and median values for each value of K.
#' @export signaturesCV
#' @import nnls
#' @import parallel
#' @importFrom glmnet cv.glmnet
#'
signaturesCV <- function( x, beta, normalize_counts = FALSE, cross_validation_entries = 0.01, cross_validation_iterations = 5, cross_validation_repetitions = 100, num_processes = Inf, verbose = TRUE ) {
    
    # check input parameters
    x <- as.matrix(x)
    if(any(colSums(x)==0)) {
        invalid_cols <- as.numeric(which(colSums(x)==0))
        for(inv_cols in invalid_cols) {
            x[sample(seq_len(length(x[,inv_cols])),size=1),inv_cols] <- 1e-05
        }
    }
    if(normalize_counts) {
        x <- (x/rowSums(x))*2500
    }

    # setting up parallel execution
    parallel <- NULL
    close_parallel <- FALSE
    if(is.na(num_processes) || is.null(num_processes)) {
        num_processes <- 1
    }
    else if(num_processes==Inf) {
        cores <- as.integer((detectCores()-1))
        num_processes <- min(cores,cross_validation_repetitions)
    }
    else {
        num_processes <- min(num_processes,cross_validation_repetitions)
    }

    if(verbose) {
        cat("Estimating the optimal number of signatures with a total of",cross_validation_repetitions,"cross validation repetitions...","\n")
        if(num_processes>1) {
            cat("Executing",num_processes,"processes in parallel...","\n")
        }
    }

    # structure to save the results
    cv_estimates <- array(NA,c(cross_validation_repetitions,length(beta)))
    rownames(cv_estimates) <- paste0("Repetition_",seq_len(cross_validation_repetitions))
    colnames(cv_estimates) <- seq_len(ncol(cv_estimates))

    # perform a total of cross_validation_repetitions repetitions of cross validation
    valid_entries <- which(x>0,arr.ind=TRUE)

    # performing inference
    if(num_processes==1) { # performing inference sequentially

        for(cv_repetitions in seq_len(cross_validation_repetitions)) {

            if(verbose) {
                cat(paste0("Performing repetition ",cv_repetitions," out of ",cross_validation_repetitions,"..."),"\n")
            }

            # randomly set the cross validation entries for the current iteration
            cv_entries <- valid_entries[sample(seq_len(nrow(valid_entries)),size=round(nrow(valid_entries)*cross_validation_entries),replace=FALSE),]
            
            # consider all the possible values of K
            cont <- 0
            for(num_signs in seq_len(length(beta))) {

                k <- nrow(beta[[num_signs]])
                if(cv_repetitions==1) {
                    colnames(cv_estimates)[num_signs] <- paste0(k,"_signatures")
                }

                if(verbose) {
                    cat(paste0("Performing estimation for K=",k,"..."),"\n")
                }
                
                # repeat the estimation for a number of cross_validation_iterations
                x_cv <- x
                for(cv_iteration in seq_len(cross_validation_iterations)) {

                    if(verbose) {
                        cat(paste0("Performing cross validation iteration ",cv_iteration," out of ",cross_validation_iterations,"..."),"\n")
                    }

                    # set a percentage of cross_validation_entries entries to 0 in order to perform cross validation
                    if(cv_iteration==1) {
                        x_cv[cv_entries] <- 0
                    }
                    else {
                        predicted_counts <- curr_results[["alpha"]]%*%curr_results[["beta"]]
                        x_cv[cv_entries] <- predicted_counts[cv_entries]
                    }

                    # columns in x_cv with all 0s are not allowed
                    if(any(colSums(x_cv)==0)) {
                        invalid_cols <- as.numeric(which(colSums(x_cv)==0))
                        for(inv_cols in invalid_cols) {
                            x_cv[sample(seq_len(length(x_cv[,inv_cols])),size=1),inv_cols] <- 1e-05
                        }
                    }

                    # perform the inference
                    curr_results <- .nmf_fit(x=x_cv,beta=beta[[num_signs]])
                    gc(verbose=FALSE)

                }

                # save an estimate of the current solution
                curr_predicted_counts <- curr_results[["alpha"]]%*%curr_results[["beta"]]
                rm(curr_results)
                curr_true_considered_counts <- as.vector(x[cv_entries])
                curr_predicted_considered_counts <- as.vector(curr_predicted_counts[cv_entries])
                error <- mean((curr_true_considered_counts-curr_predicted_considered_counts)^2)
                rm(curr_true_considered_counts)
                rm(curr_predicted_considered_counts)
                cv_estimates[cv_repetitions,paste0(k,"_signatures")] <- error
                rm(error)
                gc(verbose=FALSE)

                if(verbose) {
                    cont <- cont + 1
                    cat("Progress",paste0(round((cont/length(beta))*100,digits=3),"%..."),"\n")
                }

            }

        }

        # compute mean and median values of estimated cross validation error
        cv_mean <- NULL
        cv_median <- NULL
        for(i in seq_len(ncol(cv_estimates))) {
            cv_mean <- c(cv_mean,mean(cv_estimates[,i]))
            cv_median <- c(cv_median,median(cv_estimates[,i]))
        }
        cv_summary <- cbind(cv_mean,cv_median)
        rownames(cv_summary) <- colnames(cv_estimates)
        colnames(cv_summary) <- c("Mean cross-validation error","Median cross-validation error")

    }
    else { # performing inference in parallel

        parallel <- makeCluster(num_processes,outfile="")
        close_parallel <- TRUE
        res_clusterEvalQ <- clusterEvalQ(parallel,library("nnls",warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE))
        res_clusterEvalQ <- clusterEvalQ(parallel,library("glmnet",warn.conflicts=FALSE,quietly=TRUE,verbose=FALSE))
        clusterExport(parallel,varlist=c(".nmf_fit"),envir=environment())
        clusterExport(parallel,varlist=c("verbose","cross_validation_repetitions","cross_validation_entries"),envir=environment())
        clusterExport(parallel,varlist=c("cross_validation_iterations","valid_entries","beta","x"),envir=environment())
        clusterSetRNGStream(parallel,iseed=round(runif(1)*100000))
        rm(res_clusterEvalQ)
        gc(verbose=FALSE)
        curr_results <- parLapply(parallel,seq_len(cross_validation_repetitions),function(cv_repetitions) {

            if(verbose) {
                cat(paste0("Performing repetition ",cv_repetitions," out of ",cross_validation_repetitions,"..."),"\n")
            }

            # randomly set the cross validation entries for the current iteration
            cv_entries <- valid_entries[sample(seq_len(nrow(valid_entries)),size=round(nrow(valid_entries)*cross_validation_entries),replace=FALSE),]
            
            # consider all the possible values of K
            cv_errors <- rep(NA,length(beta))
            cv_names <- NULL
            for(num_signs in seq_len(length(beta))) {

                if(cv_repetitions==1) {
                    cv_names <- c(cv_names,paste0(nrow(beta[[num_signs]]),"_signatures"))
                }

                # repeat the estimation for a number of cross_validation_iterations
                x_cv <- x
                for(cv_iteration in seq_len(cross_validation_iterations)) {

                    # set a percentage of cross_validation_entries entries to 0 in order to perform cross validation
                    if(cv_iteration==1) {
                        x_cv[cv_entries] <- 0
                    }
                    else {
                        predicted_counts <- curr_results[["alpha"]]%*%curr_results[["beta"]]
                        x_cv[cv_entries] <- predicted_counts[cv_entries]
                    }

                    # columns in x_cv with all 0s are not allowed
                    if(any(colSums(x_cv)==0)) {
                        invalid_cols <- as.numeric(which(colSums(x_cv)==0))
                        for(inv_cols in invalid_cols) {
                            x_cv[sample(seq_len(length(x_cv[,inv_cols])),size=1),inv_cols] <- 1e-05
                        }
                    }

                    # perform the inference
                    curr_results <- .nmf_fit(x=x_cv,beta=beta[[num_signs]])
                    gc(verbose=FALSE)

                }

                # save an estimate of the current solution
                curr_predicted_counts <- curr_results[["alpha"]]%*%curr_results[["beta"]]
                rm(curr_results)
                curr_true_considered_counts <- as.vector(x[cv_entries])
                curr_predicted_considered_counts <- as.vector(curr_predicted_counts[cv_entries])
                error <- mean((curr_true_considered_counts-curr_predicted_considered_counts)^2)
                rm(curr_true_considered_counts)
                rm(curr_predicted_considered_counts)
                cv_errors[num_signs] <- error
                rm(error)
                gc(verbose=FALSE)

            }
            names(cv_errors) <- cv_names

            return(cv_errors)

        })

        # save the results from parallel computation
        colnames(cv_estimates) <- names(curr_results[[1]])
        for(par_res in seq_len(length(curr_results))) {
            cv_estimates[par_res,] <- curr_results[[par_res]]
        }
        rm(curr_results)
        gc(verbose=FALSE)

        # compute mean and median values of estimated cross validation error
        cv_mean <- NULL
        cv_median <- NULL
        for(i in seq_len(ncol(cv_estimates))) {
            cv_mean <- c(cv_mean,mean(cv_estimates[,i]))
            cv_median <- c(cv_median,median(cv_estimates[,i]))
        }
        cv_summary <- cbind(cv_mean,cv_median)
        rownames(cv_summary) <- colnames(cv_estimates)
        colnames(cv_summary) <- c("Mean cross-validation error","Median cross-validation error")

    }

    # save results
    results <- list(estimates=cv_estimates,summary=cv_summary)
    rm(cv_estimates)
    rm(cv_summary)
    gc(verbose=FALSE)
    
    # close parallel
    if(close_parallel) {
        stopCluster(parallel)
    }
    rm(close_parallel)
    rm(parallel)
    gc(verbose=FALSE)

    # return results of cross validation
    return(results)
    
}

# initialize alpha and beta for .nmf_nnls or .nmf_lasso functions
.nmf_seed <- function( model, target ) {

    # initialize alpha with an empty matrix
    alpha <- array(NA,c(nrow(target),nbasis(model)))
    rownames(alpha) <- rownames(target)
    colnames(alpha) <- paste0("S",seq_len(ncol(alpha)))

    # randomly initialize beta
    beta <- matrix(rnbinom(nbasis(model)*ncol(target),prob=0.10,size=1),nrow=nbasis(model),ncol=ncol(target))
    beta <- (beta/rowSums(beta)) # beta rows (i.e., signatures) must sum to 1
    if(any(is.nan(rowSums(beta)))) {
        beta[which(is.nan(rowSums(beta))),] <- 0
    }
    rownames(beta) <- colnames(alpha)
    colnames(beta) <- colnames(target)

    # update results
    basis(model) <- alpha
    coef(model) <- beta
    rm(alpha)
    rm(beta)
    gc(verbose=FALSE)

    # return updated model
    return(model)

}

# perform NMF by Non-negative least squares
.nmf_nnls <- function( x, seed, background = NULL ) {

    # initialization
    alpha <- basis(seed) # exposures matrix
    beta <- coef(seed) # signatures matrix
    if(!is.null(background)) { # if specified, first row of beta is background
        beta[1,] <- background
    }
    n <- nrow(x) # n is the number of observations in x, i.e., the patients
    J <- ncol(x) # J is the number of trinucleotides, i.e., 96 categories
    K <- nrow(beta) # K is the number of signatures to be fitted

    # iteratively fit alpha and beta by Non-negative least squares (nnls)
    for(i in seq_len(20)) {

        # update alpha, beta is kept fixed
        for(j in seq_len(n)) {
            alpha[j,] <- nnls(t(beta),as.vector(x[j,]))$x
        }

        # update beta, alpha is kept fixed
        for(k in seq_len(J)) {
            if(!is.null(background)) {
                # the first signature represents the background model, thus it is not changed during the fit
                beta[2:K,k] <- nnls(alpha[,2:K,drop=FALSE],as.vector(x[,k]-(alpha[,1]*beta[1,k])))$x
            }
            else {
                beta[,k] <- nnls(alpha,as.vector(x[,k]))$x
            }
        }

    }

    # normalize beta and perform final fit of alpha
    beta <- (beta/rowSums(beta))
    if(any(is.nan(rowSums(beta)))) {
        beta[which(is.nan(rowSums(beta))),] <- 0
    }
    for(j in seq_len(n)) {
        alpha[j,] <- nnls(t(beta),as.vector(x[j,]))$x
    }
    if(!is.null(background)) {
        colnames(alpha) <- c("Background",paste0("S",seq_len((ncol(alpha)-1))))
    }
    else {
        colnames(alpha) <- paste0("S",seq_len(ncol(alpha)))
    }
    rownames(beta) <- colnames(alpha)

    # update results
    basis(seed) <- alpha
    coef(seed) <- beta
    rm(x)
    rm(alpha)
    rm(beta)
    gc(verbose=FALSE)

    # return updated seed
    return(seed)

}

# perform NMF by Non-negative least squares and Non-Negative Lasso
.nmf_lasso <- function( x, seed, background = NULL ) {

    # initialization
    alpha <- basis(seed) # exposures matrix
    beta <- coef(seed) # signatures matrix
    if(!is.null(background)) { # if specified, first row of beta is background
        beta[1,] <- background
    }
    n <- nrow(x) # n is the number of observations in x, i.e., the patients
    J <- ncol(x) # J is the number of trinucleotides, i.e., 96 categories
    K <- nrow(beta) # K is the number of signatures to be fitted

    # STEP 1: we get close to a good solution by Non-negative least squares

    # iteratively fit alpha and beta by Non-negative least squares (nnls)
    for(i in seq_len(20)) {

        # update alpha, beta is kept fixed
        for(j in seq_len(n)) {
            alpha[j,] <- nnls(t(beta),as.vector(x[j,]))$x
        }

        # update beta, alpha is kept fixed
        for(k in seq_len(J)) {
            if(!is.null(background)) {
                # the first signature represents the background model, thus it is not changed during the fit
                beta[2:K,k] <- nnls(alpha[,2:K,drop=FALSE],as.vector(x[,k]-(alpha[,1]*beta[1,k])))$x
            }
            else {
                beta[,k] <- nnls(alpha,as.vector(x[,k]))$x
            }
        }

    }

    # normalize beta and perform final fit of alpha
    beta <- (beta/rowSums(beta))
    if(any(is.nan(rowSums(beta)))) {
        beta[which(is.nan(rowSums(beta))),] <- 0
    }
    for(j in seq_len(n)) {
        alpha[j,] <- nnls(t(beta),as.vector(x[j,]))$x
    }
    if(!is.null(background)) {
        colnames(alpha) <- c("Background",paste0("S",seq_len((ncol(alpha)-1))))
    }
    else {
        colnames(alpha) <- paste0("S",seq_len(ncol(alpha)))
    }
    rownames(beta) <- colnames(alpha)

    # STEP 2: we finalize the inference by Non-Negative Lasso
    
    # now regularize solutions by Non-Negative Lasso
    for(i in seq_len(10)) {

        # update alpha, beta is kept fixed
        for(j in seq_len(n)) {
            curr_beta_values <- beta
            if(nrow(curr_beta_values)>1) {
                res <- cv.glmnet(t(curr_beta_values),as.vector(x[j,]),type.measure="mse",nfolds=10,nlambda=10,family="gaussian",lower.limits=0.00)
                alpha[j,] <- as.numeric(res$glmnet.fit$beta[,ncol(res$glmnet.fit$beta)])
            }
            else {
                res_inputs <- cbind(t(curr_beta_values),matrix(rep(0,ncol(curr_beta_values)),ncol=1))
                res <- cv.glmnet(res_inputs,as.vector(x[j,]),type.measure="mse",nfolds=10,nlambda=10,family="gaussian",lower.limits=0.00)
                res <- as.numeric(res$glmnet.fit$beta[,ncol(res$glmnet.fit$beta)])
                alpha[j,] <- res[-length(res)]
            }
        }

        # update beta, alpha is kept fixed
        for(k in seq_len(J)) {
            if(!is.null(background)) {
                # the first signature represents the background model, thus it is not changed during the fit
                curr_alpha_values <- alpha[,2:K,drop=FALSE]
                if(ncol(curr_alpha_values)>1) {
                    res <- cv.glmnet(x=curr_alpha_values,y=as.vector(x[,k]-(alpha[,1]*beta[1,k])),type.measure="mse",nfolds=10,nlambda=10,family="gaussian",lower.limits=0.00)
                    beta[2:K,k] <- as.numeric(res$glmnet.fit$beta[,ncol(res$glmnet.fit$beta)])
                }
                else {
                    res_inputs <- cbind(curr_alpha_values,rep(0,nrow(curr_alpha_values)))
                    res <- cv.glmnet(x=res_inputs,y=as.vector(x[,k]-(alpha[,1]*beta[1,k])),type.measure="mse",nfolds=10,nlambda=10,family="gaussian",lower.limits=0.00)
                    res <- as.numeric(res$glmnet.fit$beta[,ncol(res$glmnet.fit$beta)])
                    beta[2:K,k] <- res[-length(res)]
                }
            }
            else {
                curr_alpha_values <- alpha
                if(ncol(curr_alpha_values)>1) {
                    res <- cv.glmnet(x=curr_alpha_values,y=as.vector(x[,k]),type.measure="mse",nfolds=10,nlambda=10,family="gaussian",lower.limits=0.00)
                    beta[,k] <- as.numeric(res$glmnet.fit$beta[,ncol(res$glmnet.fit$beta)])
                }
                else {
                    res_inputs <- cbind(curr_alpha_values,rep(0,nrow(curr_alpha_values)))
                    res <- cv.glmnet(x=res_inputs,y=as.vector(x[,k]),type.measure="mse",nfolds=10,nlambda=10,family="gaussian",lower.limits=0.00)
                    res <- as.numeric(res$glmnet.fit$beta[,ncol(res$glmnet.fit$beta)])
                    beta[,k] <- res[-length(res)]
                }
            }
        }

    }

    # normalize beta and perform final fit of alpha
    beta <- (beta/rowSums(beta))
    if(any(is.nan(rowSums(beta)))) {
        beta[which(is.nan(rowSums(beta))),] <- 0
    }
    for(j in seq_len(n)) {
        curr_beta_values <- beta
        if(nrow(curr_beta_values)>1) {
            res <- cv.glmnet(t(curr_beta_values),as.vector(x[j,]),type.measure="mse",nfolds=10,nlambda=10,family="gaussian",lower.limits=0.00)
            alpha[j,] <- as.numeric(res$glmnet.fit$beta[,ncol(res$glmnet.fit$beta)])
        }
        else {
            res_inputs <- cbind(t(curr_beta_values),matrix(rep(0,ncol(curr_beta_values)),ncol=1))
            res <- cv.glmnet(res_inputs,as.vector(x[j,]),type.measure="mse",nfolds=10,nlambda=10,family="gaussian",lower.limits=0.00)
            res <- as.numeric(res$glmnet.fit$beta[,ncol(res$glmnet.fit$beta)])
            alpha[j,] <- res[-length(res)]
        }
    }
    if(!is.null(background)) {
        colnames(alpha) <- c("Background",paste0("S",seq_len((ncol(alpha)-1))))
    }
    else {
        colnames(alpha) <- paste0("S",seq_len(ncol(alpha)))
    }
    rownames(beta) <- colnames(alpha)

    # update results
    basis(seed) <- alpha
    coef(seed) <- beta
    rm(x)
    rm(alpha)
    rm(beta)
    gc(verbose=FALSE)

    # return updated seed
    return(seed)

}

# perform fit of NMF solution by Non-negative least squares
.nmf_fit <- function( x, beta ) {

    # initialization
    alpha <- array(NA,c(nrow(x),nrow(beta)))
    rownames(alpha) <- rownames(x)
    colnames(alpha) <- rownames(beta)
    beta <- (beta/rowSums(beta))
    if(any(is.nan(rowSums(beta)))) {
        beta[which(is.nan(rowSums(beta))),] <- 0
    }
    n <- nrow(x) # n is the number of observations in x, i.e., the patients
    J <- ncol(x) # J is the number of trinucleotides, i.e., 96 categories
    
    # iteratively fit alpha and beta by Non-negative least squares (nnls)
    for(i in seq_len(20)) {

        # update alpha, beta is kept fixed
        for(j in seq_len(n)) {
            alpha[j,] <- nnls(t(beta),as.vector(x[j,]))$x
        }

        # update beta, alpha is kept fixed
        for(k in seq_len(J)) {
            beta[,k] <- nnls(alpha,as.vector(x[,k]))$x
        }

    }

    # normalize beta and perform final fit of alpha
    beta <- (beta/rowSums(beta))
    if(any(is.nan(rowSums(beta)))) {
        beta[which(is.nan(rowSums(beta))),] <- 0
    }
    for(j in seq_len(n)) {
        alpha[j,] <- nnls(t(beta),as.vector(x[j,]))$x
    }

    # update results
    results <- list(alpha=alpha,beta=beta)
    rm(x)
    rm(alpha)
    rm(beta)
    gc(verbose=FALSE)

    # return the results
    return(results)

}
