# perform the inference of K mutational signatures
.fit_nmf <- function( x, K, background ) {

    # perform the fit of the model by elastic net with LASSO penalty
    model <- .fit_regularized( x = x, seed = .fit_seed(x = x, K = K), 
        background = background )

    # evaluate goodness of fit of the inferred model
    objective <- .fit_objective( x = x, model = model )

    # save the results
    results <- list( model = model, objective = objective)

    # return the results
    return(results)

}

# initialize alpha and beta for the inference
.fit_seed <- function( x, K ) {

    # initialize alpha as an empty matrix
    alpha <- matrix(NA, nrow = nrow(x), ncol = K)
    rownames(alpha) <- rownames(x)
    colnames(alpha) <- paste0("S", seq_len(K))

    # randomly initialize beta with a negative binomial distribution
    beta <- matrix(rnbinom(n = (K * ncol(x)), size = 5, prob = 0.1), 
        nrow = K, ncol = ncol(x))
    beta <- (beta/rowSums(beta)) # rows of beta must sum to 1
    is.invalid <- is.nan(rowSums(beta))
    if (any(is.invalid)) {
        beta[is.invalid, ] <- (1/ncol(beta))
    }
    rownames(beta) <- colnames(alpha)
    colnames(beta) <- colnames(x)

    # update the seed
    seed <- list(alpha = alpha, beta = beta)

    # return the seed
    return(seed)

}

# perform the inference with regularization by elastic net with LASSO penalty
.fit_regularized <- function( x, seed, background ) {

    # initialization
    alpha <- seed$alpha # exposures matrix
    beta <- seed$beta # signatures matrix
    if (!is.null(background)) {
        # if specified, the first row of beta is the background signature
        beta[1, ] <- background
    }
    N <- nrow(x) # N is the number of samples
    M <- ncol(x) # M is the number of features
    K <- ncol(alpha) # K is the number of signatures to be fitted

    # fit alpha and beta by non-negative least squares to get close to a good solution
    for (i in seq_len(20)) {
        # update alpha, beta is kept fixed
        for (j in seq_len(N)) {
            alpha[j, ] <- nnls(A = t(beta), b = as.vector(x[j, ]))$x
        }
        # update beta, alpha is kept fixed
        for (k in seq_len(M)) {
            if (!is.null(background)) {
                # when provided, the first signature represents the background model, 
                # therefore, it is not changed during the fit
                beta[2:K, k] <- nnls(A = alpha[, 2:K, drop = FALSE], b = as.vector(x[, k]))$x
            } else {
                beta[, k] <- nnls(A = alpha, b = as.vector(x[, k]))$x
            }
        }
    }

    # normalize beta
    beta <- (beta/rowSums(beta))
    is.invalid <- is.nan(rowSums(beta))
    if (any(is.invalid)) {
        beta[is.invalid, ] <- (1/ncol(beta))
    }

    # refine the fit of alpha and beta by elastic net using LASSO penalty
    for (i in seq_len(3)) {
        # update alpha, beta is kept fixed
        for (j in seq_len(N)) {
            if (nrow(beta) > 1) {
                alpha[j, ] <- tryCatch({
                    res <- cv.glmnet(x = t(beta), y = as.vector(x[j, ]), 
                            type.measure = "mse", nfolds = 10, nlambda = 100, 
                            family = "gaussian", alpha = 1, lower.limits = 0, 
                            maxit = 1e+05)
                    res <- as.numeric(coef(res,s=res$lambda.min))
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
        # update beta, alpha is kept fixed
        for (k in seq_len(M)) {
            if (!is.null(background)) {
                # when provided, the first signature represents the background model, 
                # therefore, it is not changed during the fit
                if (ncol(alpha[, 2:K, drop = FALSE]) > 1) {
                    beta[2:K, k] <- tryCatch({
                        res <- cv.glmnet(x = alpha[, 2:K, drop = FALSE], y = as.vector(x[, k]), 
                            type.measure = "mse", nfolds = 10, nlambda = 100, 
                            family = "gaussian", alpha = 1, lower.limits = 0, 
                            upper.limits = 1, maxit = 1e+05)
                        res <- as.numeric(coef(res,s=res$lambda.min))
                        res <- res[-1]
                        res
                    }, error = function( e ) {
                        res <- 0
                        return(res)
                    }, finally = {
                        gc(verbose = FALSE)
                    })
                } else {
                    fit_inputs <- cbind(alpha[, 2:K, drop = FALSE], rep(0, nrow(alpha[, 2:K, drop = FALSE])))
                    beta[2:K, k] <- tryCatch({
                        res <- cv.glmnet(x = fit_inputs, y = as.vector(x[, k]), 
                            type.measure = "mse", nfolds = 10, nlambda = 100, 
                            family = "gaussian", alpha = 1, lower.limits = 0, 
                            upper.limits = 1, maxit = 1e+05)
                        res <- as.numeric(coef(res,s=res$lambda.min))
                        res <- res[-length(res)]
                        res <- res[-1]
                        res
                    }, error = function( e ) {
                        res <- 0
                        return(res)
                    }, finally = {
                        gc(verbose = FALSE)
                    })
                }
            } else {
                if (ncol(alpha) > 1) {
                    beta[, k] <- tryCatch({
                        res <- cv.glmnet(x = alpha, y = as.vector(x[, k]), 
                            type.measure = "mse", nfolds = 10, nlambda = 100, 
                            family = "gaussian", alpha = 1, lower.limits = 0, 
                            upper.limits = 1, maxit = 1e+05)
                        res <- as.numeric(coef(res,s=res$lambda.min))
                        res <- res[-1]
                        res
                    }, error = function( e ) {
                        res <- 0
                        return(res)
                    }, finally = {
                        gc(verbose = FALSE)
                    })
                } else {
                    fit_inputs <- cbind(alpha, rep(0, nrow(alpha)))
                    beta[, k] <- tryCatch({
                        res <- cv.glmnet(x = fit_inputs, y = as.vector(x[, k]), 
                            type.measure = "mse", nfolds = 10, nlambda = 100, 
                            family = "gaussian", alpha = 1, lower.limits = 0, 
                            upper.limits = 1, maxit = 1e+05)
                        res <- as.numeric(coef(res,s=res$lambda.min))
                        res <- res[-length(res)]
                        res <- res[-1]
                        res
                    }, error = function( e ) {
                        res <- 0
                        return(res)
                    }, finally = {
                        gc(verbose = FALSE)
                    })
                }
            }
        }
    }
    gc(verbose = FALSE)

    # normalize beta and perform final fit of alpha
    beta <- (beta/rowSums(beta))
    is.invalid <- is.nan(rowSums(beta))
    if (any(is.invalid)) {
        beta[is.invalid, ] <- (1/ncol(beta))
    }
    unexplained_mutations <- rep(0, nrow(x))
    names(unexplained_mutations) <- rownames(x)
    for (j in seq_len(N)) {
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

    # save the results
    results <- list(alpha = alpha, beta = beta, unexplained_mutations = unexplained_mutations)
    rm(x)
    rm(seed)
    rm(alpha)
    rm(beta)
    gc(verbose = FALSE)

    # return the results
    return(results)

}

# compute an estimate of goodness of fit for a given model
.fit_objective <- function( x, model ) {

    # compute cosine similarities comparing observations and predictions
    predictions <- ((model$alpha %*% model$beta) + model$unexplained_mutations)
    cosine_similarities <- rep(NA, nrow(x))
    names(cosine_similarities) <- rownames(x)
    for(i in seq_len(nrow(x))) {
        cosine_similarities[i] <- as.numeric(cosine(x[i,],predictions[i,]))
    }

    # compute goodness of fit for a given model
    goodness_fit <- mean(cosine_similarities,na.rm=TRUE)

    # return the results
    results <- list(goodness_fit = goodness_fit, cosine_similarities = cosine_similarities)
    return(results)

}

# estimate mean cosine similarity between two models
.fit_stability <- function( model1, model2 ) {

    # consider the signatures from the first model
    mean_cos_sig <- rep(NA, nrow(model1))
    search_model <- seq_len(nrow(model2))
    for(i in seq_len(nrow(model1))) {
        # compare each signatures to the ones from the second model
        cos_sig <- rep(NA, length(search_model))
        pos <- 0
        for(j in search_model) {
            # compute cosine similarity between the two signatures
            pos <- pos + 1
            cos_sig[pos] <- as.numeric(cosine(model1[i,],model2[j,]))
        }
        best_match <- which.max(cos_sig)[1]
        search_model <- search_model[-best_match]
        mean_cos_sig[i] <- cos_sig[best_match]
    }
    mean_cos_sig <- mean(mean_cos_sig, na.rm = TRUE)

    # return the mean cosine similarity between the two models
    return(mean_cos_sig)

}

# perform the inference given an estimate of beta
.fit_model <- function( x, beta ) {

    # initialization
    alpha <- matrix(NA, nrow = nrow(x), ncol = nrow(beta))
    rownames(alpha) <- rownames(x)
    colnames(alpha) <- rownames(beta)
    N <- nrow(x) # N is the number of samples
    M <- ncol(x) # M is the number of features

    # fit alpha and beta by non-negative least squares to get close to a good solution
    for (i in seq_len(10)) {
        # update alpha, beta is kept fixed
        for (j in seq_len(N)) {
            alpha[j, ] <- nnls(A = t(beta), b = as.vector(x[j, ]))$x
        }
        # update beta, alpha is kept fixed
        for (k in seq_len(M)) {
            beta[, k] <- nnls(A = alpha, b = as.vector(x[, k]))$x
        }
    }

    # normalize beta and perform final fit of alpha
    beta <- (beta/rowSums(beta))
    is.invalid <- is.nan(rowSums(beta))
    if (any(is.invalid)) {
        beta[is.invalid, ] <- (1/ncol(beta))
    }
    for (j in seq_len(N)) {
        alpha[j, ] <- nnls(A = t(beta), b = as.vector(x[j, ]))$x
    }

    # save the results
    results <- list(alpha = alpha, beta = beta)
    rm(x)
    rm(alpha)
    rm(beta)
    gc(verbose = FALSE)

    # return the results
    return(results)

}

# iteratively estimate alpha until a given level of cosine similarity is reached
.signatures_significance <- function( x, beta, cosine_thr ) {
    
    # initialization
    alpha <- matrix(NA, nrow = 1, ncol = nrow(beta))
    rownames(alpha) <- rownames(x)
    colnames(alpha) <- rownames(beta)
    if (nrow(beta) > 1) {
        alpha[1, ] <- tryCatch({
            res <- cv.glmnet(x = t(beta), y = as.vector(x[1, ]), 
                    type.measure = "mse", nfolds = 10, nlambda = 100, 
                    family = "gaussian", alpha = 1, lower.limits = 0, 
                    maxit = 1e+05)
            res <- as.numeric(coef(res,s=res$lambda.min))
            res <- res[-1]
            res
        }, error = function( e ) {
            res <- (sum(x[1,])/ncol(alpha))
            return(res)
        }, finally = {
            gc(verbose = FALSE)
        })
    } else {
        fit_inputs <- cbind(t(beta), rep(0, ncol(beta)))
        alpha[1, ] <- tryCatch({
            res <- cv.glmnet(x = fit_inputs, y = as.vector(x[1, ]), 
                    type.measure = "mse", nfolds = 10, nlambda = 100, 
                    family = "gaussian", alpha = 1, lower.limits = 0, 
                    maxit = 1e+05)
            res <- as.numeric(coef(res,s=res$lambda.min))
            res <- res[-length(res)]
            res <- res[-1]
            res
        }, error = function( e ) {
            res <- (sum(x[1,])/ncol(alpha))
            return(res)
        }, finally = {
            gc(verbose = FALSE)
        })
    }
    gc(verbose = FALSE)

    # iteratively include signatures into the model until a given level of cosine similarity is reached
    sigs <- names(sort(alpha[, which(alpha[1, ] > 0)], decreasing = TRUE))
    alpha <- matrix(0, nrow = 1, ncol = nrow(beta))
    rownames(alpha) <- rownames(x)
    colnames(alpha) <- rownames(beta)
    for (i in seq_len(length(sigs))) {
        # consider the current set of signatures and perform fit of alpha
        curr_alpha <- matrix(NA, nrow = 1, ncol = i)
        curr_beta <- beta[sigs[seq_len(i)], , drop = FALSE]
        unexplained_mutations <- 0
        if (nrow(curr_beta) > 1) {
            curr_alpha[1, ] <- tryCatch({
                res <- cv.glmnet(x = t(curr_beta), y = as.vector(x[1, ]), 
                        type.measure = "mse", nfolds = 10, nlambda = 100, 
                        family = "gaussian", alpha = 1, lower.limits = 0, 
                        maxit = 1e+05)
                res <- as.numeric(coef(res,s=res$lambda.min))
                unexplained_mutations <- res[1]
                res <- res[-1]
                res
            }, error = function( e ) {
                res <- (sum(x[1,])/ncol(curr_alpha))
                return(res)
            }, finally = {
                gc(verbose = FALSE)
            })
        } else {
            fit_inputs <- cbind(t(curr_beta), rep(0, ncol(curr_beta)))
            curr_alpha[1, ] <- tryCatch({
                res <- cv.glmnet(x = fit_inputs, y = as.vector(x[1, ]), 
                        type.measure = "mse", nfolds = 10, nlambda = 100, 
                        family = "gaussian", alpha = 1, lower.limits = 0, 
                        maxit = 1e+05)
                res <- as.numeric(coef(res,s=res$lambda.min))
                unexplained_mutations <- res[1]
                res <- res[-length(res)]
                res <- res[-1]
                res
            }, error = function( e ) {
                res <- (sum(x[1,])/ncol(curr_alpha))
                return(res)
            }, finally = {
                gc(verbose = FALSE)
            })
        }
        # estimate goodness of fit
        curr_predicted_counts <- ((curr_alpha %*% curr_beta) + unexplained_mutations)
        curr_goodness_fit <- as.numeric(cosine(as.numeric(x),as.numeric(curr_predicted_counts)))
        if (curr_goodness_fit > cosine_thr) {
            break
        }
    }
    alpha[, rownames(curr_beta)] <- curr_alpha
    gc(verbose = FALSE)

    # return the estimated signatures exposure
    return(alpha)

}
