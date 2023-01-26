# perform the inference for a given rank
.fit_nmf <- function( x, K, background = NULL ) {

    # generate the seed to start the inference
    seed <- .fit_seed(x = x, K = K)

    # perform the fit by elastic net with LASSO penalty
    model <- .fit_regularized( x = x, seed = seed, background = background )

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

# perform regularized fit by elastic net with LASSO penalty
.fit_regularized <- function( x, seed, background = NULL ) {

    # initialization
    alpha <- seed$alpha # exposures matrix
    beta <- seed$beta # signatures matrix
    if (!is.null(background)) {
        # if specified, the first row of beta is the background signature
        beta[1, ] <- background
    }
    N <- nrow(x) # N is the number of samples
    M <- ncol(x) # M is the number of features
    K <- ncol(alpha) # K is the number of latent variables to be fitted

    # iteratively fit alpha and beta by elastic net using LASSO penalty
    for (i in seq_len(10)) {
        # update alpha, beta is kept fixed
        for (j in seq_len(N)) {
            if (nrow(beta) > 1) {
                alpha[j, ] <- tryCatch({
                    res <- cv.glmnet(x = t(beta), y = as.vector(x[j, ]), 
                            type.measure = "mse", nfolds = 10, nlambda = 100, 
                            family = "gaussian", alpha = 1, lower.limits = 0, 
                            maxit = 1e+05)
                    res <- as.numeric(coef(res,s=res$lambda.min))
                    res <- (res[1]+res[-1])
                    is.invalid <- (res<0)
                    if(any(is.invalid)) {
                        res[is.invalid] <- 0
                    }
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
                    res <- (res[1]+res[-1])
                    is.invalid <- (res<0)
                    if(any(is.invalid)) {
                        res[is.invalid] <- 0
                    }
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
                # therefore it is not changed during the fit
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
    for (j in seq_len(N)) {
        if (nrow(beta) > 1) {
            alpha[j, ] <- tryCatch({
                res <- cv.glmnet(x = t(beta), y = as.vector(x[j, ]), 
                        type.measure = "mse", nfolds = 10, nlambda = 100, 
                        family = "gaussian", alpha = 1, lower.limits = 0, 
                        maxit = 1e+05)
                res <- as.numeric(coef(res,s=res$lambda.min))
                res <- (res[1]+res[-1])
                is.invalid <- (res<0)
                if(any(is.invalid)) {
                    res[is.invalid] <- 0
                }
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
                res <- (res[1]+res[-1])
                is.invalid <- (res<0)
                if(any(is.invalid)) {
                    res[is.invalid] <- 0
                }
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
    results <- list(alpha = alpha, beta = beta)
    rm(x)
    rm(seed)
    rm(alpha)
    rm(beta)
    gc(verbose = FALSE)

    # return the results
    return(results)

}

# compute goodness of fit estimates for a given model
.fit_objective <- function( x, model ) {

    # compute cosine similarities comparing observations and predictions
    predictions <- (model$alpha %*% model$beta)
    cosine_similarities <- rep(NA, nrow(x))
    names(cosine_similarities) <- rownames(x)
    for(i in seq_len(nrow(x))) {
        cosine_similarities[i] <- as.numeric(cosine(x[i,],predictions[i,]))
    }

    # compute goodness of fit for a given model
    goodness_fit <- median(cosine_similarities,na.rm=TRUE)

    # return the results
    results <- list(goodness_fit = goodness_fit, cosine_similarities = cosine_similarities)
    return(results)

}

# perform fit of a model solution by elastic net with LASSO penalty
.fit_model <- function( x, beta ) {

    # initialization
    alpha <- matrix(NA, nrow = nrow(x), ncol = nrow(beta))
    rownames(alpha) <- rownames(x)
    colnames(alpha) <- rownames(beta)
    N <- nrow(x) # N is the number of samples
    M <- ncol(x) # M is the number of features

    # iteratively fit alpha and beta by elastic net using LASSO penalty
    for (i in seq_len(10)) {
        # update alpha, beta is kept fixed
        for (j in seq_len(N)) {
            if (nrow(beta) > 1) {
                alpha[j, ] <- tryCatch({
                    res <- cv.glmnet(x = t(beta), y = as.vector(x[j, ]), 
                            type.measure = "mse", nfolds = 10, nlambda = 100, 
                            family = "gaussian", alpha = 1, lower.limits = 0, 
                            maxit = 1e+05)
                    res <- as.numeric(coef(res,s=res$lambda.min))
                    res <- (res[1]+res[-1])
                    is.invalid <- (res<0)
                    if(any(is.invalid)) {
                        res[is.invalid] <- 0
                    }
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
                    res <- (res[1]+res[-1])
                    is.invalid <- (res<0)
                    if(any(is.invalid)) {
                        res[is.invalid] <- 0
                    }
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
    gc(verbose = FALSE)

    # normalize beta and perform final fit of alpha
    beta <- (beta/rowSums(beta))
    is.invalid <- is.nan(rowSums(beta))
    if (any(is.invalid)) {
        beta[is.invalid, ] <- (1/ncol(beta))
    }
    for (j in seq_len(N)) {
        if (nrow(beta) > 1) {
            alpha[j, ] <- tryCatch({
                res <- cv.glmnet(x = t(beta), y = as.vector(x[j, ]), 
                        type.measure = "mse", nfolds = 10, nlambda = 100, 
                        family = "gaussian", alpha = 1, lower.limits = 0, 
                        maxit = 1e+05)
                res <- as.numeric(coef(res,s=res$lambda.min))
                res <- (res[1]+res[-1])
                is.invalid <- (res<0)
                if(any(is.invalid)) {
                    res[is.invalid] <- 0
                }
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
                res <- (res[1]+res[-1])
                is.invalid <- (res<0)
                if(any(is.invalid)) {
                    res[is.invalid] <- 0
                }
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
    results <- list(alpha = alpha, beta = beta)
    rm(x)
    rm(seed)
    rm(alpha)
    rm(beta)
    gc(verbose = FALSE)

    # return the results
    return(results)

}

# iteratively estimate alpha coefficients until a given level of cosine similarity is reached
.signatures_significance <- function( x, beta, cosine_thr = 0.95 ) {
    
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
            res <- (res[1]+res[-1])
            is.invalid <- (res<0)
            if(any(is.invalid)) {
                res[is.invalid] <- 0
            }
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
            res <- (res[1]+res[-1])
            is.invalid <- (res<0)
            if(any(is.invalid)) {
                res[is.invalid] <- 0
            }
            res
        }, error = function( e ) {
            res <- (sum(x[1,])/ncol(alpha))
            return(res)
        }, finally = {
            gc(verbose = FALSE)
        })
    }

    # iteratively include signatures into the fit until a given level of cosine similarity is reached
    sigs <- colnames(alpha[, which(alpha[1, ] > 0),
        drop = FALSE])[sort.int(alpha[, which(alpha[1, ] > 0)], decreasing = TRUE, index.return = TRUE)$ix]
    alpha <- matrix(0, nrow = 1, ncol = nrow(beta))
    rownames(alpha) <- rownames(x)
    colnames(alpha) <- rownames(beta)
    for (i in seq_len(length(sigs))) {

        # consider the current set of signatures and perform fit of alpha
        curr_alpha <- matrix(NA, nrow = 1, ncol = i)
        curr_beta <- beta[sigs[seq_len(i)], , drop = FALSE]
        if (nrow(curr_beta) > 1) {
            curr_alpha[1, ] <- tryCatch({
                res <- cv.glmnet(x = t(curr_beta), y = as.vector(x[1, ]), 
                        type.measure = "mse", nfolds = 10, nlambda = 100, 
                        family = "gaussian", alpha = 1, lower.limits = 0, 
                        maxit = 1e+05)
                res <- as.numeric(coef(res,s=res$lambda.min))
                res <- (res[1]+res[-1])
                is.invalid <- (res<0)
                if(any(is.invalid)) {
                    res[is.invalid] <- 0
                }
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
                res <- res[-length(res)]
                res <- (res[1]+res[-1])
                is.invalid <- (res<0)
                if(any(is.invalid)) {
                    res[is.invalid] <- 0
                }
                res
            }, error = function( e ) {
                res <- (sum(x[1,])/ncol(curr_alpha))
                return(res)
            }, finally = {
                gc(verbose = FALSE)
            })
        }

        # estimate goodness of fit
        curr_predicted_counts <- curr_alpha %*% curr_beta
        curr_goodness_fit <- as.numeric(cosine(as.numeric(x),as.numeric(curr_predicted_counts)))
        if (curr_goodness_fit > cosine_thr) {
            break
        }

    }
    alpha[, rownames(curr_beta)] <- curr_alpha

    # return the estimated signatures exposure
    return(alpha)

}
