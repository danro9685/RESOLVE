#' Perform signatures based clustering given the signatures assignments matrix alpha.
#'
#' @examples
#' data(sbs_assignments)
#' set.seed(12345)
#' norm_alpha = (sbs_assignments$alpha / rowSums(sbs_assignments$alpha))
#' sbs_clustering = signaturesClustering(alpha = norm_alpha, num_clusters = 1:3, num_processes = 1, verbose = FALSE)
#'
#' @title signaturesClustering
#' @param alpha Signatures assignments matrix.
#' @param num_clusters Range of number of clusters to be considered.
#' @param num_processes Number of processes to be used during parallel execution. To execute in single process mode,
#' this parameter needs to be set to either NA or NULL.
#' @param verbose Boolean. Shall I print information messages?
#' @return A list a clusters assignments for each number within the num_clusters range.
#' @export signaturesSignificance
#' @import cluster
#' @import lsa
#' @import parallel
#'
signaturesClustering <- function( alpha, num_clusters = 1:10, num_processes = Inf, verbose = TRUE ) {

    # setting up parallel execution
    parallel <- NULL
    close_parallel <- FALSE
    if (is.na(num_processes) || is.null(num_processes)) {
        num_processes <- 1
    } else if (num_processes == Inf) {
        cores <- as.integer((detectCores() - 1))
        num_processes <- min(cores, length(num_clusters))
    } else {
        num_processes <- min(num_processes, length(num_clusters))
    }

    if (verbose) {
        message("Estimating signatures-based clusters...", "\n")
        if (num_processes > 1) {
            message("Executing ", num_processes, " processes in parallel...", "\n")
        }
    }

    # computing the signatures-based distance
    if (verbose) {
        message("Computing a signatures-based distance...", "\n")
    }
    cos_similarity <- cosine(t(alpha))
    invalid <- which(is.na(cos_similarity))
    if(length(invalid)>0) {
        cos_similarity[invalid] <- 0
    }
    signatures_distance <- (1-cos_similarity)
    invalid <- which(signatures_distance<0)
    if(length(invalid)>0) {
        signatures_distance[invalid] <- 0
    }
    invalid <- which(signatures_distance>1)
    if(length(invalid)>0) {
        signatures_distance[invalid] <- 1
    }
    sig_dist <- as.dist(signatures_distance)

    # performing the inference
    if (num_processes == 1) { # performing the inference sequentially
        sig_clustering <- list()
        num_clust_pos <- 0
        for(num_clust in num_clusters) {
            num_clust_pos <- num_clust_pos + 1
            sig_clustering[[num_clust_pos]] <- pam( x = sig_dist,
                                                  k = num_clust,
                                                  diss = TRUE,
                                                  medoids = "random",
                                                  nstart = 100,
                                                  cluster.only = TRUE,
                                                  do.swap = TRUE,
                                                  keep.diss = FALSE,
                                                  keep.data = FALSE)
        }
    } else { # performing the inference in parallel
        parallel <- makeCluster(num_processes, outfile = "")
        close_parallel <- TRUE
        res_clusterEvalQ <- clusterEvalQ(parallel, library("cluster", warn.conflicts = FALSE,
            quietly = TRUE, verbose = FALSE))
        clusterExport(parallel, varlist = "sig_dist", envir = environment())
        clusterSetRNGStream(parallel, iseed = round(runif(1) * 1e+05))
        sig_clustering <- parLapply(parallel, num_clusters, function(num_clust) {
            if (verbose) {
                message("Performing clustering with k=", num_clust, "...", "\n")
            }
            res_clustering <- pam( x = sig_dist,
                                  k = num_clust,
                                  diss = TRUE,
                                  medoids = "random",
                                  nstart = 100,
                                  cluster.only = TRUE,
                                  do.swap = TRUE,
                                  keep.diss = FALSE,
                                  keep.data = FALSE)
            return(res_clustering)
        })
    }

    # close parallel
    if (close_parallel) {
        stopCluster(parallel)
    }
    rm(close_parallel)
    rm(parallel)
    gc(verbose = FALSE)

    # save the results
    results <- sig_clustering
    names(results) <- paste0("k=",num_clusters)

    # return the estimeted signatures-based clusters
    return(results)

}
