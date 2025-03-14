#' Perform the estimation of mutational signatures associated to prognosis.
#'
#' @examples
#' data(association_survival)
#' set.seed(12345)
#' clinical_data = association_survival$clinical_data
#' normalized_alpha = association_survival$normalized_alpha
#' prognosis_associations = associationPrognosis(clinical_data = clinical_data, signatures = normalized_alpha)
#'
#' @title associationPrognosis
#' @param clinical_data Matrix with two columns proving survival times (first column) and status (second column).
#' @param signatures Matrix with the estimated exposures to signatures (columns) for each patient (rows).
#' @return A vector of coefficients resulting from regularized Cox regression.
#' @export associationPrognosis
#' @import glmnet
#' @import survival
#'
associationPrognosis <- function( clinical_data, signatures ) {
    
    # perform the analysis
    x <- as.matrix(signatures)
    y <- Surv(as.numeric(clinical_data[,1]),as.numeric(clinical_data[,2]))
    cv.fit <- cv.glmnet(x = x, y = y, family = "cox", nfolds = 10, maxit = 1e+05)
    coeff_fit <- coef(cv.fit,s=cv.fit$lambda.min)[colnames(signatures), 1, drop = TRUE]
    coeff_values <- as.numeric(coeff_fit)
    names(coeff_values) <- names(coeff_fit)
    prognosis_associations <- coeff_values

    # return the associations
    return(prognosis_associations)

}
