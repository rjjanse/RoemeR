#' Shrink coefficients
#'
#' @description
#' A helper function to apply shrinkage to a model's coefficients. Depending on the method used, different arguments
#' are required.
#'
#' @inheritParams pm_develop
#' @param output a list object containing the likelihood-ratio statistic and unpenalized coefficients.
#'
#' @return
#' @export
#'
#' @examples
.shrinkage <- function(formula, data, imputation, params, output){
    # Likelihood-based uniform shrinkage (LU)
    if(params[["shrink_method"]] == "lu") s <- .oc_lu(output)

    # Bootstrap-based uniform shrinkage (BU)
    if(params[["shrink_method"]] == "bu") s <- .oc_bu(formula, data, imputation, params)

    # Multiply coefficients with calibration slope for penalized coefficients
    penalized_coefs <- output[["coefs"]] * s

    # Create output list
    output <- list(penalized_coefs = penalized_coefs,
                   s = s)

    # Return penalized coefficients
    return(output)
}










