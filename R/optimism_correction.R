.optimism_correction <- function(method, formula, data, imputation, params, coefficients){
    # Likelihood-based uniform shrinkage (LU)
    if(method == "lu") s <- .oc_lu(formula, data, imputation, params)

    # Bootstrap-based uniform shrinkage (BU)
    if(method == "bu") s <- .oc_bu(formula, data, imputation, params)

    # Multiply coefficients with calibration slope for penalized coefficients
    penalized_coefs <- coefficients * s

    # Return penalized coefficients
    return(penalized_coefs)
}










