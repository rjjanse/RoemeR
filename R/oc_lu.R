#' LU shrinkage factor
#'
#' @description
#' A helper function to calculate the shrinkage factor using likelihood-based uniform shrinkage (LU). This function
#' is called by [RoemeR::.shrinkage()].
#'
#' @inheritParams .shrinkage
#'
#' @return a single value representing the shrinkage factor
#' @export
#'
#' @examples
.oc_lu <- function(output){
    # Calculate shrinkage factor
    s <- (output[["lrs"]] - nrow(output[["coefs"]])) / output[["lrs"]]

    # Return shrinkage factor
    return(s)
}
