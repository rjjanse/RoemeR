#' BU shrinkage factor
#'
#' @description
#' A helper function to calculate the shrinkage factor using bootstrap-based uniform shrinkage (BU). This function
#' is called by [RoemeR::.optimism_correction()].
#'
#' @inheritParams pm_develop
#' @param imputation an integer indicating the \emph{nth} imputation to be used from `data`.
#'
#' @return
#' @export
#'
#' @examples
.oc_bu <- function(formula, data, imputation, params){
    # Change formula to calibration slope formula
    formula_cs <- stringr::str_replace(formula, "~.*", "~lps")

    # Create temporary data
    dat_tmp <- dplyr::filter(data, .imp == imputation)

    # Bootstrapping
    straps <- lapply(1:params[["bootstraps"]], \(x){
        # Set seed
        set.seed(x)

        # Create bootstrap sample
        dat_boot <- dat_tmp[sample.int(nrow(dat_tmp), nrow(dat_tmp), replace = TRUE), ]

        # Parameters for bootstrap sample
        params_boot <- params

        # Change optimism correction to false to avoid infinite looping
        params_boot[["shrinkage"]] <- FALSE

        # Develop model on bootstrap sample
        fit_boot <- .develop_model(formula, dat_boot, imputation, params_boot, return_lps = TRUE)

        # Add linear predictors to data
        dat_tmp[["lps"]] <- fit_boot[["lps"]]

        # Calculate calibration slope
        cs <- cslope(dat_tmp, formula_cs, params[["model"]])[["slope"]]
    })

    # Calculate mean calibration slope (absolute value to not change predictor direction)
    cs <- abs(do.call("mean", straps))

    # Return calibration slope
    return(cs)
}
