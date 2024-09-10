cslope <- function(data, formula, model){
    # Linear model
    if(model == "linear"){
        # Fit model
        fit <- stats::glm(formula = stats::as.formula(formula), family = "gaussian", data = data)

        # Get calibration intercept
        ci <- fit[["coefficients"]][["(Intercept)"]]

        # Get calibration slope
        cs <- fit[["coefficients"]][["lps"]]
    }

    # Logistic model
    if(model == "logistic"){
        fit <- stats::glm(formula = stats::as.formula(formula), family = "binomial", data = data)

        # Get calibration intercept
        ci <- fit[["coefficients"]][["(Intercept)"]]

        # Get calibration slope
        cs <- fit[["coefficients"]][["lps"]]
    }

    # Poisson model
    if(model == "poisson"){
        fit <- stats::glm(formula = stats::as.formula(formula), family = "poisson", data = data)

        # Get calibration intercept
        ci <- fit[["coefficients"]][["(Intercept)"]]

        # Get calibration slope
        cs <- fit[["coefficients"]][["lps"]]
    }

    # Accelerated failure time model
    if(model == "aft"){
        # Get calibration slope
        fit <- survival::survreg(formula = stats::as.formula(formula), x = TRUE, data = data)

        # Get calibration intercept
        ci <- fit[["coefficients"]][["(Intercept)"]]

        # Get calibration slope
        cs <- fit[["coefficients"]][["lps"]]
    }

    # Cox model
    if(model == "cox"){
        # Get calibration slope
        cs <- survival::coxph(formula = stats::as.formula(formula), x = TRUE, data = data)[["coefficients"]][["lps"]]
    }

    # Fine-Gray model
    if(model == "fine-gray"){
        # Get left-hand side of formula
        lhs <- stringr::str_extract(formula, ".*(?=~)")

        # Fit model
        fit <- survival::coxph(formula = survival::Surv(fgstart, fgstop, fgstatus) ~ lps, x = TRUE, weight = fgwt,
                               data = survival::finegray(stats::as.formula(paste0(lhs, "~ . ")), data = data))

        # Get calibration slope
        cs <- fit[["coefficients"]][["lps"]]
    }

    # Determine output based on presence of calibration intercept
    if(exists("ci")) output <- list("intercept" = ci, "slope" = cs) else output <- list("slope" = cs)

    # Return calibration slope
    return(output)
}
