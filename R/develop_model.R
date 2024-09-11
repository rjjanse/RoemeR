#' Develop a model on a single imputation
#'
#' @description
#' A helper function to develop a model on a single imputation. This function is called by [RoemeR::pm_develop()].
#'
#' @inheritParams pm_develop
#' @param imputation an integer indicating the \emph{nth} imputation to be used from `data`.
#' @param return_lps a logical indicating whether the linear predictors for each observation should be returned. These are used
#' by [RoemeR::.shrinkage()] when using bootstrap-based uniform shrinkage.
#'
#' @return a list object containing model information of a single imputation's model fit.
#' @export
#'
#' @examples
.develop_model <- function(formula, data, imputation, params, return_lps = FALSE){
    # Create temporary data
    dat_tmp <- dplyr::filter(data, .imp == imputation)

    # If shrinkage by LU, set return_lrs to TRUE
    if(params[["shrinkage"]] & params[["shrink_method"]] == "lu") return_lrs <- TRUE else return_lrs <- FALSE

    # Model fits ----
    # Linear model
    if(params[["model"]] == "linear") fit <- stats::glm(formula = stats::as.formula(formula), family = "gaussian",
                                                        data = dat_tmp)

    # Logistic model
    if(params[["model"]] == "logistic") fit <- stats::glm(formula = stats::as.formula(formula), family = "binomial",
                                                          data = dat_tmp)

    # Poisson model
    if(params[["model"]] == "poisson") fit <- stats::glm(formula = stats::as.formula(formula), family = "poisson",
                                                         data = dat_tmp)

    # Accelerated failure time model
    if(params[["model"]] == "aft") fit <- survival::survreg(formula = stats::as.formula(formula),
                                                            dist = params[["aft_distribution"]],
                                                            data = dat_tmp)

    # Cox model
    if(params[["model"]] == "cox") fit <- survival::coxph(formula = stats::as.formula(formula), x = TRUE, data = dat_tmp)

    # Fine-Gray model
    if(params[["model"]] == "fine-gray"){
        # Get left-hand side of formula
        lhs <- stringr::str_extract(formula, ".*(?=~)")

        # Get right-hand side of formula
        rhs <- stringr::str_extract(formula, "(?<=~).*")

        # Fit model
        fit <- survival::coxph(formula = stats::as.formula(paste0("survival::Surv(fgstart, fgstop, fgstatus) ~ ", rhs)),
                               x = TRUE, weight = fgwt, data = survival::finegray(stats::as.formula(paste0(lhs, "~ . ")),
                                                                                  data = dat_tmp))
    }

    # Get model information ----
    # Linear and logistic models
    if(params[["model"]] %in% c("linear", "logistic")){
        # Get coefficients
        coefs <- stats::coef(fit) %>%
            # Change to data frame
            as.data.frame() %>%
            # Transpose to make calculating mean later easier
            t() %>%
            # Rename intercept
            magrittr::set_colnames(c("a", colnames(.)[2:length(colnames(.))]))

        # Create list for output
        output <- list(coefs = coefs)
    }

    # AFT model
    if(params[["model"]] == "aft"){
        # Get coefficients
        coefs <- stats::coef(fit) %>%
            # Change to data frame
            as.data.frame() %>%
            # Transpose to make calculating mean later easier
            t() %>%
            # Rename intercept
            magrittr::set_colnames(c("a", colnames(.)[2:length(colnames(.))]))

        # Get scale parameter
        aft_scale <- fit[["scale"]]

        # Create list for output
        output <- list(coefs = coefs,
                       aft_scale = aft_scale)
    }

    # Poisson model
    if(params[["model"]] == "poisson"){
        # Get coefficients
        coefs <- stats::coef(fit) %>%
            # Change to data frame
            as.data.frame() %>%
            # Transpose to make calculating mean later easier
            t()

        # Create list for output
        output <- list(coefs = coefs)
    }

    # Cox and Fine-Gray model
    if(params[["model"]] %in% c("cox", "fine-gray")){
        # Get baseline hazard
        bh <- survival::basehaz(fit) %>%
            # Change to data frame
            as.data.frame() %>%
            # Remove everything after prediction horizon
            dplyr::filter(time <= params[["horizon"]]) %>%
            # Keep only last observation
            dplyr::slice_tail(n = 1L) %>%
            # Keep baseline hazard
            magrittr::extract2("hazard")

        # Get coefficients
        coefs <- stats::coef(fit) %>%
            # Change to data frame
            as.data.frame() %>%
            # Transpose to make calculating mean later easier
            t() %>%
            # Add baseline hazard
            dplyr::bind_cols(as.data.frame(bh), .)

        # Create list for output
        output <- list(coefs = coefs)

        # Export all baseline hazards if needed
        if(params[["all_hazards"]]){
            # All times
            times <- 1:params[["horizon"]]

            # Get baseline hazards for each timepoint
            bhs <- do.call("c", c(0, lapply(times[2:length(times)], \(x){
                # All baseline hazards
                survival::basehaz(fit) %>%
                    # Hazards before timepoint of interest
                    dplyr::filter(time <= x) %>%
                    # Keep hazard only
                    magrittr::extract2("hazard") %>%
                    # Keep last hazard only
                    dplyr::last()})))

            # Add baseline hazards to output
            output[["hazards"]] <- bhs
        }
    }

    # Add likelihood-ratio statistic if requested
    if(return_lrs) output[["lrs"]] <- lrtest(fit)[["Chisq"]][[2]]

    # Return linear predictors if requested
    if(return_lps) output[["lps"]] <- fit[["linear.predictors"]]

    # Shrink coefficients if requested
    if(params[["shrinkage"]]){
        # Rename unshrunk coefficients
        output[["unpenalized_coefs"]] <- output[["coefs"]]

        # Perform shrinkage
        shrink <- .shrinkage(formula, data, imputation, params, output)

        # Add shrunken coefficients to output
        output[["coefs"]] <- shrink[["penalized_coefs"]]

        # Add shrinkage factor to output
        output[["shrinkage_factor"]] <- shrink[["s"]]
    }

    # Return output
    return(output)
}
