.develop_model <- function(formula, data, imputation, params, return_lps = FALSE, return_lrs = FALSE){
    # Create temporary data
    dat_tmp <- dplyr::filter(data, .imp == imputation)

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

        # Export likelihood-ratio statistic if requested
        if(return_lrs) output[["lrs"]] <- summary(fit)[["logtest"]][["test"]]

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

    # Return linear predictors if requested
    if(return_lps) output[["lps"]] <- fit[["linear.predictors"]]

    # Penalize coefficients if requested
    if(params[["optimism_correction"]]){
        # Rename unpenalized coefficients
        output[["unpenalized_coefs"]] <- output[["coefs"]]

        # Add optimism adjusted coefficients to data
        output[["coefs"]] <- .optimism_correction(formula, data, imputation, params, output[["coefs"]])
    }

    # Return output
    return(output)
}
