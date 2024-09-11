#' Develop a prediction model
#'
#' @description
#' `pm_develop()` develops a prediction model according to the specified regression method.
#'
#' @param formula a string representing a formula object. The outcome should be on the left of a tilde (~)
#' and the predictors should be on the right. The outcome should adhere to the specification of the used
#' regression model (e.g. using Surv() for a Cox regression model).
#' @param data a matrix, data frame, or tibble containing the outcome variable(s) and predictors.
#' @param params a list of parameters for fitting the model:
#' * `model`: the regression model to be fitted. Should be one of 'linear', 'logistic', 'poisson', 'aft'
#' (accelerated failure time), 'cox', or 'fine-gray'.
#' * `imputation_column`: a string indicating the column containing imputation numbers. This defaults to .imp as
#' created by \pkg{mice}.
#' * `horizon`: prediction horizon for survival models (aft, cox, and fine-gray).
#' * `aft_distribution`: a string indicating the distribution with which the accelerated failure time model should be fitted.
#' * `all_hazards`: a logical indicating whether all hazards should be extracted and returned, instead of only the
#' baseline hazard.
#' * `shrinkage`: a logical indicating whether coefficients should be shrunk (i.e. optimism corrected/penalized) (see Details).
#' * `shrink_method`: a string indicating the method for shrinkage. This can be 'BU' (bootstrap-based uniform shrinkage)
#' or 'LU' (likelihood-based uniform shrinkage).
#' * `bootstraps`: the number of bootstraps to be used for shrinkage by bootstrap-based uniform shrinkage.
#'
#' @details
#' All parameters, except for `imputation_column` are case-insensitive.
#'
#' Shrinkage is performed according to the specified methods as in Van Calster \emph{et al.} 2020. For bootstrap-based
#' uniform shrinkage, by default 500 bootstraps are used. Likelihood-based uniform shrinkage uses \link[lmtest]{lrtest} to
#' derive the likelihood-ratio statistic and degrees of freedom. Note, as also mentioned in the study, improved performance
#' is not guaranteed.
#'
#' @return A list object of class 'pm' containing the model type, model coefficients, supplied parameters, and if shrinkage
#' was applied, the unpenalized coefficients.
#'.
#' @export
#' @examples
#' # Load survival and dplyr packages
#' library(survival); library(dplyr)
#'
#' # Change status of data to 0 and 1
#' dat <- mutate(lung, status = status - 1)
#'
#' # Set parameters for Cox model with outcome at 1 year
#' params <- list(model = "cox",
#'                horizon = 365.25)
#'
#' # Create model formula
#' form <- "Surv(time, status) ~ age + sex"
#'
#' # Develop model
#' fit <- pm_develop(formula, dat, params)
#'
#' @references
#' Van Calster B, van Smeden M, De Cock B, Steyerberg EW. Regression shrinkage methods for clinical prediction models
#' do not guarantee improved performance: Simulation study. Stat Methods Med Res. 2020 Nov;29(11):3166-3178. doi:
#' \href{https://doi.org/10.1177/0962280220921415}{10.1177/0962280220921415}.
#'
pm_develop <- function(formula, data, params){
    # Params argument checks ----
    # Params should be of type list
    if(!is.list(params)){
        stop("params should be of type list")
    }

    # Model should be defined
    if(is.null(params[["model"]])){
        stop("Define 'model' in params")
    }

    # Standardization of params to lower case ----
    for(i in 1:length(params)){
        # Imputation column should not be made lower case as it refers to a column in the data
        params[[i]] <- ifelse(is.character(params[[i]]) & names(params[i]) != "imputation_column",
                              tolower(params[[i]]), params[[i]])
    }

    # Further argument checks ----
    # The model should be one of the available models
    if(!(params[["model"]] %in% c("linear", "logistic", "poisson", "aft", "cox", "fine-gray"))){
        stop("'model' should be one of 'linear', 'logistic', 'poisson', 'aft', 'cox', or 'fine-gray'")
    }

    # If the model is a survival model, horizon should be defined
    if(params[["model"]] %in% c("fine-gray", "cox", "aft") & is.null(params[["horizon"]])){
        stop("Define 'horizon' in params for survival models")
    }

    # If the model is an accelerated failure time model, distribution should be defined
    if(params[["model"]] == "aft" & is.null(params[["aft_distribution"]])){
        stop("Define 'aft_distribution' in params for accelerated failure time models")
    }

    # Data should be either a matrix or data frame
    if(!is.data.frame(data) & !is.matrix(data)){
        stop("data should be of type matrix, data frame, or tibble")
    }

    # Set argument values ----
    # If data is matrix, change to data frame
    if(is.matrix(data)) data <- as.data.frame(data)

    # Export all baseline hazards for Cox & Fine-Gray model
    if(is.null(params[["all_hazards"]])) params[["all_hazards"]] <- FALSE

    # If no imputation column is given, set to no imputation of .imp is also absent
    if(is.null(params[["imputation_column"]])){
        # If .imp is not present, assume data was not imputed and add singular imputation column
        if(is.null(data[[".imp"]])) data[[".imp"]] <- 1
    }

    # If an imputation column is given, set to .imp
    if(!is.null(params[["imputation_column"]])) data[[".imp"]] <- data[[params[["imputation_column"]]]]

    # If a Cox model is requested, add package specification to Surv()
    if(params[["model"]] == "cox" & !stringr::str_detect(formula, "survival::Surv")){
        formula <- stringr::str_replace(formula, "^Surv", "survival::Surv")
    }

    # If shrinkage is not specified, set to false
    if(is.null(params[["shrinkage"]])) params[["shrinkage"]] <- FALSE

    # If shrinkage is specified, but method is not, set to BU
    if(params[["shrinkage"]] & is.null(params[["shrink_method"]])) params[["shrink_method"]] <- "bu"

    # If optimism correction is requested but no bootstraps are requested, set to 500
    if(params[["shrinkage"]] & params[["shrink_method"]] == "bu" & is.null(params[["bootstraps"]])){
        params[["bootstraps"]] <- 500
    }

    # Get final model ----
    # Fit models for each imputation
    models <- lapply(unique(data[[".imp"]]), \(x) .develop_model(formula, data, x, params))

    # Get coefficients
    coefs <- do.call("rbind", lapply(models, \(x) x[["coefs"]])) %>%
        # Change to data frame
        as.data.frame() %>%
        # Final values are mean of each column
        dplyr::mutate(dplyr::across(tidyselect::everything(), mean)) %>%
        # Keep only first row
        dplyr::slice(1L) %>%
        # Transpose
        t() %>%
        # Change back to data frame
        as.data.frame() %>%
        # Rename coefficients
        dplyr::rename(coefficient = 1)

    # Final output list with model information
    output <- list(model = params[["model"]],
                   coefficients = coefs)

    # If survival model, add time horizon
    if(params[["model"]] %in% c("aft", "cox", "fine-gray")) output[["horizon"]] <- params[["horizon"]]

    # If model is Cox or Fine-Gray and all hazards were requested, get those too
    if(params[["model"]] %in% c("cox", "fine-gray") & params[["all_hazards"]]){
        all_bhs <- do.call("cbind", lapply(models, \(x) as.data.frame(x[["all_hazards"]]))) %>%
            # Transpose
            t() %>%
            # Change to data frame
            as.data.frame() %>%
            # Final values are mean of each column
            dplyr::mutate(dplyr::across(tidyselect::everything(), mean)) %>%
            # Keep only first row
            dplyr::slice(1L) %>%
            # Transpose back
            t() %>%
            # Change to data frame
            as.data.frame() %>%
            # Add time
            dplyr::mutate(time = 1:params[["horizon"]]) %>%
            # Change name of first column
            dplyr::rename(bh = 1)

        # Add all baseline hazards to output
        output[["hazards"]] <- all_bhs
    }

    # If model is AFT, get scale too
    if(params[["model"]] == "aft"){
        # Get mean scale
        aft_scale <- do.call("rbind", lapply(models, \(x) x[["aft_scale"]])) %>%
            # Take mean
            mean()

        # Add scale to output
        output[["scale"]] <- aft_scale

        # Add distribution to output
        output[["distribution"]] <- params[["aft_distribution"]]
    }

    # If optimism correction, add unpenalized coefficients too
    if(params[["shrinkage"]]){
        # Get coefficients
        coefs <- do.call("rbind", lapply(models, \(x) x[["unpenalized_coefs"]])) %>%
            # Change to data frame
            as.data.frame() %>%
            # Final values are mean of each column
            dplyr::mutate(dplyr::across(tidyselect::everything(), mean)) %>%
            # Keep only first row
            dplyr::slice(1L) %>%
            # Transpose
            t() %>%
            # Change back to data frame
            as.data.frame() %>%
            # Rename coefficients
            dplyr::rename(coefficient = 1)

        # Add to output
        output[["unpenalized_coefficients"]] <- coefs

        # Add mean shrinkage factor to output
        output[["shrinkage_factor"]] <- mean(sapply(models, \(x) x[["shrinkage_factor"]]))
    }

    # Add call parameters to model output
    output[["params"]] <- params

    # Add formula to model output
    output[["formula"]] <- formula

    # Set output class
    class(output) <- "pm"

    # Return output
    return(output)
}
