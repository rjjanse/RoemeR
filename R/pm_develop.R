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
        stop("'model' should be one of 'linear', 'logistic', 'poisson', 'aft', 'cox', or 'fine-gray'.")
    }

    # If the model is a survival model, horizon should be defined
    if(params[["model"]] %in% c("fine-gray", "cox", "aft") & is.null(params[["horizon"]])){
        stop("Define 'horizon' in params for survival models")
    }

    # If the model is an accelerated failure time model, distribution should be defined
    if(params[["model"]] == "aft" & is.null(params[["aft_distribution"]])){
        stop("Define 'aft_distribution' in params for accelerated failure time models")
    }

    # Set argument values ----
    # Export all baseline hazards for Cox & Fine-Gray model
    if(is.null(params[["all_baseline_hazards"]])) params[["all_baseline_hazards"]] <- FALSE

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

    # Get final model ----
    # Fit models for each imputation
    models <- lapply(unique(data[[".imp"]]), \(x) develop_model(formula, data, x, params))

    # Return
    return(models)
}
