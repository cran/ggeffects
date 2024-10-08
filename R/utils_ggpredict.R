.validate_type_argument <- function(model, type, marginaleffects = FALSE, emmeans_call = FALSE) {
  # marginaleffects supports the predict-method types
  # we need a different approach to validation here
  if (marginaleffects) {
    # for zero-inflation models, we need to find the correct name
    # for the type argument...
    is_zero_inflated <- insight::model_info(model)$is_zero_inflated
    if (is_zero_inflated) {
      if (inherits(model, "glmmTMB")) {
        types <- c("conditional", "zprob")
      } else {
        types <- c("count", "zero")
      }
    }
    # first, we overwrite the "default"
    if (type == "fixed") {
      if (is_zero_inflated) {
        type <- types[1]
      } else if (class(model)[1] %in% .default_type$class) {
        type <- .default_type$type[.default_type$class == class(model)[1]]
      } else {
        type <- "response"
      }
    } else if (type %in% c("zi", "zero_inflated", "fe.zi")) {
      type <- "response"
    } else if (type %in% c("zi.prob", "zi_prob")) {
      type <- types[2]
    }
    # check which types are supported by the model's predict-method
    type_options <- .typedic$type[.typedic$class == class(model)[1]]
    if (!type %in% c("response", type_options)) {
      insight::format_error(sprintf(
        "`type = \"%s\"` is not supported. Please use %s%s.",
        type,
        if (length(type_options) > 1) "one of " else "",
        toString(paste0("`", type_options, "`"))
      ))
    }
    return(type)
  }

  # if we call "predict()" or "emmeans()", we have these different options
  if (emmeans_call) {
    type_choices <- c(
      "fe", "fixed", "count", "re", "random", "fe.zi", "zero_inflated",
      "re.zi", "zi_random", "zero_inflated_random", "zi.prob", "zi_prob"
    )
  } else {
    type_choices <- c(
      "fe", "fixed", "count", "re", "random", "fe.zi", "zero_inflated", "re.zi",
      "zi_random", "zero_inflated_random", "zi.prob", "zi_prob", "sim",
      "simulate", "surv", "survival", "cumhaz", "cumulative_hazard", "sim_re",
      "simulate_random", "debug"
    )
  }
  type <- match.arg(type, choices = type_choices)

  switch(type,
    fe = ,
    count = "fixed",
    re = "random",
    zi = ,
    fe.zi = "zero_inflated",
    re.zi = ,
    zi_random = "zero_inflated_random",
    zi.prob = "zi_prob",
    surv = "survival",
    cumhaz = "cumulative_hazard",
    sim = "simulate",
    sim_re = "simulate_random",
    type
  )
}


.retrieve_type_option <- function(model) {
  # retrieve model object's predict-method prediction-types (if any)
  predict_method <- .safe(lapply(
    class(model), function(i) {
      utils::getS3method("predict", i)
    }
  ))
  # check whether model class has a predict method
  if (!is.null(predict_method)) {
    predict_method <- predict_method[!vapply(predict_method, is.null, TRUE)][[1]]
  }
  # retrieve model object's predict-method prediction-types (if any)
  .safe(suppressWarnings(eval(formals(predict_method)$type)))
}


.back_transform_response <- function(model, mydf, back.transform, response.name = NULL, verbose = TRUE) {
  # skip if no information available
  if (is.null(model) && is.null(response.name)) {
    return(mydf)
  }

  # skip for multivariate response models
  if (insight::is_multivariate(model)) {
    # tell user
    if (verbose) {
      insight::format_alert("Back-transforming response variables is not carried out for multivariate-response models.")
    }
    return(mydf)
  }

  # we need the string of the response variable, to get information about transformation
  if (is.null(response.name)) {
    rv <- insight::find_terms(model)[["response"]]
  } else {
    # for pool_predictions(), we have no model object, but the response-string
    rv <- response.name
  }

  # for pool_predictions(), we have no model object, but rather the response-string
  if (is.null(model)) {
    # find possible transformations from response-string
    transformation <- insight::find_transformation(rv)
  } else {
    # find possible transformations from model
    transformation <- insight::find_transformation(model)
  }

  # skip if no information available, or no transformation applies
  if (is.null(transformation) || identical(transformation, "identity")) {
    return(mydf)
  }

  # transformed response, but no back-transform requested?
  # Tell user and return untransformed predictions
  if (!back.transform) {
    if (verbose) {
      insight::format_alert(
        paste0("Model has ", transformation, " transformed response. Predictions are on transformed scale.")
      )
    }
    return(mydf)
  }

  # for pool_predictions(), we have no model object, but rather the response-string
  if (is.null(model)) {
    # get inverse transformation function response-string
    trans_fun <- insight::get_transformation(rv, verbose = verbose)$inverse
  } else {
    # get inverse transformation function from model
    trans_fun <- insight::get_transformation(model, verbose = verbose)$inverse
  }

  # Tell user about transformation
  if (verbose && !is.null(trans_fun)) {
    insight::format_alert(
      paste0("Model has ", transformation, "-transformed response. Back-transforming predictions to original response scale. Standard errors are still on the transformed scale.") # nolint
    )
  }

  if (startsWith(transformation, "sqrt")) {
    # handle sqrt-transformed response separately - might be "sqrt(x + 1)"
    plus_minus <- eval(parse(text = gsub("sqrt\\(([^,\\+)]*)(.*)\\)", "\\2", rv)))
    if (is.null(plus_minus)) plus_minus <- 0
    mydf$predicted <- mydf$predicted^2 - plus_minus
    if (all(c("conf.low", "conf.high") %in% colnames(mydf))) {
      mydf$conf.low <- mydf$conf.low^2 - plus_minus
      mydf$conf.high <- mydf$conf.high^2 - plus_minus
    }
  } else if (startsWith(transformation, "log(") && transformation != "log-log") {
    # handle log-transformed response separately - might be "log(x + 1)"
    plus_minus <- eval(parse(text = gsub("log\\(([^,\\+)]*)(.*)\\)", "\\2", rv)))
    if (is.null(plus_minus)) plus_minus <- 0
    mydf$predicted <- exp(mydf$predicted) - plus_minus
    if (all(c("conf.low", "conf.high") %in% colnames(mydf))) {
      mydf$conf.low <- exp(mydf$conf.low) - plus_minus
      mydf$conf.high <- exp(mydf$conf.high) - plus_minus
    }
  } else if (!is.null(trans_fun)) {
    mydf$predicted <- trans_fun(mydf$predicted)
    if (all(c("conf.low", "conf.high") %in% colnames(mydf))) {
      mydf$conf.low <- trans_fun(mydf$conf.low)
      mydf$conf.high <- trans_fun(mydf$conf.high)
    }
  }

  mydf
}


.back_transform_data <- function(model, mydf, back.transform, response.name = NULL) {
  # skip if no information available
  if (is.null(mydf)) {
    return(NULL)
  }
  if (back.transform) {
    return(mydf)
  }

  # check if outcome is log-transformed, and if so,
  # back-transform predicted values to response scale
  if (is.null(response.name)) {
    rv <- insight::find_terms(model)[["response"]]
  } else {
    rv <- response.name
  }

  # sanity check
  if (!"response" %in% colnames(mydf)) {
    return(mydf)
  }

  # find transformation
  if (is.null(model)) {
    # find possible transformations from response-string
    transformation <- insight::find_transformation(rv)
  } else {
    # find possible transformations from model
    transformation <- insight::find_transformation(model)
  }

  # get transformation function
  if (is.null(model)) {
    # get inverse transformation function response-string
    trans_fun <- insight::get_transformation(rv, verbose = FALSE)$transformation
  } else {
    # get inverse transformation function from model
    trans_fun <- insight::get_transformation(model, verbose = FALSE)$transformation
  }

  if (!is.null(transformation) && !identical(transformation, "identity") && !is.null(trans_fun)) {
    if (startsWith(transformation, "sqrt")) {
      # handle sqrt-transformed response separately - might be "sqrt(x + 1)"
      plus_minus <- eval(parse(text = gsub("sqrt\\(([^,\\+)]*)(.*)\\)", "\\2", rv)))
      if (is.null(plus_minus)) plus_minus <- 0
      mydf$response <- sqrt(mydf$response) + plus_minus
    } else if (startsWith(transformation, "log(") && transformation != "log-log") {
      # handle log-transformed response separately - might be "log(x + 1)"
      plus_minus <- eval(parse(text = gsub("log\\(([^,\\+)]*)(.*)\\)", "\\2", rv)))
      if (is.null(plus_minus)) plus_minus <- 0
      mydf$response <- log(mydf$response) + plus_minus
    } else if (!is.null(trans_fun)) {
      mydf$response <- trans_fun(mydf$response)
    }
  }

  mydf
}
