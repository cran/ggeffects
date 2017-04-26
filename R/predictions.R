# select prediction method, based on model-object
select_prediction_method <- function(fun, model, expanded_frame, ci.lvl, type, binom_fam, ...) {
  if (fun == "svyglm") {
    # survey-objects -----
    fitfram <- get_predictions_svyglm(model, expanded_frame, ci.lvl, ...)
  } else if (fun == "svyglm.nb") {
    # survey-glm.nb-objects -----
    fitfram <- get_predictions_svyglmnb(model, expanded_frame, ci.lvl, ...)
  } else if (fun == "lrm") {
    # lrm-objects -----
    fitfram <- get_predictions_lrm(model, expanded_frame, ci.lvl, ...)
  } else if (fun %in% c("glm", "glm.nb")) {
    # glm-objects -----
    fitfram <- get_predictions_glm(model, expanded_frame, ci.lvl, ...)
  } else if (fun == "glmmTMB") {
    # glmTMB-objects -----
    fitfram <- get_predictions_glmmTMB(model, expanded_frame, ci.lvl, ...)
  } else if ((fun == "lmer" || (fun == "glmer" && binom_fam)) && !is.na(ci.lvl)) {
    # merMod-objects -----
    fitfram <- get_predictions_merMod(model, expanded_frame, ci.lvl, type, fun, ...)
  } else if (fun %in% c("lmer", "nlmer", "glmer")) {
    # merMod-objects, variant -----
    fitfram <- get_predictions_merMod2(model, expanded_frame, type, ...)
  } else if (fun == "gam") {
    # gam-objects -----
    fitfram <- get_predictions_gam(model, expanded_frame, ci.lvl, ...)
  } else if (fun == "vgam") {
    # gam-objects -----
    fitfram <- get_predictions_vgam(model, expanded_frame, ci.lvl, ...)
  } else if (fun == "lm") {
    # lm-objects -----
    fitfram <- get_predictions_lm(model, expanded_frame, ci.lvl, ...)
  } else if (fun %in% c("lme", "gls", "plm")) {
    # lme-objects -----
    fitfram <- get_predictions_lme(model, expanded_frame, ...)
  } else if (fun == "gee") {
    # gee-objects -----
    fitfram <- get_predictions_gee(model, expanded_frame, ...)
  } else {
    # general-objects -----
    fitfram <- get_predictions_generic(model, expanded_frame, ...)
  }

  fitfram
}


# predictions for survey objects

get_predictions_svyglm <- function(model, fitfram, ci.lvl, ...) {
  # does user want standard errors?
  se <- !is.null(ci.lvl) && !is.na(ci.lvl)

  prdat <-
    stats::predict(
      model,
      newdata = fitfram,
      type = "response",
      se.fit = se,
      level = ci.lvl,
      ...
    )

  # check if user wants standard errors
  if (se) {
    # get variance matrix for standard errors. "survey" stores the information
    # somewhat different from classical predict function
    vv <- attr(prdat, "var")
    # compute standard errors
    if (is.matrix(vv))
      prdat <- as.data.frame(cbind(prdat, sqrt(diag(vv))))
    else
      prdat <- as.data.frame(cbind(prdat, sqrt(vv)))

    # consistent column names
    colnames(prdat) <- c("fit", "se.fit")

    # copy predictions
    fitfram$predicted <- prdat$fit
    fitfram$std.error <- prdat$se.fit

    # calculate CI
    fitfram$conf.low <- prdat$fit - stats::qnorm(.975) * prdat$se.fit
    fitfram$conf.high <- prdat$fit + stats::qnorm(.975) * prdat$se.fit
  } else {
    # copy predictions
    fitfram$predicted <- as.vector(prdat)

    # no CI and no SE
    fitfram$std.error <- NA
    fitfram$conf.low <- NA
    fitfram$conf.high <- NA
  }

  fitfram
}


# predictions for glm

get_predictions_glm <- function(model, fitfram, ci.lvl, ...) {
  # does user want standard errors?
  se <- !is.null(ci.lvl) && !is.na(ci.lvl)

  prdat <-
    stats::predict.glm(
      model,
      newdata = fitfram,
      type = "response",
      se.fit = se,
      level = ci.lvl,
      ...
    )
  # copy predictions
  fitfram$predicted <- prdat$fit

  # did user request standard errors? if yes, compute CI
  if (se) {
    # standard error
    fitfram$std.error <- prdat$se.fit

    # calculate CI
    fitfram$conf.low <- prdat$fit - stats::qnorm(.975) * prdat$se.fit
    fitfram$conf.high <- prdat$fit + stats::qnorm(.975) * prdat$se.fit
  } else {
    # standard error
    fitfram$std.error <- NA

    # No CI
    fitfram$conf.low <- NA
    fitfram$conf.high <- NA
  }

  fitfram
}


# predictions for lrm

#' @importFrom stats plogis
get_predictions_lrm <- function(model, fitfram, ci.lvl, ...) {
  # does user want standard errors?
  se <- !is.null(ci.lvl) && !is.na(ci.lvl)

  prdat <-
    stats::predict(
      model,
      newdata = fitfram,
      type = "lp",
      se.fit = se,
      ...
    )
  # copy predictions
  fitfram$predicted <- stats::plogis(prdat$linear.predictors)

  # did user request standard errors? if yes, compute CI
  if (se) {
    # standard error
    fitfram$std.error <- prdat$se.fit

    # calculate CI
    fitfram$conf.low <- stats::plogis(prdat$linear.predictors - stats::qnorm(.975) * prdat$se.fit)
    fitfram$conf.high <- stats::plogis(prdat$linear.predictors + stats::qnorm(.975) * prdat$se.fit)
  } else {
    # standard error
    fitfram$std.error <- NA

    # No CI
    fitfram$conf.low <- NA
    fitfram$conf.high <- NA
  }

  fitfram
}


# predictions for svyglm.nb

get_predictions_svyglmnb <- function(model, fitfram, ci.lvl, ...) {
  # does user want standard errors?
  se <- !is.null(ci.lvl) && !is.na(ci.lvl)

  prdat <-
    stats::predict(
      model,
      newdata = fitfram,
      type = "response",
      se.fit = se,
      level = ci.lvl,
      ...
    )
  # copy predictions
  fitfram$predicted <- prdat$fit

  # did user request standard errors? if yes, compute CI
  if (se) {
    # standard error
    fitfram$std.error <- prdat$se.fit

    # calculate CI
    fitfram$conf.low <- prdat$fit - stats::qnorm(.975) * prdat$se.fit
    fitfram$conf.high <- prdat$fit + stats::qnorm(.975) * prdat$se.fit
  } else {
    # standard error
    fitfram$std.error <- NA

    # No CI
    fitfram$conf.low <- NA
    fitfram$conf.high <- NA
  }

  fitfram
}


# predictions for merMod

#' @importFrom merTools predictInterval
get_predictions_merMod <- function(model, fitfram, ci.lvl, type, fun, ...) {
  # prediction intervals from merMod-package only work for linear or
  # binary logistic multilevel models
  prdat <- suppressWarnings(
    merTools::predictInterval(
      model, newdata = fitfram, which = ifelse(type == "fe", "fixed", "full"),
      type = ifelse(fun == "lmer", "linear.prediction", "probability"),
      level = ci.lvl,
      ...
    ))
  # copy predictions
  fitfram$predicted <- prdat$fit

  # no standard error
  fitfram$std.error <- NA

  # calculate CI
  fitfram$conf.low <- prdat$lwr
  fitfram$conf.high <- prdat$upr

  fitfram
}


# predictions for glmmTMB

get_predictions_glmmTMB <- function(model, fitfram, ci.lvl, ...) {
  # does user want standard errors?
  se <- !is.null(ci.lvl) && !is.na(ci.lvl)

  prdat <-
    stats::predict(model,
                   newdata = fitfram,
                   zitype = "response",
                   se.fit = se,
                   ...)

  # did user request standard errors? if yes, compute CI
  if (se) {
    # copy predictions
    fitfram$predicted <- prdat$fit

    # standard error
    fitfram$std.error <- prdat$se.fit

    # calculate CI
    fitfram$conf.low <- prdat$fit - stats::qnorm(.975) * prdat$se.fit
    fitfram$conf.high <- prdat$fit + stats::qnorm(.975) * prdat$se.fit
  } else {
    # copy predictions
    fitfram$predicted <- as.vector(prdat)

    # no standard error
    fitfram$std.error <- NA

    # no CI
    fitfram$conf.low <- NA
    fitfram$conf.high <- NA
  }

  fitfram
}


# predictions for merMod, alternative

get_predictions_merMod2 <- function(model, fitfram, type, ...) {
  # for all other kinds of glmer or nlmer, we need the predict-function from
  # lme4, however, without ci-bands
  if (type == "fe")
    fitfram$predicted <-
      stats::predict(model,
                     newdata = fitfram,
                     type = "response",
                     re.form = NA,
                     ...)
  else
    fitfram$predicted <-
      stats::predict(model,
                     newdata = fitfram,
                     type = "response",
                     re.form = NULL,
                     ...)
  # no SE and CI for lme4-predictions
  fitfram$std.error <- NA
  fitfram$conf.low <- NA
  fitfram$conf.high <- NA

  fitfram
}


# predictions for gam

#' @importFrom prediction prediction
get_predictions_gam <- function(model, fitfram, ci.lvl, ...) {
  # call prediction
  prdat <- prediction::prediction(model, data = fitfram, at = NULL, type = "response", ...)

  # copy predictions
  fitfram$predicted <- prdat$fitted

  # no standard error
  fitfram$std.error <- NA

  fitfram
}


# predictions for vgam

#' @importFrom prediction prediction
get_predictions_vgam <- function(model, fitfram, ci.lvl, ...) {
  prdat <- prdat <-
    stats::predict(
      model,
      type = "response",
      ...
    )

  # copy predictions
  fitfram$predicted <- prdat$fitted

  # no standard error
  fitfram$std.error <- NA

  fitfram
}


# predictions for lm

get_predictions_lm <- function(model, fitfram, ci.lvl, ...) {
  # does user want standard errors?
  se <- !is.null(ci.lvl) && !is.na(ci.lvl)

  prdat <-
    stats::predict(
      model,
      newdata = fitfram,
      type = "response",
      se.fit = se,
      level = ci.lvl,
      ...
    )

  # did user request standard errors? if yes, compute CI
  if (se) {
    # copy predictions
    fitfram$predicted <- prdat$fit

    # standard error
    fitfram$std.error <- prdat$se.fit

    # calculate CI
    fitfram$conf.low <- prdat$fit - stats::qnorm(.975) * prdat$se.fit
    fitfram$conf.high <- prdat$fit + stats::qnorm(.975) * prdat$se.fit
  } else {
    # copy predictions
    fitfram$predicted <- as.vector(prdat)

    # no standard error
    fitfram$std.error <- NA

    # no CI
    fitfram$conf.low <- NA
    fitfram$conf.high <- NA
  }

  fitfram
}


# predictions for lme

get_predictions_lme <- function(model, fitfram, ...) {
  prdat <-
    stats::predict(
      model,
      newdata = fitfram,
      type = "response",
      ...
    )
  # copy predictions
  fitfram$predicted <- as.vector(prdat)

  # No standard error
  fitfram$std.error <- NA

  # No CI
  fitfram$conf.low <- NA
  fitfram$conf.high <- NA

  fitfram
}


# predictions for gee

get_predictions_gee <- function(model, fitfram, ...) {
  prdat <- prdat <-
    stats::predict(
      model,
      type = "response",
      ...
    )
  # copy predictions
  fitfram$predicted <- as.vector(prdat)

  # No standard error
  fitfram$std.error <- NA

  # No CI
  fitfram$conf.low <- NA
  fitfram$conf.high <- NA

  fitfram
}


# predictions for generic models

#' @importFrom prediction prediction
#' @importFrom tibble as_tibble
#' @importFrom sjmisc var_rename
get_predictions_generic <- function(model, fitfram, ...) {
  prdat <-
    prediction::prediction(
      model,
      data = fitfram,
      type = "response",
      ...
    ) %>%
    tibble::as_tibble() %>%
    sjmisc::var_rename(
      fitted = "predicted",
      se.fitted = "std.error"
    )

  # copy predictions
  fitfram$predicted <- prdat$fit

  # No standard error
  fitfram$std.error <- NA

  # No CI
  fitfram$conf.low <- NA
  fitfram$conf.high <- NA

  fitfram
}
