utils::globalVariables(c("std.error", "conf.low", "conf.high", "datacol", "models", "avgpred"))

#' @rdname ggpredict
#' @importFrom dplyr group_by_ summarise arrange_ mutate ungroup select_
#' @importFrom stats predict loess lm
#' @importFrom tidyr nest unnest
#' @importFrom sjmisc var_rename
#' @export
ggaverage <- function(model, terms, ci.lvl = .95, type = c("fe", "re"), typical = c("mean", "median"), ...) {
  # get predictions for full data
  dat <- ggpredict(model, terms, ci.lvl, type, full.data = TRUE, typical, ...)
  # and for unique data, to get proper SE / CI for plotting
  dat2 <- ggpredict(model, terms, ci.lvl, type, full.data = FALSE, typical, ...)


  # remove columns with obs and resid
  dat <- dplyr::select_(dat, "-observed", "-residuals")


  # do summary, depending on third group
  if (tibble::has_name(dat, "facet")) {
    dat <- dplyr::group_by_(dat, "x", "group", "facet")
  } else {
    dat <- dplyr::group_by_(dat, "x", "group")
  }

  # is x a factor?
  xif <- attr(dat, "x.is.factor", exact = TRUE)
  x_is_factor <- !is.null(xif) && xif == "1"


  # check if x is a factor. We can then simply take mean values for
  # standard error and confidence intervals. For continuous variables,
  # we have no linear trend for the predicted values, and hence no "proper"
  # confidence bands that can be plotted with geom_ribbon. So we do some
  # workaround in such cases (see below)
  if (x_is_factor) {
    zus <- get_average_values(dat)
  } else {
    zus <- get_smoothed_avg(dat, dat2, ci.lvl)
  }


  # to tibble
  zus <- tibble::as_tibble(zus)


  # add back attributes. therefore, we first copy the attributes from the
  # original data frame and delete the common attributes, like class etc.
  # and then add attributes to our final df
  a <- attributes(dplyr::ungroup(dat))
  a[c("names","row.names","class","dim","dimnames")] <- NULL
  attributes(zus) <- c(attributes(zus), a)
  # no full data for averages
  attr(zus, "full.data") <- "0"

  # add class attribute
  class(zus) <- c("ggeffects", class(zus))

  zus
}


#' @rdname ggpredict
#' @export
ame <- function(model, terms, ci.lvl = .95, type = c("fe", "re"), typical = c("mean", "median"), ...) {
  ggaverage(model, terms, ci.lvl, type, typical, ...)
}


# this method simply computes the mean values of predictions, se and ci
# to get average marginal effects for categorical variables
# As categorical variables do not necessarily need to follow a "linear"
# (or other distribution based) trend, we can simply calculate the mean
get_average_values <- function(dat) {
  # get average values for predictions, SE and CI
  zus <- dplyr::summarise(
    dat,
    predicted = mean(predicted),
    std.error = mean(std.error, na.rm = TRUE),
    conf.low = mean(conf.low, na.rm = TRUE),
    conf.high = mean(conf.high, na.rm = TRUE)
  ) %>%
    dplyr::ungroup()

  # sort columns
  if (tibble::has_name(dat, "facet")) {
    zus <- zus[, c(1, 4:7, 2:3)]
  } else {
    zus <- zus[, c(1, 3:6, 2)]
  }

  zus
}


# this method prepares the data to get smoothed predictions for
# average effects
get_smoothed_avg <- function(dat, dat2, ci.lvl) {
  # get average prediction. this is quite strait forward...
  zus <- dplyr::summarise(dat, predicted = mean(predicted)) %>% dplyr::ungroup()


  # remove unused columns for 2nd data frame. we will bind some
  # columns from this df to our final df
  dat2 <- dplyr::select_(dat2, "-predicted", "-conf.low", "-conf.high")


  # arrange columns, so we have equal row order in both data frames
  # then, join data frames, so we have proper standard errors. we cannot
  # average standard errors and conf. int. because of variation within groups.
  # hence, we have no linear trend of predictions, nor for SE and CI. this
  # would result in unproper CI-bands when plotting with geom_ribbon.
  # so instead, we compute a regression on the predictions, based by the
  # variable "x" (the same family and link function as in the orignal model),
  # take the fitted values of predictions and compute the CI based on these
  # values.
  if (tibble::has_name(dat, "facet")) {
    zus <- dplyr::arrange_(zus, "x", "group", "facet")
    dat2 <- dplyr::arrange_(dat2, "x", "group", "facet")
    zus <- dplyr::left_join(zus, dat2, by = c("x", "group", "facet"))
    zus <- dplyr::group_by_(zus, "group", "facet")
  } else {
    zus <- dplyr::arrange_(zus, "x", "group")
    dat2 <- dplyr::arrange_(dat2, "x", "group")
    zus <- dplyr::left_join(zus, dat2, by = c("x", "group"))
    zus <- dplyr::group_by_(zus, "group")
  }

  get_smoothed_predictions(zus, dat, ci.lvl)
}


# this method computes the standard errors and confidence intervals
# for continuous variables. since the averaged predicted values follow a specific
# "trend" (linear, curvilinear), the standard errors and ci should do so, too.
# however, due to different spread of data points, se and ci may follow a "zigzag" curve.
# to prevent this, we "predict" the (curvi-)linear trend of the average predicted
# values, and then take the standard errors from the marginal effects at the mean
# and calculate ci manually. This allows us to plot a "geom_ribbon()", without
# needing to use a smoother for the confidence bands.
get_smoothed_predictions <- function(zus, dat, ci.lvl) {
  # now we need to compute confidence intervals, by predicting
  # values from the averaged predictions. This is how to get
  # "smoothed" predictions from averages
  fam <- attr(dat, "family", exact = TRUE)
  link <- attr(dat, "link", exact = TRUE)


  # since average marginal effects may vary in their slopes by group,
  # we need to predict the values for each group separately.
  if (fam == "gaussian" && link == "identity") {
    # for linear models, compute linear trend
    zus <- zus %>%
      tidyr::nest(.key = datacol) %>%
      dplyr::mutate(
        models = purrr::map(datacol, ~stats::lm(
          formula = predicted ~ x,
          data = .x
        )))
  } else {
    # furthermore, we can't compute a glm on predicted values of a glm - so we use
    # instead a local smoother to achieve predicted values for average effects
    zus <- zus %>%
      tidyr::nest(.key = datacol) %>%
      dplyr::mutate(
        models = purrr::map(datacol, ~stats::loess(
          formula = predicted ~ x,
          family = "symmetric",
          degree = 1,
          data = .x
        )))
  }


  # now compute new predictions of average marginal effects
  zus <- zus %>%
    dplyr::mutate(avgpred = purrr::map(models, ~as.vector(stats::predict(.x)))) %>%
    dplyr::select_("-models") %>%
    tidyr::unnest() %>%
    dplyr::ungroup()


  # if user does not want CI, we don't need to compute the predictions
  if (!is.null(ci.lvl) && !is.na(ci.lvl)) {
    # bind std.error and CI to final data frame
    zus <- zus %>%
      dplyr::mutate(
        conf.low = avgpred - stats::qnorm(.975) * std.error,
        conf.high = avgpred + stats::qnorm(.975) * std.error
      )
  } else {
    zus$std.error <- NA
    zus$conf.low <- NA
    zus$conf.high <- NA
  }

  # remove standard predicted, and replace with avg predicted
  zus <- zus %>%
    dplyr::select_("-predicted") %>%
    sjmisc::var_rename(avgpred = "predicted")

  # re-arrange columns
  if (tibble::has_name(dat, "facet"))
    zus <- zus[, c(4, 3, 5:7, 1:2)]
  else
    zus <- zus[, c(3, 2, 4:6, 1)]

  zus
}
