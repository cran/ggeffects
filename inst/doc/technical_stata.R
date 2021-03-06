## ----set-options, echo = FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", dev = "png", fig.width = 7, fig.height = 5, message = FALSE, warning = FALSE)
if (!requireNamespace("magrittr", quietly = TRUE)) {
  knitr::opts_chunk$set(eval = FALSE)
}

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(magrittr)
set.seed(5)

data <- data.frame(
  outcome = rbinom(100, 1, 0.5),
  var1 = rbinom(100, 1, 0.1),
  var2 = rnorm(100, 10, 7)
)

m <- glm(
  outcome ~ var1 * var2, 
  data = data, 
  family = binomial(link = "logit")
)

## ----eval=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  use data.dta, clear
#  quietly logit outcome c.var1##c.var2
#  quietly margins, at(var2 = (-8(0.5)28) var1 = (0 1))
#  marginsplot

## ----out.width="100%", echo=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::include_graphics("vignette-stata-1.png", dpi = 72)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(ggeffects)
ggpredict(m, c("var2", "var1")) %>% plot()

