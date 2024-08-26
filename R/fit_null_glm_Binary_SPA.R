#' Fit generalized linear model under the null hypothesis for imbalanced case-control unrelated samples.
#'
#' The \code{fit_null_glm_Binary_SPA} function is a wrapper of the \code{\link{glm}} function from the
#' \code{\link{stats}} package that fits a regression model under the null hypothesis for
#' imbalanced case-control unrelated samples, which provides the preliminary step for subsequent
#' variant-set tests in whole genome sequencing data analysis. See \code{\link{glm}} for more details.
#' @param fixed an object of class \code{\link{formula}} (or one that can be coerced to that class):
#' a symbolic description of the fixed effects model to be fitted.
#' @param data a data frame or list (or object coercible by \code{\link{as.data.frame}} to a data frame)
#' containing the variables in the model.
#' @param family a description of the error distribution and link function to be used
#' in the model. This can be a character string naming a family function, a family
#' function or the result of a call to a family function. (See \code{\link{family}} for details of family functions).
#' Can be either "gaussian" for continuous phenotype or "binomial" for binary phenotype.
#' @param ... additional arguments that could be passed to \code{\link{glm}}.
#' @return The function returns an object of the model fit from \code{\link{glm}} (\code{obj_nullmodel}),
#' with additional elements indicating the samples are unrelated
#' (\code{obj_nullmodel$relatedness = FALSE}), and indicating the samples are under imbalanced case-control design (obj_nullmodel$use_SPA = TRUE).
#' See \code{\link{glm}} for more details.
#' @references Dey, R., et al. (2017). A fast and accurate algorithm to test for binary phenotypes
#' and its application to PheWAS. \emph{The American Journal of Human Genetics}, \emph{101}(1), 37-49.
#' (\href{https://doi.org/10.1016/j.ajhg.2017.05.014}{pub})
#' @export

fit_null_glm_Binary_SPA <- function(fixed, data, family = binomial(link = "logit"), ...){
  obj_nullmodel <- glm(formula = fixed, data = data, family = family, ...)
  obj_nullmodel$relatedness <- FALSE

  X <- model.matrix(obj_nullmodel)
  working <- obj_nullmodel$weights

  obj_nullmodel$XW <- t(X*working)
  obj_nullmodel$XXWX_inv <- X%*%solve(t(X*working)%*%X)

  obj_nullmodel$use_SPA <- TRUE

  return(obj_nullmodel)
}

