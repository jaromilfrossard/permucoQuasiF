#' Integration of the quasi-F statistics into the \code{aovperm} function from the \code{permuco} package.
#'
#' @description  Provides p-values for omnibus tests based on permutations for quasi-F statistics
#' @param formula A formula object. The formula for the quasi-F statistic should be written using the notation
#' \code{+ Error(id/(A*C)) + Error(id/(B*C))} to specify random effect associate to subjects and items, where \code{A} is a within \code{id} factor, \code{B} is a within \code{item} factor and \code{C} is within \code{id} and \code{item} factors.
#' @param data similar to the \code{permuco} package.
#' @param np similar to the \code{permuco} package. Default value is \code{5000}.
#' @param method similar to the \code{permuco} package. The method is set to \code{"terBraak"} and for the quasi-F statistic.
#' @param ... Futher arguments, see details.
#'
#' @return A \code{lmperm} object containing most of the objects given in an \link{lm} object, an ANOVA table with parametric and permutation p-values, the test statistics and the permutation distributions.
#'
#' @details similar to the \code{permuco} package.
#' @author jaromil.frossard@unige.ch
#' @importFrom stats terms contr.sum model.frame terms
#' @export
aovperm<-function(formula, data=NULL, np = 5000, method = NULL,...){
  #method <- pmatch(method)

  if(is.null(data)){data <- model.frame(formula = formula)}

  #Formula CHECK
  Terms <- terms(formula, special = "Error", data = data)
  indError <- attr(Terms, "specials")$Error

  #check for intercept
  if(!attr(Terms,"intercept")){warning("Intercept should be specified in the formula")}

  #dotargs
  dotargs=list(...)

  ###switch fix effet
  if (is.null(indError)) {
    result <- permuco:::aovperm_fix( formula = formula, data = data, method = method, np = np, coding_sum = dotargs$coding_sum, P = dotargs$P,
                           rnd_rotation = dotargs$rnd_rotation, new_method = dotargs$new_method)
  } else if (!is.null(indError)&(length(indError)==1)){
    result <- permuco:::aovperm_rnd( formula = formula, data = data, method = method, np = np, coding_sum = dotargs$coding_sum, P = dotargs$P,
                           rnd_rotation = dotargs$rnd_rotation, new_method = dotargs$new_method)
    }else if (!is.null(indError)&(length(indError)==2)){
      result <- aovperm_quasif(formula = formula, data = data, method = method, np = np, coding_sum = dotargs$coding_sum,  P = dotargs$P, S = dotargs$S,
                               rnd_rotation = dotargs$rnd_rotation, new_method = dotargs$new_method)
}

  ###output
  return(result)
}
