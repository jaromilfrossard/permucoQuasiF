#' Integration of the quasi-F statistics into the \code{clusterlm} function from the \code{permuco} package.
#'
#' @description Provides cluster-mass tests based on permutations for quasi-F statistics
#' @param formula A formula object. The formula for the quasi-F statistic should be written using the notation
#' \code{+ Error(id/(A*C)) + Error(id/(B*C))} to specify random effect associate to subjects and items, where \code{A} is a within \code{id} factor, \code{B} is a within \code{item} factor and \code{C} is within \code{id} and \code{item} factors.
#' @param data similar to the \code{permuco} package.
#' @param np The number of permutations. Default value is \code{5000}.
#' @param method similar to the \code{permuco} package. The method is set to \code{"terBraak"} and \code{"terBraak_logp"} are available for the quasi-F statistic.
#' @param test similar to the \code{permuco} package. the argument does not influence the quasi-F statistic.
#' @param threshold similar to the \code{permuco} package.
#' @param aggr_FUN similar to the \code{permuco} package.
#' @param multcomp similar to the \code{permuco} package.
#' @param effect a numeric specifying the effect to test. The value correspond to the "assign" attribute of the model.matrix argument of the fixed effect. The default is NULL tests all effects.
#' @param ... Futher arguments, see details.
#' @return A list containing : a table of the clusters, or a \code{multcomp} object for the other multiple comparison procedures. Use the \link{plot.clusterlm} method to have a quick overview of the results.
#' @details
#' Similar to the \code{permuco} package.
#'@author jaromil.frossard@unige.ch
#'@import Matrix
#'@export
clusterlm <- function(formula, data=NULL, np = 5000, method = NULL, test = "fisher", threshold = NULL, aggr_FUN = NULL,
                      multcomp = "clustermass", effect = NULL,...){

  cl = match.call()
  if(is.null(data)){data <- model.frame(formula = formula)}



  ############
  #Formula CHECK
  Terms <- terms(formula, special = "Error", data = data)
  indError <- attr(Terms, "specials")$Error

  #dotargs
  dotargs = list(...)

  ####other parameters

  if(is.null(dotargs$alpha)){
    dotargs$alpha = 0.05
  }

  if(is.null(dotargs$p_scale)){
    dotargs$p_scale = F
  }

  if(is.null(dotargs$H)){
    switch(test,
           "t" = {dotargs$H = 2},
           "fisher" = {dotargs$H = 1})
  }

  if(is.null(dotargs$E)){
    dotargs$E = 0.5
  }


  if(is.null(dotargs$ndh)){
    dotargs$ndh = 500
  }

  if(is.null(dotargs$return_distribution)){
    dotargs$return_distribution = F
  }

  # if(is.null(threshold)){
  #   switch(test,
  #          "t" = {threshold = 2},
  #          "fisher" = {threshold = 4})
  # }

  if(is.null(dotargs$new_method)){
    dotargs$new_method = F
  }

  if(is.null(dotargs$coding_sum)){
    switch(test,
           "t" = {dotargs$coding_sum = F},
           "fisher" = {dotargs$coding_sum = T})
  }

  ###switch fix effet
  if (is.null(indError)) {
    result <- permuco:::clusterlm_fix( formula = formula, data = data, method = method, test = test, np = np,
                             P = dotargs$P, rnd_rotation = dotargs$rnd_rotation, aggr_FUN = aggr_FUN,
                             E = dotargs$E, H = dotargs$H, threshold = threshold,
                             return_distribution = dotargs$return_distribution, cl = cl, multcomp = multcomp,
                             alpha = dotargs$alpha, p_scale = dotargs$p_scale, coding_sum = dotargs$coding_sum,ndh = dotargs$ndh,
                             new_method = dotargs$new_method)
  } else if (!is.null(indError)&(length(indError)==1)){
    if(test!="fisher"){
      warning("Random effects model only accept fisher statistics. Test statistic is set to fisher.")
      test="fisher"}
    result <- permuco:::clusterlm_rnd( formula = formula, data = data, method = method, test = test, np = np,
                             P = dotargs$P, rnd_rotation = dotargs$rnd_rotation, aggr_FUN = aggr_FUN,
                             E = dotargs$E, H = dotargs$H, threshold = threshold,
                             return_distribution = dotargs$return_distribution, cl = cl, multcomp = multcomp,
                             alpha = dotargs$alpha, p_scale = dotargs$p_scale, coding_sum = dotargs$coding_sum,ndh = dotargs$ndh,
                             new_method = dotargs$new_method)}
  else if (!is.null(indError)&(length(indError)==2)) {
    if (test != "fisher") {
      warning("Random effect model only accept fisher type test statitics. test statistic set to fisher")
      test = "fisher"
    }
    result <- clusterlm_quasif(formula = formula, data = data,
                               method = method, test = test, np = np, P = dotargs$P,S = dotargs$S,
                               rnd_rotation = dotargs$rnd_rotation, aggr_FUN = aggr_FUN,
                               E = dotargs$E, H = dotargs$H, threshold = threshold,
                               return_distribution = dotargs$return_distribution,
                               cl = cl, multcomp = multcomp, alpha = dotargs$alpha,
                               p_scale = dotargs$p_scale, coding_sum = dotargs$coding_sum, ndh = dotargs$ndh,
                               new_method = dotargs$new_method, effect = effect)
  }

  ###output
  return(result)
}