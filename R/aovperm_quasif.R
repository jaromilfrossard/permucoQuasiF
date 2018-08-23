#' @importFrom  stats contrasts contrasts<- model.response delete.response model.matrix as.formula
#' @importFrom permuco Pmat
aovperm_quasif = function(formula, data, method, np, P, S, coding_sum, rnd_rotation, new_method = NULL,effect = NULL) {
  if (is.null(coding_sum)) {
    coding_sum = T
  }
  if (is.null(new_method)) {
    new_method = F
  }
  if (is.null(method)) {
    method = "terBraak"
  }
  if (!new_method) {
    method = match.arg(method, c("terBraak",
                                 "terBraak_logp"))
  }


  switch(method,
         terBraak_logp ={funP = function(...){quasif_terBraak_logp(...)}},
         terBraak ={funP = function(...){quasif_terBraak(...)}},{
           warning(paste("the method", method, "is not defined."))
           funP <- function(...) {
             eval(parse(text = paste("quasif_", method, "(...)",
                                     sep = "", collpase = "")))}})

  terms <- terms(formula, special = "Error", data = data)
  ind_error <- attr(terms, "specials")$Error
  error_term1 <- attr(terms, "variables")[[1 + ind_error[1]]]
  error_term2 <- attr(terms, "variables")[[1 + ind_error[2]]]

  formula_f <- update.formula(formula, paste(". ~ .-", deparse(error_term1,
                                                       width.cutoff = 500L, backtick = TRUE),"-",deparse(error_term2,
                                                                                                         width.cutoff = 500L, backtick = TRUE) ))
  e_term1 <- deparse(error_term1[[2L]], width.cutoff = 500L,
                     backtick = TRUE)
  e_term2 <- deparse(error_term2[[2L]], width.cutoff = 500L,
                     backtick = TRUE)
  formula_allfixed <- as.formula(paste(c(formula_f[[2]], "~",
                                         formula_f[[3]], "+", e_term1,"+",e_term2), collapse = ""))
  formula_within1 <- formula(paste("~", e_term1, collapse = ""))
  formula_within1 <- formula(paste("~", deparse(error_term1[[2]][[3]]),
                                   collapse = ""))
  formula_within2 <- formula(paste("~", e_term2, collapse = ""))
  formula_within2 <- formula(paste("~", deparse(error_term2[[2]][[3]]),
                                   collapse = ""))

  formula_id1 <- formula(paste("~", deparse(error_term1[[2]][[2]]),
                               collapse = ""))
  formula_id2 <- formula(paste("~", deparse(error_term2[[2]][[2]]),
                               collapse = ""))


  #formulas <<- list(id1 = formula_id1, id2 = formula_id2, fe = formula_f,within1 = formula_within1,within2 =formula_within2)


  ####model.frame
  mf <- model.frame(formula = formula_allfixed, data = data)
  if (coding_sum) {
    mf <- permuco:::changeContrast(mf, contr = contr.sum)
  }
  mf_f <- model.frame(formula = formula_f, data = mf)
  mf_id1 <- model.frame(formula = formula_id1, data = as.data.frame(lapply(mf,
                                                                           function(col) {
                                                                             col = as.factor(col)
                                                                             contrasts(col) = contr.sum
                                                                             col
                                                                           })))
  mf_id2 <- model.frame(formula = formula_id2, data = as.data.frame(lapply(mf,
                                                                           function(col) {
                                                                             col = as.factor(col)
                                                                             contrasts(col) = contr.sum
                                                                             col
                                                                           })))
  y <- model.response(mf)

  ###link
  link = link(formula_f = formula_f, formula_within1 = formula_within1,formula_within2 = formula_within2)

  mm_f <- model.matrix(attr(mf_f, "terms"), data = mf_f)
  #mm_id1 <- model.matrix(attr(mf_id1, "terms"), data = mf_id1)[,
  #                                                          -1, drop = F]
  #mm_id2 <- model.matrix(attr(mf_id2, "terms"), data = mf_id2)[,
  #                                                            -1, drop = F]

  lf <- list(model_frame = mf,formula_f = formula_f,formula_within = formula_within1,formula_id = formula_id1)

  mm_id1 <- mm_z(model_frame = mf,formula_f = formula_f,formula_within = formula_within1,formula_id = formula_id1)
  mm_id2 <- mm_z(model_frame = mf,formula_f = formula_f,formula_within = formula_within2,formula_id = formula_id2)



  mf0 = mf


  for(i in 1:NCOL(mf0)){
    if(is.factor(mf0[,i])){
      contrasts(mf0[,i],how.many = length(levels(mf0[,i]))) = contrasts(mf0[,i],contrasts=F)
    }else{
      mf0[,i]<-mf0[,i]
    }

  }

  mm0 = model.matrix(formula_f,mf0)
  mm0_id1 = model.matrix(formula_id1,mf0)[,-1]
  mm0_id2 = model.matrix(formula_id2,mf0)[,-1]


  name <- colnames(mm_f)
  permuco:::checkBalancedData(fixed_formula = formula_f, data = cbind(y,
                                                                      mf))

  tf = delete.response(terms(update.formula(formula_f, ~.)))

  zm = Zmat(mm0 = mm0, mm = mm_f, link = link,mm_id1 = mm_id1, mm_id2 = mm_id2,
            mm0_id1= mm0_id1, mm0_id2= mm0_id2, terms_f= tf)


  if (is.null(S)) {
    S = Pmat(np = np, n = ncol(zm$z0),type = "coinflip")
  }
  np = permuco:::np.Pmat(S)

  pry = as.numeric(zm$coding%*%(qr.coef(zm$qr_xz,y)[-c(1:ncol(mm_f))]))*S
  pry = zm$z0%*%Matrix(as.matrix(pry))
  args <- list(y = y, mm = mm_f, mm_id1 = mm_id1, mm_id2 = mm_id2,
               link = link, S = S, zm = zm, pry = pry)



  # args <- list(y = y, mm = mm_f, mm_id1 = mm_id1, mm_id2 = mm_id2, link = link,
  #              P = P,mm0 = mm0, mm0_id1 = mm0_id1, mm0_id2 = mm0_id2,terms_f = tf)
  # ag <<- args

  if(is.null(effect)){
    effect = 1:max(attr(mm_f, "assign"))
  }else{
    effect = sort(unique(effect))
  }


  distribution <- sapply(effect, function(i) {
    args$i = i
    funP(args = args)
  })

  distribution = matrix(distribution,nrow=np)


  colnames(distribution) = attr(attr(mf_f, "terms"), "term.labels")[effect]
  permuco:::check_distribution(distribution = distribution, digits = 10,
                               n_unique = 300)

  table = anova_table_quasif(args)

  rownames(table) = attr(attr(mf_f, "terms"), "term.labels")
  permutation_pvalue = apply(distribution, 2, function(d) {
    permuco:::compute_pvalue(distribution = d, laterality = "bilateral",
                             na.rm = T)
  })

  table$"permutation P(>F)" = NA
  table$"permutation P(>F)"[effect] = permutation_pvalue

  attr(table, "type") <- paste("Permutation test using", method,
                               "to handle noise variable and", np, "permutations.")

  out = list()
  out$y = y
  out$pry = pry
  out$zm =zm
  out$model.matrix = mm_f
  out$model.matrix_id1 = mm_id1
  out$model.matrix_id2 = mm_id2
  out$link = link
  out$P = P
  out$S = S
  out$np = np
  out$table = table
  out$distribution = distribution
  out$data = mf
  out$method = method
  class(out) <- "lmperm"
  return(out)


}
