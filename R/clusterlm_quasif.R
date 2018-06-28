clusterlm_quasif = function (formula, data, method, test, threshold, np,
                             P, S, rnd_rotation, aggr_FUN, E, H, cl, multcomp, alpha, p_scale,
                             return_distribution, ndh, coding_sum, new_method) {
  if (is.null(method)) {
    method = "terBraak"
  }

  if (is.null(aggr_FUN)) {
    fun_name = "the sum"
    aggr_FUN = function(x) sum(x)
    } else {
    fun_name = "a user-defined function"
  }


  switch(method,
         terBraak_logp ={funP = function(...){cluster_quasif_terBraak_logp(...)}},{
           warning(paste("the method", method, "is not defined."))
           funP = function(...) {
             eval(parse(text = paste("cluster_quasif_", method, "(...)",
                                     sep = "", collpase = "")))}})


  if (!(class(formula[[2]]) == "matrix")) {
    formula[[2]] <- call("as.matrix", formula[[2]])
  }

  terms <- terms(formula, special = "Error", data = data)
  ind_error <- attr(terms, "specials")$Error
  error_term1 <- attr(terms, "variables")[[1 + ind_error[1]]]
  error_term2 <- attr(terms, "variables")[[1 + ind_error[2]]]

  formula_f <- update(formula, paste(". ~ .-", deparse(error_term1,
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


  ####withour reponse
  formula_allfixed_design = delete.response(terms(update(formula_allfixed, ~.)))
  formula_f_design = delete.response(terms(update(formula_f, ~.)))



  ####model.frame
  mf <- model.frame(formula = formula_allfixed, data = data)
  mf_design <- model.frame(formula = formula_allfixed_design, data = data)
  if (coding_sum) {
    mf_design <- permuco:::changeContrast(mf_design, contr = "contr.sum")
  }

  mf_f <- model.frame(formula = formula_f_design, data = mf_design)


  mf_id1 <- model.frame(formula = formula_id1, data = as.data.frame(lapply(mf_design,
                                                                           function(col) {
                                                                             col = as.factor(col)
                                                                             contrasts(col) = contr.sum
                                                                             col
                                                                           })))
  mf_id2 <- model.frame(formula = formula_id2, data = as.data.frame(lapply(mf_design,
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

  mm0 = model.matrix(formula_f_design,mf0)
  mm0_id1 = model.matrix(formula_id1,mf0)[,-1]
  mm0_id2 = model.matrix(formula_id2,mf0)[,-1]


  name <- colnames(mm_f)
  permuco:::checkBalancedData(fixed_formula = formula_f_design, data = mf)

  tf = delete.response(terms(update(formula_f, ~.)))


  zm = Zmat(mm0 = mm0, mm = mm_f, link = link,mm_id1 = mm_id1, mm_id2 = mm_id2,
            mm0_id1= mm0_id1, mm0_id2= mm0_id2, terms_f= tf)


  if (is.null(P)) {
    P = Pmat(np = np, n = ncol(mm_f))
  }

  if (is.null(S)) {
    S = Pmat(np = np, n = ncol(zm$z0),type = "coinflip")
  }

  np = permuco:::np.Pmat(S)

  gamma = zm$coding%*%((qr.coef(zm$qr_xz,y)[-c(1:ncol(mm_f)),]))

  args <- list(y = y, mm = mm_f, mm_id1 = mm_id1, mm_id2 = mm_id2, link = link,
               P = P, S = S, zm = zm, gamma = gamma)


  multiple_comparison <- list()
  length(multiple_comparison) <- max(attr(mm_f, "assign"))
  names(multiple_comparison) <- attr(attr(mf_f, "terms"), "term.labels")

  for (i in 1:max(attr(mm_f, "assign"))) {
    args$i = i
    distribution = funP(args = args)

    pvalue <- apply(distribution, 2, function(col) permuco:::compute_pvalue(distribution = col))
    multiple_comparison[[i]]$uncorrected = list(main = cbind(statistic = distribution[1,],
                                                             pvalue = pvalue))
    if (return_distribution){
      multiple_comparison[[i]]$uncorrected$distribution = distribution
    }

    multiple_comparison[[i]] = c(multiple_comparison[[i]],
                                 permuco:::switch_multcomp(multcomp = c("clustermass",multcomp),
                                                           distribution = distribution, threshold = threshold,
                                                           aggr_FUN = aggr_FUN, laterality = "bilateral",
                                                           E = E, H = H, ndh = ndh, pvalue = pvalue, alpha = alpha))
  }
  multiple_comparison = multiple_comparison[order(link[3, ],
                                                  link[1, ])]
  cluster_table <- permuco:::cluster_table(multiple_comparison)

  attr(cluster_table, "type") <- paste("Permutation test using ",
                                       method, " to handle noise variable and ", np, " permutations.")


  out = list()
  out$y = y
  out$model.matrix = mm_f
  out$model.matrix_id1 = mm_id1
  out$model.matrix_id2 = mm_id2
  out$zm = zm
  out$link = link
  out$P = P
  out$S = S
  out$cluster_table = cluster_table
  out$multiple_comparison = multiple_comparison
  out$data = mf
  out$method = method
  out$alpha = alpha
  out$multcomp = multcomp
  out$threshold = threshold
  out$test = test
  out$fun_name <- fun_name
  class(out) <- "clusterlm"
  return(out)
}
