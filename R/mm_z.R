#' Create the Z matrices based with sum contrasts
#'
#'@description Compute the z matrice for the random effects and reduced its rank and dimension with respect to the between effects.
#'
#'@param model_frame a model.frame object containint the facotors of interest
#'@param formula_f a formula of all the effects as fixed
#'@param formula_within a formula indicating the within effect
#'@param formula_id the formula indicating the sampling units
#'
#'@return a matrix with its contrasts as attributs.
#'@importFrom stats formula update.formula
mm_z = function(model_frame, formula_f, formula_within, formula_id){
  tf = terms(formula_f)
  tw = terms(formula_within)
  tid = terms(formula_id)

  term.labels_between =  (attr(tf,"term.labels")[attr(tf,"order")==1])[
    !(attr(tf,"term.labels")[attr(tf,"order")==1])%in%(attr(tw,"term.labels")[attr(tw,"order")==1])]

  formula_between = formula(paste("~ 0",paste(term.labels_between, collapse = ":"),collapse = "+"))
  formula_id = update.formula(old = formula_id,new = ~ 0 + .)

  mm0_id = model.matrix(formula_id,model_frame)
  mm0_b = model.matrix(formula_between,model_frame)
  #debug
  # z = lapply(1:NCOL(mm0_b),function(i){
  #   which_row = which(mm0_b[,i] == 1)
  #   which_col = which(colSums(mm0_id[which_row,])!=0)
  #
  #   z0i = mm0_id[which_row,which_col]%*%contr.sum(length(which_col))
  #   zi = matrix(0, ncol= NCOL(z0i),nrow= NROW(mm0_id))
  #   zi[which_row,] = z0i
  #   attr(zi,"contrasts") = contr.sum(length(which_col))
  #   zi
  # })
  #
  # c = as.matrix(bdiag(lapply(z,function(zi)attr(zi,"contrasts"))))
  # z = do.call("cbind",z)
  # attr(z,"contrasts") = c
  z = lapply(1:NCOL(mm0_b),function(i){
    which_row = which(mm0_b[,i] == 1)
    which_col = which(colSums(mm0_id[which_row,])!=0)

    z0i = mm0_id[which_row,which_col]%*%contr.sum(length(which_col))
    zi = matrix(0, ncol= NCOL(z0i),nrow= NROW(mm0_id))
    zi[which_row,] = z0i
    contr = matrix(0,nrow = ncol(mm0_id),ncol = length(which_col)-1)
    contr[which_col,] = contr.sum(length(which_col))
    attr(zi,"contrasts") = contr
    zi
  })

  lapply(z,function(zi)attr(zi,"contrasts"))

  c = do.call("cbind",lapply(z,function(zi)attr(zi,"contrasts")))
  z = do.call("cbind",z)
  attr(z,"contrasts") = c
  attr(z,"label") = attr(terms(formula_id),"term.labels")

  return(z)
}




























