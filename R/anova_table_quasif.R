anova_table_quasif = function(args){
  table <- t(sapply(1:max(attr(args$mm, "assign")), function(i) {
    args$i = i
    effect_anova_quasif(args = args)
  }))
  table = as.data.frame(table)
  class(table) = c("lmpermutation_table", "data.frame")
  return(table)
}


effect_anova_quasif = function(args){
  link = args$link
  mm = Matrix(args$mm,sparse = T)
  mm_id1 = Matrix(args$mm_id1,sparse = T)
  mm_id2 = Matrix(args$mm_id2,sparse = T)

  assign = attr(args$mm, "assign")
  select_x = assign == args$i
  select_within1 = assign == (args$link[3, args$i])
  select_within2 = assign == (args$link[4, args$i])
  qr_mm = qr(mm)
  z1 = t(KhatriRao(X = t(mm_id1), Y = t(mm[, select_within1, drop = F])))
  z2 = t(KhatriRao(X = t(mm_id2), Y = t(mm[, select_within2, drop = F])))

  mm_id12 <- t(KhatriRao(X = t(mm_id1), Y = t(mm_id2)))

  if(sum(select_within1&select_within2)==0){
    z12 <- mm_id12
  }else{
    z12 <- t(KhatriRao(X = t(mm_id12), Y = t(mm[, select_within1&select_within2, drop = F])))}

  #z12 = qr.resid(qr_mm, z12)


  qr_d = qr(mm[, !select_x, drop = F])
  rdx = mm[, select_x, drop = F]


  qr_rdx = qr(rdx)
  qr_z12 = qr(z12)
  qr_z1 = qr(z1)
  qr_z2 = qr(z2)

  rank_x = sum(select_x)
  rank_z1 = ncol(z1)
  rank_z2 = ncol(z2)
  rank_z12 = ncol(z12)


  SSn1 = sum(qr.fitted(qr_rdx, args$y)^2)
  SSn2 = sum(qr.fitted(qr_z12, args$y)^2)
  SSd1 = sum(qr.fitted(qr_z1, args$y)^2)
  SSd2 = sum(qr.fitted(qr_z2, args$y)^2)

  MSn1 = SSn1/rank_x
  MSn2 = SSn2/rank_z12
  MSd1 = SSd1/rank_z1
  MSd2 = SSd2/rank_z2

  dfn = (MSn1+MSn2)^2/(MSn1^2/rank_x+MSn2^2/rank_z12)
  dfd = (MSd1+MSd2)^2/(MSd1^2/rank_z1+MSd2^2/rank_z2)

  qf = (MSn1+MSn2)/(MSd1+MSd2)

  out = c(SSn1 = SSn1, SSn2 = SSn2, SSd1 = SSd1, SSd2 = SSd2,
          dfn1 = rank_x,dfn2 = rank_z12,dfd1 = rank_z1, dfd2 = rank_z2,dfn = dfn ,dfd = dfd,
          `quasi F` = qf, `parametric P(>F)` = 1 -
            pf(q = qf, df1 = dfn, df2 = dfd))


  out
}
