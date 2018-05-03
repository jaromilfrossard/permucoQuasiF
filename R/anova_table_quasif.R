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
  mm = args$mm
  assign = attr(mm, "assign")
  select_x = assign == args$i
  select_within1 = assign == (args$link[3, args$i])
  select_within2 = assign == (args$link[4, args$i])
  qr_mm = qr(mm)
  z1 = permuco:::khatrirao(a = args$mm_id1, b = mm[, select_within1, drop = F])
  #z1 = qr.resid(qr_mm, z1)

  z2 = permuco:::khatrirao(a = args$mm_id2, b = mm[, select_within2, drop = F])
  #z2 = qr.resid(qr_mm, z2)

  z12 = permuco:::khatrirao(a = z1, b = z2)
  z12 = permuco:::khatrirao(b = z12,a = mm[,select_x,drop=F])

  #z12 = qr.resid(qr_mm, z12)


  qr_d = qr(mm[, !select_x, drop = F])
  rdx = qr.resid(qr_d, mm[, select_x, drop = F])
  rdx = mm[, select_x, drop = F]
  qr_rdx = qr(rdx)


  qr_z12 = qr(z12)
  qr_z1 = qr(z1)
  qr_z2 = qr(z2)


  SSn1 = sum(qr.fitted(qr_rdx, args$y)^2)
  SSn2 = sum(qr.fitted(qr_z12, args$y)^2)
  SSd1 = sum(qr.fitted(qr_z1, args$y)^2)
  SSd2 = sum(qr.fitted(qr_z2, args$y)^2)

  MSn1 = SSn1/qr_rdx$rank
  MSn2 = SSn2/qr_z12$rank
  MSd1 = SSd1/qr_z1$rank
  MSd2 = SSd2/qr_z2$rank

  dfn = (MSn1+MSn2)^2/(MSn1^2/qr_rdx$rank+MSn2^2/qr_z12$rank)
  dfd = (MSd1+MSd2)^2/(MSd1^2/qr_z1$rank+MSd2^2/qr_z2$rank)

  qf = (MSn1+MSn2)/(MSd1+MSd2)

  out = c(SSn1 = SSn1, SSn2 = SSn2, SSd1 = SSd1, SSd2 = SSd2,
          dfn1 = qr_rdx$rank,dfn2 = qr_z12$rank,dfd1 = qr_z1$rank, dfd2 = qr_z2$rank,dfn = dfn ,dfd = dfd,
          `quasi F` = qf, `parametric P(>F)` = 1 -
            pf(q = qf, df1 = dfn, df2 = dfd))
  out
}
