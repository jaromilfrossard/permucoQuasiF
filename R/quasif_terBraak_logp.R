quasif_terBraak_logp = function(args){
  link = args$link
  mm = args$mm
  assign = attr(mm, "assign")
  assign0 = attr(args$mm0, "assign")

  select_x = assign == args$i
  select_within1 = assign == (args$link[3, args$i])
  select_within2 = assign == (args$link[4, args$i])
  qr_mm = qr(mm)
  z1 = permuco:::khatrirao(a = args$mm_id1, b = mm[, select_within1, drop = F])

  z2 = permuco:::khatrirao(a = args$mm_id2, b = mm[, select_within2, drop = F])

  mm_id12 = permuco:::khatrirao(a = args$mm_id1, b = args$mm_id2)
  attr(mm_id12,"contrasts") = attr(args$mm_id1,"contrasts")%x%attr(args$mm_id2,"contrasts")


  z12 = permuco:::khatrirao(a = mm_id12, b = mm[,select_x,drop=F])


  qr_d = qr(mm[, !select_x, drop = F])
  qr_x = qr(mm[, select_x, drop = F])


  qr_z12 = qr(z12)
  qr_z1 = qr(z1)
  qr_z2 = qr(z2)


  #print(sry[,1]+mm%*%(qr.coef(qr_xz,args$y)[1:ncol(mm)])-args$y)

  py = as.matrix(args$pry) + qr.fitted(qr_d,args$y)
  num1 = colSums(qr.fitted(qr_x, py)^2)/qr_x$rank
  num2 = colSums(qr.fitted(qr_z12, py)^2)/qr_z12$rank
  den1 = colSums(qr.fitted(qr_z1, py)^2)/qr_z1$rank
  den2 = colSums(qr.fitted(qr_z2, py)^2)/qr_z2$rank


  num1[1] = sum(qr.fitted(qr_x, args$y)^2)/qr_x$rank
  num2[1] = sum(qr.fitted(qr_z12, args$y)^2)/qr_z12$rank
  den1[1] = sum(qr.fitted(qr_z1, args$y)^2)/qr_z1$rank
  den2[1] = sum(qr.fitted(qr_z2, args$y)^2)/qr_z2$rank



  dfn = (num1+num2)^2/(num1^2/qr_x$rank+num2^2/qr_z12$rank)
  dfd = (den1+den2)^2/(den1^2/qr_z1$rank+den2^2/qr_z2$rank)


  return(-pf(q = c((num1+num2)/(den1+den2)),dfn,dfd,log.p = T,lower.tail = F))
}
