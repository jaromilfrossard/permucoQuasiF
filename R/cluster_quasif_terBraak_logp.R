#' @importMethodsFrom Matrix colMeans colSums
cluster_quasif_terBraak_logp = function(args){
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


  qr_xM = qr(Matrix((mm[, select_x, drop = F])[,qr_x$pivot[1:qr_x$rank]],sparse = T))
  qr_z1M = qr(Matrix(z1[,qr_z1$pivot[1:qr_z1$rank]],sparse = T))
  qr_z2M = qr(Matrix(z2[,qr_z2$pivot[1:qr_z2$rank]],sparse = T))
  qr_z12M = qr(Matrix(z12[,qr_z12$pivot[1:qr_z12$rank]],sparse = T))

  hdy = qr.fitted(qr_d,args$y)


  qf = t(apply(as.matrix(args$S),2,function(pii){
    pyi = args$zm$z0%*%(args$gamma*pii)+hdy
    num1 = colSums(qr.fitted(qr_xM, pyi)^2)/qr_x$rank
    num2 = colSums(qr.fitted(qr_z12M, pyi)^2)/qr_z12$rank
    den1 = colSums(qr.fitted(qr_z1M, pyi)^2)/qr_z1$rank
    den2 = colSums(qr.fitted(qr_z2M, pyi)^2)/qr_z2$rank
    dfn = (num1+num2)^2/(num1^2/qr_x$rank+num2^2/qr_z12$rank)
    dfd = (den1+den2)^2/(den1^2/qr_z1$rank+den2^2/qr_z2$rank)
    abs(pf(q = as.numeric(c((num1+num2)/(den1+den2))),dfn,dfd,log.p = T,lower.tail = F))
  }))



  num1 = colSums(qr.fitted(qr_xM, args$y)^2)/qr_x$rank
  num2 = colSums(qr.fitted(qr_z12M, args$y)^2)/qr_z12$rank
  den1 = colSums(qr.fitted(qr_z1M, args$y)^2)/qr_z1$rank
  den2 = colSums(qr.fitted(qr_z2M, args$y)^2)/qr_z2$rank

  dfn = (num1+num2)^2/(num1^2/qr_x$rank+num2^2/qr_z12$rank)
  dfd = (den1+den2)^2/(den1^2/qr_z1$rank+den2^2/qr_z2$rank)



  qf[1,] = abs(pf(q = c((num1+num2)/(den1+den2)),dfn,dfd,log.p = T,lower.tail = F))

  qf


}






