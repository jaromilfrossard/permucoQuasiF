
cluster_quasif_terBraak_logp = function(args){
  link = args$link
  mm = args$mm
  assign = attr(mm, "assign")
  assign0 = attr(args$mm0, "assign")

  select_x = assign == args$i
  select_within1 = assign == (args$link[3, args$i])
  select_within2 = assign == (args$link[4, args$i])
  qr_mm = qr(mm)
  z1 <- permuco:::khatrirao(a = args$mm_id1, b = mm[, select_within1, drop = F])

  z2 <- permuco:::khatrirao(a = args$mm_id2, b = mm[, select_within2, drop = F])

  mm_id12 <- permuco:::khatrirao(a = args$mm_id1, b = args$mm_id2)
  attr(mm_id12,"contrasts") = attr(args$mm_id1,"contrasts")%x%attr(args$mm_id2,"contrasts")
  z12 <- permuco:::khatrirao(a = mm_id12, b = mm[,select_within1&select_within2,drop=F])



  ##base qr for rank and pivot
  qr_d = qr(mm[, !select_x, drop = F])

  rank_x = sum(select_x)
  rank_z1 = ncol(z1)
  rank_z2 = ncol(z2)
  rank_z12 = ncol(z12)




  qr_xM = Matrix::qr(Matrix((mm[, select_x, drop = F]),sparse = T))
  qr_z1M = Matrix::qr(Matrix(z1,sparse = T))
  qr_z2M = Matrix::qr(Matrix(z2,sparse = T))
  qr_z12M = Matrix::qr(Matrix(z12,sparse = T))

  hdy = qr.fitted(qr_d,args$y)



  qf = t(apply(as.matrix(args$S),2,function(pii){
    pyi <- args$zm$z0%*%(args$gamma*pii)+hdy
    num1 = Matrix::colSums(Matrix::qr.fitted(qr_xM, pyi)^2)/rank_x
    num2 = Matrix::colSums(Matrix::qr.fitted(qr_z12M, pyi)^2)/rank_z12
    den1 = Matrix::colSums(Matrix::qr.fitted(qr_z1M, pyi)^2)/rank_z1
    den2 = Matrix::colSums(Matrix::qr.fitted(qr_z2M, pyi)^2)/rank_z2
    dfn = (num1+num2)^2/(num1^2/rank_x +num2^2/rank_z12)
    dfd = (den1+den2)^2/(den1^2/rank_z1+den2^2/rank_z2)
    abs(pf(q = as.numeric(c((num1+num2)/(den1+den2))),dfn,dfd,log.p = T,lower.tail = F))
    # out = abs(pf(q = as.numeric(c((num1+num2)/(den1+den2))),dfn,dfd,log.p = T,lower.tail = F))
    # rm(num1,num2,den1,den2)
    # out
  }))



  num1 = Matrix::colSums(Matrix::qr.fitted(qr_xM, args$y)^2)/rank_x
  num2 = Matrix::colSums(Matrix::qr.fitted(qr_z12M, args$y)^2)/rank_z12
  den1 = Matrix::colSums(Matrix::qr.fitted(qr_z1M, args$y)^2)/rank_z1
  den2 = Matrix::colSums(Matrix::qr.fitted(qr_z2M, args$y)^2)/rank_z2



  dfn = (num1+num2)^2/(num1^2/rank_x +num2^2/rank_z12)
  dfd = (den1+den2)^2/(den1^2/rank_z1+den2^2/rank_z2)


  qf[1,] = abs(pf(q = c((num1+num2)/(den1+den2)),dfn,dfd,log.p = T,lower.tail = F))

  qf


}






