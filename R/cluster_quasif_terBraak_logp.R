
cluster_quasif_terBraak_logp = function(args){
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

  ##qr
  qr_d = qr(mm[, !select_x, drop = F])

  rank_x = sum(select_x)
  rank_z1 = ncol(z1)
  rank_z2 = ncol(z2)
  rank_z12 = ncol(z12)





  qr_x = Matrix::qr(mm[, select_x, drop = F])
  qr_z1 = Matrix::qr(z1)
  qr_z2 = Matrix::qr(z2)
  qr_z12 = Matrix::qr(z12)


  hdy = qr.fitted(qr_d,args$y)



  qf = t(apply(as.matrix(args$S),2,function(pii){
    pyi <- args$zm$z0%*%(args$gamma*pii)+hdy
    num1 = Matrix::colSums(Matrix::qr.fitted(qr_x, pyi)^2)/rank_x
    num2 = Matrix::colSums(Matrix::qr.fitted(qr_z12, pyi)^2)/rank_z12
    den1 = Matrix::colSums(Matrix::qr.fitted(qr_z1, pyi)^2)/rank_z1
    den2 = Matrix::colSums(Matrix::qr.fitted(qr_z2, pyi)^2)/rank_z2
    dfn = (num1+num2)^2/(num1^2/rank_x +num2^2/rank_z12)
    dfd = (den1+den2)^2/(den1^2/rank_z1+den2^2/rank_z2)
    abs(pf(q = as.numeric(c((num1+num2)/(den1+den2))),dfn,dfd,log.p = T,lower.tail = F))
    # out = abs(pf(q = as.numeric(c((num1+num2)/(den1+den2))),dfn,dfd,log.p = T,lower.tail = F))
    # rm(num1,num2,den1,den2)
    # out
  }))



  num1 = Matrix::colSums(Matrix::qr.fitted(qr_x, args$y)^2)/rank_x
  num2 = Matrix::colSums(Matrix::qr.fitted(qr_z12, args$y)^2)/rank_z12
  den1 = Matrix::colSums(Matrix::qr.fitted(qr_z1, args$y)^2)/rank_z1
  den2 = Matrix::colSums(Matrix::qr.fitted(qr_z2, args$y)^2)/rank_z2



  dfn = (num1+num2)^2/(num1^2/rank_x +num2^2/rank_z12)
  dfd = (den1+den2)^2/(den1^2/rank_z1+den2^2/rank_z2)


  qf[1,] = abs(pf(q = c((num1+num2)/(den1+den2)),dfn,dfd,log.p = T,lower.tail = F))

  qf


}






