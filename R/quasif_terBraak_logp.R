#'@importFrom stats pf
quasif_terBraak_logp = function(args){
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

  py = as.matrix(args$pry) + qr.fitted(qr_d,args$y)


  num1 = Matrix::colSums(Matrix::qr.fitted(qr_x, py)^2)/ rank_x
  num2 = Matrix::colSums(Matrix::qr.fitted(qr_z12, py)^2)/rank_z12
  den1 = Matrix::colSums(Matrix::qr.fitted(qr_z1, py)^2)/rank_z1
  den2 = Matrix::colSums(Matrix::qr.fitted(qr_z2, py)^2)/rank_z2



  num1[1] = sum(Matrix::qr.fitted(qr_x, args$y)^2)/ rank_x
  num2[1] = sum(Matrix::qr.fitted(qr_z12, args$y)^2)/rank_z12
  den1[1] = sum(Matrix::qr.fitted(qr_z1, args$y)^2)/rank_z1
  den2[1] = sum(Matrix::qr.fitted(qr_z2, args$y)^2)/rank_z2



  dfn = (num1+num2)^2/(num1^2/ rank_x+num2^2/rank_z12)
  dfd = (den1+den2)^2/(den1^2/rank_z1+den2^2/rank_z2)


  return(-pf(q = c((num1+num2)/(den1+den2)),dfn,dfd,log.p = T,lower.tail = F))
}
