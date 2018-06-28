quasif_terBraak = function(args){
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


  mm_id12 <- permuco:::khatrirao(a = args$mm_id1, b = args$mm_id2)
  attr(mm_id12,"contrasts") = attr(args$mm_id1,"contrasts")%x%attr(args$mm_id2,"contrasts")
  z12 <- permuco:::khatrirao(a = mm_id12, b = mm[,select_within1&select_within2,drop=F])


  ##qr
  qr_d = qr(mm[, !select_x, drop = F])

  rank_x = sum(select_x)
  rank_z1 = ncol(z1)
  rank_z2 = ncol(z2)
  rank_z12 = ncol(z12)




  qr_x = Matrix::qr(Matrix((mm[, select_x, drop = F]),sparse = T))
  qr_z1 = Matrix::qr(Matrix(z1,sparse = T))
  qr_z2 = Matrix::qr(Matrix(z2,sparse = T))
  qr_z12 = Matrix::qr(Matrix(z12,sparse = T))

  py = as.matrix(args$pry) + qr.fitted(qr_d,args$y)

  #print(class(qr.fitted(qr_x, py)^2))

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

  return(c((num1+num2)/(den1+den2)))
}
