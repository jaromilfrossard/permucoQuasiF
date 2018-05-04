Zmat = function(mm0, mm, link, mm_id1, mm_id2, mm0_id1, mm0_id2, terms_f){


  mm_id12 = permuco:::khatrirao(a = mm_id1, b = mm_id2)
  attr(mm_id12,"contrasts") = attr(mm_id1,"contrasts")%x%attr(mm_id2,"contrasts")
  mm0_id12 = permuco:::khatrirao(a = mm0_id1, b = mm0_id2)


  assign = attr(mm, "assign")
  assign0 = attr(mm0, "assign")

  which_w1 = sort(unique(link[3, ]))
  which_w2 = sort(unique(link[4, ]))
  which_w12  = which_w1[which_w1%in%which_w2]


  zs1 <- lapply(which_w1,function(i){
    z = permuco:::khatrirao(a = mm_id1, b = mm[, assign==i, drop = F])
    z0 = permuco:::khatrirao(a = mm0_id1, b = mm0[, assign0==i, drop = F])
    if(i==0){
      cod = 1

    }else{
      ai = which(attr(terms_f,"factors")[,i] ==1)
      cod <- lapply(ai,function(aa)contr.sum(sum(assign==aa)+1))
      cod = Reduce(f = "%x%",cod[length(cod):1])

    }
    coding = attr(mm_id1,"contrasts")%x%cod
    return(list(z=z,z0=z0,coding = coding))
  })
  names(zs1)[which_w1!=0]= paste(attr(mm_id1,"label"),colnames(link)[which_w1[which_w1!=0]],sep=":")
  names(zs1)[which_w1==0] = attr(mm_id1,"label")



  zs2 <- lapply(which_w2,function(i){
    z = permuco:::khatrirao(a = mm_id2, b = mm[, assign==i, drop = F])
    z0 = permuco:::khatrirao(a = mm0_id2, b = mm0[, assign0==i, drop = F])
    if(i==0){
      cod = 1

    }else{
      ai = which(attr(terms_f,"factors")[,i] ==1)
      cod <- lapply(ai,function(aa)contr.sum(sum(assign==aa)+1))
      cod = Reduce(f = "%x%",cod[length(cod):1])

    }
    coding = attr(mm_id2,"contrasts")%x%cod
    return(list(z=z,z0=z0,coding = coding))
  })

  names(zs2)[which_w2!=0]= paste(attr(mm_id2,"label"),colnames(link)[which_w2[which_w2!=0]],sep=":")
  names(zs2)[which_w2==0] = attr(mm_id2,"label")

  zs12 <- lapply(which_w12,function(i){
    z = permuco:::khatrirao(a = mm_id12, b = mm[, assign==i, drop = F])
    z0 = permuco:::khatrirao(a = mm0_id12, b = mm0[, assign0==i, drop = F])
    if(i==0){
      cod = 1

    }else{
      ai = which(attr(terms_f,"factors")[,i] ==1)
      cod <- lapply(ai,function(aa)contr.sum(sum(assign==aa)+1))
      cod = Reduce(f = "%x%",cod[length(cod):1])

    }
    coding = attr(mm_id12,"contrasts")%x%cod
    return(list(z=z,z0=z0,coding = coding))
  })

  names(zs12)[which_w12!=0]= paste(paste(attr(mm_id1,"label"),attr(mm_id2,"label"),sep=":"),colnames(link)[which_w12[which_w12!=0]],sep=":")
  names(zs12)[which_w12==0] = paste(attr(mm_id1,"label"),attr(mm_id2,"label"),sep=":")


  zs <- c(zs1,zs2,zs12)

  z0 = Matrix(do.call("cbind",lapply(zs,function(x)x$z0)),sparse = T)
  z = Matrix(do.call("cbind",lapply(zs,function(x)x$z)),sparse = T)
  coding  = bdiag(lapply(zs,function(x)x$coding))


  qr_xz = qr(cbind(Matrix(mm,sparse = T),z))

  return(list(z0 = z0, z=z, coding = coding, zs = zs, qr_xz = qr_xz, zs=zs))

}
