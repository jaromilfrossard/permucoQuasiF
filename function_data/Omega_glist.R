# Zlist = zl_s;scale = scale_s;type = "correlated"
# Zlist = zl_i;scale = scale_i;type = "correlated"

Omega_glist = function(Zlist,scale, type = "homoscedastic"){
  
  z_rv = Zlist$rv_list[[1]]
  corr_bloc2 = 0.2
  corr_bloc = 0.5
  corr_max = 0.1
  
  
  # z0 = do.call("cBind",Zlist$x_list)
  
  pz = as.numeric(apply(z_rv,2,function(col)which(col==1)))
  temp = lapply(sapply(Zlist$x_list,function(x)ncol(x)),function(i){
    sort(rep(1:ncol(z_rv),i))
  })
  
  psigma = order(unlist(temp))

  
  # as.numeric(apply(z_rv,2,function(x)which(x==1)))
  # x1 = t(KhatriRao(t(do.call("cBind",Zlist$x_list)),t(z_rv)))
  # x = t(KhatriRao(t(z_rv),t(do.call("cBind",Zlist$x_list))))
  # id_order = as.numeric(apply(z_rv,2,function(x)which(x==1)))
  # 
  # x0 = t(KhatriRao(t(z_rv[id_order,]),t(do.call("cBind",Zlist$x_list)[id_order,])))
  

  switch(type,
         "homoscedastic" = {
           Sigma_glist = mapply(function(z_x,s){
             Diagonal(n = ncol(z_x),x = s^2 )}, 
             z_x = Zlist$x_list,s = scale, SIMPLIFY = F)
           Sigma = kronecker(Diagonal(n = ncol(z_rv)), bdiag(Sigma_glist))[psigma,psigma]},
         "heteroscedastic" = {
           Sigma_glist = mapply(function(z_x,s){
             heteroscale = seq(from= 0.75, to = 1.5,length.out = ncol(z_x))
             heteroscale = heteroscale-mean(heteroscale)+1
             Diagonal(n = ncol(z_x),x = heteroscale*s^2 )}, 
             z_x = Zlist$x_list,s = scale, SIMPLIFY = F)
           Sigma = kronecker(Diagonal(n = ncol(z_rv)), bdiag(Sigma_glist))[psigma,psigma]},
         "bloc" = {
           Sigma_glist = mapply(function(z_x,s){
             Diagonal(n = ncol(z_x),x = corr_bloc*s^2 )+(1-corr_bloc)*s^2}, 
             z_x = Zlist$x_list,s = scale, SIMPLIFY = F)
           Sigma = kronecker(Diagonal(n = ncol(z_rv)), bdiag(Sigma_glist))[psigma,psigma]},
         "correlated" = {
           Sigma_glist = mapply(function(z_x,s){
             Diagonal(n = ncol(z_x),x = corr_bloc*s^2 )+(1-corr_bloc)*s^2}, 
             z_x = Zlist$x_list,s = scale, SIMPLIFY = F)
           Sigma_g = bdiag(Sigma_glist)
           Sigma_g[1,-1] = corr_bloc2*scale[1]^2
           Sigma_g[-1,1] = corr_bloc2*scale[1]^2
           Sigma_g = Sigma_g + 0.5*corr_bloc2*scale[1]^2
           Sigma = kronecker(Diagonal(n = ncol(z_rv)), Sigma_g)[psigma,psigma]})
  
         

  ## Sigma = correlation for one id
  
  Z =  t(KhatriRao(t(z_rv),t(do.call("cBind",Zlist$x_list))))
  #product on reorder
  Omega  = ((Z[pz,])%*%(Sigma[order(psigma),order(psigma)])%*%t(Z[pz,]))[order(pz),order(pz)]

  return(list(Z = Z, Omega = Omega, Sigma = Sigma, psigma = psigma,pz = pz))
}

















# switch(type,
#        "homoscedastic" = {
#          Sigma_glist = mapply(function(z_x,s){
#            Sigma0 = Diagonal(n = ncol(z_x),x = s^2 )
#            kronecker(Diagonal(n = ncol(z_rv)), Sigma0)
#          }, z_x = Zlist$x_list,s = scale, SIMPLIFY = F)},
#        "heteroscedastic" = {
#          Sigma_glist = mapply(function(z_x,s){
#            heteroscale = seq(from= 0.75, to = 1.5,length.out = ncol(z_x))
#            heteroscale = heteroscale-mean(heteroscale)+1
#            Sigma0 = Diagonal(n = ncol(z_x),x = heteroscale*s^2 )
#            kronecker(Diagonal(n = ncol(z_rv)), Sigma0)
#          }, z_x = Zlist$x_list,s = scale, SIMPLIFY = F)},
#        "bloc" = {
#          Sigma_glist = mapply(function(z_x,s){
#            Sigma0 = Diagonal(n = ncol(z_x),x = corr_bloc*s^2 )+(1-corr_bloc)*s^2
#            kronecker(Diagonal(n = ncol(z_rv)), Sigma0)
#          }, z_x = Zlist$x_list,s = scale, SIMPLIFY = F)},
#        "maximal" = {
#          # Sigma_glist = mapply(function(z_x,s){
#          #   Sigma0 = Diagonal(n = ncol(z_x),x = corr_bloc*s^2 )+(1-corr_bloc)*s^2
#          #   kronecker(Diagonal(n = ncol(z_rv)), Sigma0)
#          # }, z_x = Zlist$x_list,s = scale, SIMPLIFY = F)
#          # Sigma_glist = list(bdiag(Sigma_glist)+corr_max)
#          # Zlist$x_list = list(do.call("cBind",Zlist$x_list))
#          
#        },
#        "correlated" = {
#          
#          Sigma_glist = mapply(function(z_x,s){
#            Sigma0 <- Diagonal(n = ncol(z_x),x = corr_bloc*s^2 )+(1-corr_bloc)*s^2
#            kronecker(Diagonal(n = ncol(z_rv)), Sigma0)
#          }, z_x = Zlist$x_list,s = scale, SIMPLIFY = F)
#          
#          Sigma_glist = list(bdiag(Sigma_glist))
#          corr = Matrix(0,ncol= sum(sapply(Zlist$x_list,ncol)),
#                        nrow= sum(sapply(Zlist$x_list,ncol)))
#          corr[1,-1] = corr_bloc2*scale[1]^2
#          corr[-1,1] = corr_bloc2*scale[1]^2
#          corr = kronecker(Diagonal(n = ncol(z_rv)), corr)
#          Sigma_glist[[1]] = Sigma_glist[[1]]+corr[id_order,id_order]
#          Zlist$x_list = list(do.call("cBind",Zlist$x_list))
#        })



