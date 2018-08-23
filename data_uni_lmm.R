
item = paste("I", sprintf("%02d",1:ni),sep="")
subject = paste("S", sprintf("%02d",1:ns),sep="")


fS = paste("s_",letters[c(1:na$S)],sep="")
fI = paste("i_",letters[c(1:na$I)],sep="")
fC = paste("c_",letters[c(1:na$C)],sep="")
fSI = paste("si_",letters[c(1:na$SI)],sep="")
fO = paste("o_",letters[c(1:na$O)],sep="")



# ### table subject
df_s = data.frame(subject=subject,fS=rep(fS,ceiling(ns/na$S))[1:ns])
# ### table item
df_i = data.frame(item=item,fI=rep(fI,ceiling(ni/na$I))[1:ni])
# ### table condition
df_c = data.frame(fC = fC)




# ### table subject
# df_s = data.frame(subject=subject,fS=rep(fS,ceiling(ns/na$S))[1:ns], xS = as.numeric(scale(rnorm(ns),scale = F)))
# ### table item
# df_i = data.frame(item=item,fI=rep(fI,ceiling(ni/na$I))[1:ni],xI = as.numeric(scale(rnorm(ni),scale = F)))
# ### table condition
# df_c = data.frame(fC = fC, xC =  as.numeric(scale(rnorm(na$C),scale = F)))


### table interaction
df = merge(df_s,df_i)
df_si = merge(rep(fSI,ceiling(ni/na$SI))[1:ni],rep(fSI,ceiling(ns/na$SI))[1:ns])

set.seed(fseed)
df$fSI = interaction(df_si$x,df_si$y)[sample(nrow(df))]

# inter = interaction(df$fS,df$fI)
# sel= inter == "s_a.i_a"
# df[sel,]

##permute level for full balancedment
levels(df$fSI) = do.call("c",permn(fSI))

# df$xSI = as.numeric(scale(rnorm(ns*ni),scale = F))

##add contition
df = merge(df, df_c)


### create Obsevational factor
stemp = df$subject
levels(stemp) = rep(fO,ceiling(ns/na$O))[1:ns]
itemp = df$item
levels(itemp) = rep(fO,ceiling(ni/na$O))[1:ni]
pfo = permn(fO)
pnfo = permn(length(pfo))


df = do.call("rbind",lapply(1:na$R,function(r){
  dfr = df
  dfr$replic = r
  dfr$fO = interaction(stemp,itemp)
  rmodulo = r%%length(pnfo)+1
  order = c(rmodulo:length(pnfo),0:(rmodulo-1))
  order = order[order>0]

  levelnames = lapply(1:length(pnfo),function(i){
    pi = (pnfo[order])[[i]]
    do.call("c",pfo[pi])
  })
  
  if(length(levelnames) == 1){
    levels(dfr$fO) = levelnames[[1]]}else{
    levels(dfr$fO) = do.call("c",levelnames)}
  return(dfr)}))


# df$xO = as.numeric(scale(rnorm(nrow(df)),scale = F))
df$subject_item = interaction(df$subject,df$item)
contrasts(df$subject_item) = contr.sum

fnames = paste("f",names(na),sep="")
fnames[6] = "replic"


df = df[,!(names(df) %in%(fnames[as.numeric(na)<=1]))]



#contr.scaled.helmert = function(n)scale(contr.helmert(n),center = F)
# my_contr = function(n)matrix(contr.poly(n),nrow=n)
# df = permuco:::changeContrast(df,contr = my_contr)
# 
# fX  = paste(c("~1",fnames[as.numeric(na)[-6]>1]),collapse = "+")
# 
# 
# df = data.frame(df, model.matrix(formula(fX),df)[,-1])
df = permuco:::changeContrast(df,contr = contr.sum)

rm(list = c("df_c", "df_i", "df_s", "df_si", 
            "fC","fI", "fO","fS" ,"fSI", "item","itemp", 
            "pfo","pnfo", "stemp", "subject"))







###### with formula


# 
# 
# # sl = as.factor(deg)
# # levels(sl) = 1/c(1:length(levels(sl)))
# # sl = as.numeric(as.character(sl))*2
# sl = rep(2,length(deg))
# sl[length(sl)]=max(sl)
# names(sl) = names(deg)
# sl = c(sl,epsilon =max(sl)*1.5)
# ## increasing variable through interaction
# # sl = 1/sl
# # sl = sl/max(sl)*2
# # sl[length(sl)]=max(sl)*1.5
# 
# ## random effect covariance matrice
# ### diagonal covariance matrix
# # sigmas = lapply(zl,function(z)diag(ncol(z)))
# # sigmas = mapply(function(sig,s)sig*s^2,sig=sigmas,s=sl,SIMPLIFY = F)
# 
# ### heretoscedastic random effect
# 
# sigmas = lapply(zl,function(z)Diagonal(n = ncol(z)))
# 
# 
# 
# 
# sigmas = mapply(function(sig,s)sig*s^2,sig=sigmas,s=sl,SIMPLIFY = F)
# 
# omegas = mapply(function(z,sig)Matrix(z%*%sig%*%t(z),sparse = T),z=zl,sig=sigmas,SIMPLIFY = F)
# 
# 
# 
# omega = Reduce("+",omegas)
# 
# 
# 
# 
# 
# 
# 
# 
# xtabs(~subject+fSI,)
# 
# 
# df$
# # xtabs(~subject+item,df)
# # 
# # xtabs(~subject+fI,df)
# # xtabs(~subject+fS,df)
# # xtabs(~subject+fSI,df)
# # xtabs(~subject+fC,df)
# # xtabs(~subject+fO,df)
# # 
# # xtabs(~item+fS,df)
# # xtabs(~item+fI,df)
# # xtabs(~item+fSI,df)
# # xtabs(~item+fC,df)
# # xtabs(~item+fO,df)
# # 
# # xtabs(~fI+interaction(subject,item),df)
# # xtabs(~fS+interaction(subject,item),df)
# # xtabs(~interaction(subject,item)+fSI,df)
# # xtabs(~fC+interaction(subject,item),df)
# # xtabs(~fO+interaction(subject,item),df)
# 
# # xtabs(~fSI+fC,df)
# # xtabs(~fSI+fO,df)
# # xtabs(~fO+replic,df)
# 
# 
# #xtabs(~subject+interaction(fI,fC,fSI),df)
# 
# 
# # f_s = ~subject/(fI*fC*fSI*fO)
# # f_i = ~item/(fS*fC*fSI*fO)
# # f_si = ~(subject*item)/(fC*fO)
# 
# 
# # f_s = ~subject/(fI*fC*fSI*fO*xI*xC*xSI*xO)
# # f_i = ~item/(fS*fC*fSI*fO*xS*xC*xSI*xO)
# # f_si = ~(subject*item)/(fC*fO*xC*xO)
# 
# 
# # f_s = ~subject/(fI*fC*fSI*xO)
# # f_i = ~item/(fC*fSI*xS*xO)
# # f_si = ~(subject*item)/(fC*xO)




