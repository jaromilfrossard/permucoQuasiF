rm(list=ls())
library(permuco)
library(Matrix)
library(MASS)
library(combinat)
library(profvis)

devtools::install_github("jaromilfrossard/permucoQuasiF")
#devtools::install_github("jaromilfrossard/gANOVA")
Sys.setenv(LANG = "en")

dir= ""

#### simulation type
## decreasing variance through interaction
## increasing variable through interaction
## Large variance for id:item
##



t0 = proc.time()
#setwd(dir)


lf= list.files(paste(dir,"function_data",sep=""))

for(lfi in lf){
  print(lfi)
  source(paste(dir,"function_data/",lfi,sep=""))
}


ni= 18 #18
ns = 18#16
na =list(S = 3,I =3, C =2, SI = 1, O = 1, R = 1)


fseed = 42

source(paste("data_uni_lmm.R",sep=""))




f_s = subject~(fI*fC)
f_i = item~(fS*fC)
f_si = (subject_item)~(fC)

zl_s = Zlist(f_s,df)
scale_s = scalelist(Zlist = zl_s,type="flat")
omega_glist_s = Omega_glist(Zlist = zl_s,scale = scale_s, type = "homoscedastic")


zl_i = Zlist(f_i,df)
scale_i = scalelist(Zlist = zl_i,type="flat")
omega_glist_i = Omega_glist(Zlist = zl_i,scale = scale_i, type = "homoscedastic")


if(!is.null(f_si)){
  zl_si = Zlist(f_si,df)
  scale_si = scalelist(Zlist = zl_si,type="decreasing")
  omega_glist_si = Omega_glist(Zlist = zl_si,scale = scale_i, type = "homoscedastic")
  Omega = omega_glist_s$Omega + omega_glist_i$Omega + omega_glist_si$Omega+
    Diagonal(n = nrow(df))*(max(scale_s,scale_i,scale_si)*1.5)^2}else{
      Omega = omega_glist_s$Omega + omega_glist_i$Omega +
        Diagonal(n = nrow(df))*(max(scale_s,scale_i)*1.5)^2}

err = MASS:::mvrnorm(n=1, mu= rep(0,NCOL(Omega)), Sigma=Omega)
df$y=err

#matrixcalc:::is.positive.definite(as.matrix(Omega))


df =permuco:::changeContrast(df,contr = contr.sum)




fqf=y ~ fS * fC * fI + Error(subject/(fI*fC)) + Error(item/(fS*fC))



lf= list.files(paste(dir,"R",sep=""))

for(lfi in lf){
  print(lfi)
  source(paste(dir,"R/",lfi,sep=""))
}

qf_p = aovperm(fqf, data = df, np = 4000,method = "terBraak",effect=1)
qf_p = aovperm(fqf, data = df, np = 4000,method = "terBraak_logp",effect=1)

anova_table_quasif(arg)


arg$mm
rankMatrix(z12)

profvis({
  qf_p = aovperm(fqf, data = df, np = 2,method = "terBraak",effect=1)})





