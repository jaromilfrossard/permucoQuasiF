rm(list=ls())
library(lme4)
library(MASS)
library(devtools)
library(Matrix)
library(abind)
#install_github("jaromilfrossard/permuco")
#install_github("jaromilfrossard/permucoQuasiF")
library(permuco)
library(permucoQuasiF)

data("signal")
data("design")



lf = list.files("R")
for(i in 1:length(lf)){
  source(paste("R/", lf,sep="")[i])
}


modelQuasif <- clusterlm(signal ~A*B*C + Error(id/(B*C))+ Error(item/(A*C)),data = design,np=3)


library(profvis)
Sys.setenv(LANG = "en")


library()




# save(signal,file="signal_18i_20s.rda")
# save(design,file="design_18i_20s.rda")

dir=  "../../article/quasif/article_quasif/2017_09_15_cluster_quasif"
#dir="../../../Dropbox/Uni/article/quasif/article_quasif/2017_09_15_cluster_quasif"

lf = list.files(paste(dir, "/function_data_signal/",sep=""))
for(i in 1:length(lf)){
  source(paste(dir, "/function_data_signal/", lf,sep="")[i])
}







ni= 18#18 #
ns = 20 #20
t = 10
na =list(A=2,B=3,C=2)
#aggr_FUN = function(x){if(length(x)==0){return(-Inf)}else{sum(log(1-x))}}
aggr_FUN = sum
np = 4000


fseed = 42000


set.seed(fseed)
source(paste(dir, "/data_fix_signal_quasif.R",sep=""))
err = rand_n(zl = zl,corl = corl,sl = sl)
y = x%*%beta+err

#df$y=y[,1]
fqf = y~A*B*C + Error(id/(B*C))+ Error(item/(A*C))



dir=  "R"


lf = list.files(paste(dir,sep=""))
for(i in 1:length(lf)){
  source(paste(dir,"/", lf,sep="")[i])
}


np=400


qf_p = clusterlm(fqf,df,np=np,method = "terBraak_logp",aggr_FUN = aggr_FUN,threshold = abs(log(1-0.95)),return_distribution=T,effect=c(1))



profvis({
  qf_p = clusterlm(fqf,df,np=np,method = "terBraak_logp",aggr_FUN = aggr_FUN,threshold = abs(log(1-0.95)),return_distribution=T,effect=c(1))})



Rprof(tf <- "rprof.log", memory.profiling=TRUE,interval = 0.002)

qf_p = clusterlm(fqf,df,np=np,method = "terBraak_logp",aggr_FUN = aggr_FUN,threshold = abs(log(1-0.95)),return_distribution=T,effect=c(1))

Rprof(NULL)
summaryRprof(tf)




profvis({
  qf_p = clusterlm(fqf,df,np=np,method = "terBraak_logp",aggr_FUN = aggr_FUN,threshold = abs(log(1-0.95)),return_distribution=T,effect=c(1))})



qf_p = clusterlm(fqf,df,np=np,method = "terBraak_logp",aggr_FUN = aggr_FUN,threshold = abs(log(1-0.95)),return_distribution=T,effect=1)









qf_p = aovperm(fqf,df,np=np,method = "terBraak_logp",aggr_FUN = aggr_FUN,effect=c(3,2))
qf_p = aovperm(fqf,df,np=np,method = "terBraak_logp",aggr_FUN = aggr_FUN,effect=c(2,3))

qf_p$distribution



qf_p2 = permucoQuasiF::clusterlm(fqf,df,np=np,method = "terBraak_logp",aggr_FUN = aggr_FUN,threshold = abs(log(1-0.95)))

exp(qf_p2$multiple_comparison$A$uncorrected$main[,1])
plot(qf_p2)
plot(qf_p)


df2 =df
df2$y=y[,1]
aovp= permucoQuasiF::aovperm(fqf,df2,np=np,method = "terBraak_logp",aggr_FUN = aggr_FUN,threshold = abs(log(1-0.95)))
aovp= aovperm(fqf,df2,np=np,method = "terBraak_logp",aggr_FUN = aggr_FUN,threshold = abs(log(1-0.95)))
aovp= aovperm(fqf,df2,np=np,method = "terBraak",aggr_FUN = aggr_FUN,threshold = abs(log(1-0.95)))


plot(qf_p2$multiple_comparison$A$uncorrected$main[,1],
     qf_p$multiple_comparison$A$uncorrected$main[,1])


qf_p$multiple_comparison$A$uncorrected

dim(z1)
dim(z2)


ag$


qf_p$model.matrix_id1

dim(z12)
dim(z12bis)


np=4000
profvis({
qf_p = clusterlm(fqf,df,np=np,method = "terBraak_logp",aggr_FUN = aggr_FUN,threshold = abs(log(1-0.95)))
})







