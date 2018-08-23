Zlist = function(formula,data){
  if(is.null(formula)){return(NULL)}
  ## set overparametrize model
  data = droplevels(data)
  for(i in 1:NCOL(data)){
    if(is.factor(data[,i])){
      if(length(levels(data[,i]))>=2 ){
        contrasts(data[,i],how.many = length(levels(data[,i]))) = contrasts(data[,i],contrasts=F)}else{
          data[,i]<-data[,i]
        }
    }else{
      data[,i]<-data[,i]
    }
  }
  #transform formula
  f_x = f_rv = formula
  f_rv[[3]] = f_rv[[2]] 
  f_rv[[2]] = NULL
  f_x[[2]] = NULL
  f_rv = update(f_rv,~.+0)
  
  rv_list = list(Matrix(model.matrix(f_rv,data),sparse=TRUE))
  #compute matrix
  mm_x = model.matrix(f_x,data)
  x_list = lapply(0:max(attr(mm_x,"assign")),function(ai){
    mmi = mm_x[,attr(mm_x,"assign")==ai, drop=F]
    Matrix(mmi,sparse=TRUE)
  })
  names(x_list) = c("(intercept)",attr(terms(f_x),"term.labels"))
  names(rv_list) = c(attr(terms(f_rv),"term.labels"))
  
  return(list(rv_list = rv_list, x_list = x_list))
}



# Zlist2 = function(formula,data){
#   if(is.null(formula)){return(NULL)}
#   ## set overparametrize model
#   data = droplevels(data)
#   for(i in 1:NCOL(data)){
#     if(is.factor(data[,i])){
#       if(length(levels(data[,i]))>=2 ){
#         contrasts(data[,i],how.many = length(levels(data[,i]))) = contrasts(data[,i],contrasts=F)}else{
#           data[,i]<-data[,i]
#         }
#     }else{
#       data[,i]<-data[,i]
#     }
#   }
#   #compute matrix
#   mm = model.matrix(formula,data)
#   zl = lapply(1:max(attr(mm,"assign")),function(ai){
#     mmi = mm[,attr(mm,"assign")==ai]
#     Matrix(mmi,sparse=TRUE)
#   })
#   names(zl) = attr(terms(formula),"term.labels")
#   return(zl)
# }


