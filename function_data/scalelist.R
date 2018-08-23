#' compute the scales value from a list of random effect
#' @param Zlist a list of random effect
#' @param mx_sigma a scalar indicating the maximal value of sigma
#' @param type the type of structure,  choose between \code{"flat"}, \code{"decreasing"} and \code{"increasing"}
#' @return a vector a scale
scalelist = function(Zlist, mx_sigma = 2, type = "flat"){
  names = names(Zlist$x_list)
  neff = length(Zlist$x_list)
  degree = sapply(Zlist$x_list,function(zi)mean(stringr::str_count(colnames((zi)),":")))+1
  degree[names=="(intercept)"] = 0
  max_degree = max(max(degree),1)

  switch(type,
         "flat" = {scale = rep(mx_sigma, neff)},
         "increasing" = {scale = seq(from=1, to = max_degree+1,length.out = length(unique(degree)))[degree+1]/max_degree*mx_sigma},
         "decreasing" = {scale = seq(from=max_degree+1, to = 1,length.out = length(unique(degree)))[degree+1]/max_degree*mx_sigma})
  return(scale)
}

