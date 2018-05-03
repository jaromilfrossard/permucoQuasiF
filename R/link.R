link = function(formula_f,formula_within1,formula_within2){
  t_f = terms(formula_f)
  t_within1 = terms(formula_within1)
  t_within2 = terms(formula_within2)

  factors_f = attr(t_f, "factors")[-1, , drop = F]
  factors_between = colSums(factors_f[rownames(factors_f) %in%
                                        c(attr(t_within1, "term.labels"),attr(t_within2, "term.labels")), , drop = F]) == 0
  factors_link1 = factors_f
  factors_link1[!rownames(factors_f) %in% attr(t_within1, "term.labels"),
                ] = 0
  factors_link1 = apply(factors_link1, 2, function(l) {
    out = which(apply(factors_f, 2, function(f) identical(as.logical(f),
                                                          as.logical(l))))
    if (length(out) == 0) {
      out = 0
    }
    out
  })
  factors_link2 = factors_f
  factors_link2[!rownames(factors_f) %in% attr(t_within2, "term.labels"),
                ] = 0
  factors_link2 = apply(factors_link2, 2, function(l) {
    out = which(apply(factors_f, 2, function(f) identical(as.logical(f),
                                                          as.logical(l))))
    if (length(out) == 0) {
      out = 0
    }
    out
  })



  link = rbind(order = attr(t_f, "order"), between_effect = factors_between *
                 c(1:length(factors_between)), within1 = factors_link1,within2 = factors_link2)
  link}
