TransMeta <-
function(beta, var, K, r.all = c(0, 0.09, 0.25, 1), method = "davies"){
Q.all = SKAT_Optimal_Get_Q(beta, var, K, r.all)
pval = SKAT_Optimal_Get_Pvalue(Q.all, var, K, r.all, method)
return(pval)
}
