Calculate_TransMeta <-
function(i, study, n.study, forK, K.st, Beta, Var, r.all, method){ 

observed = study[which(forK[i,] != 0)]
k = K.st[observed, observed]
res = TransMeta(beta = Beta[i, observed], var = Var[i,observed], K =k, r.all = r.all, method = method) ##snp_id and p-value and opt.rho......
res1 = c( res$p.value,res$p.val.each, res$opt.rho )

return(res1)
}
