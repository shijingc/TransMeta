SKAT_Optimal_Param <-
function(Sigma.Hat, K, r.all){

p.m<-dim(Sigma.Hat)[1]
r.n<-length(r.all)

Sigma.neg.half = Sigma.Hat %^% (-0.5)
      Sigma.inv = Sigma.Hat %^% (-1)

Z = Sigma.neg.half %*% rep(1,p.m)
M = Z %*% solve(t(Z) %*% Z) %*% t(Z)

if(p.m == 1){
lambda = MuQ = VarQ = W3.3.item = 0
KerQ = 12
Df = 1
} else{
Z.item1 = Sigma.neg.half %*% (diag(1, p.m) - M)  
Z.item2 = Sigma.neg.half %*% M %*% Sigma.neg.half

# W3.2 Term : mixture chisq
W3.2.t <- t(Z.item1) %*% K %*% Z.item1
lambda <-  SKAT:::Get_Lambda(W3.2.t)

# W3.3 Term : variance of remaining ...
W3.3.item<-sum(diag(Z.item2 %*% K %*% Z.item1 %*% Sigma.neg.half %*% K)) * 4

# Mixture Parameters
MuQ <- sum(lambda)
VarQ <- sum(lambda^2) *2 + W3.3.item
KerQ <-  sum(lambda^4)/(sum(lambda^2))^2 * 12
Df <- 12/KerQ
}


# W3.1 Term : tau1 * chisq_1
tau = rep(0,r.n)
for(i in 1:r.n){
r.corr<-r.all[i]
a = solve(t(Z) %*% Z) 
b = t(Z) %*% Sigma.neg.half %*% K %*% Sigma.neg.half %*% Z
tau[i] <- (r.corr  + (a^2)* b * (1-r.corr)) /a
}

out<-list(MuQ=MuQ,VarQ=VarQ,KerQ=KerQ,lambda=lambda,VarRemain=W3.3.item,Df=Df,tau=tau)
return(out)
}
