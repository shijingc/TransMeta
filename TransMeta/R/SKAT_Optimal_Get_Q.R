SKAT_Optimal_Get_Q <-
function(beta, var, K, r.all, n.Resampling = 0, res.moments=NULL, Q.sim=NULL){

n.r <- length(r.all)
p.m <- length(beta)

Beta.Hat = as.numeric(beta)
Sigma.Hat = diag(var, p.m)


Q.r <- rep(0,n.r)
Q.r.res <- NULL

Sigma.inv = Sigma.Hat %^% (-1)
one = rep(1, p.m)

for(i in 1:n.r){
r.corr <- r.all[i]
Q1 <- (1-r.corr) * t(Beta.Hat) %*% Sigma.inv %*% K %*% Sigma.inv %*% Beta.Hat
Q2 <- r.corr * t(Beta.Hat) %*% Sigma.inv %*% one %*% t(one) %*% Sigma.inv %*% Beta.Hat
Q.r[i] <- Q1 + Q2
}


  if(n.Resampling > 0){

#temp<-t(res.out) %*% Z1
Q.r.res<-matrix(rep(0,n.Resampling *n.r),ncol=n.r)
for(i in 1:n.r){
r.corr<-r.all[i]
Q1<-(1-r.corr) * t(Beta.Hat) %*% Sigma.inv %*% K %*% Sigma.inv %*% Beta.Hat
Q2<-r.corr * t(Beta.Hat) %*% Sigma.inv %*% one %*% t(one) %*% Sigma.inv %*% Beta.Hat
Q.r[i]<-Q1 + Q2
}

  }

#if(!is.null(res.moments) && is.null(Q.sim)){

#temp<-t(res.moments) %*% Z1
#n.moments<-dim(res.moments)[2]
#Q.sim<-matrix(rep(0,n.moments *n.r),ncol=n.r)
#for(i in 1:n.r){
#r.corr<-r.all[i]
#Q1<-(1-r.corr) * rowSums(temp^2)
#Q2<-r.corr * p.m^2 * rowMeans(temp)^2
#Q.sim[,i]<-Q1 + Q2
#}

#}

re = t(as.matrix(Q.r))
#re = as.numeric(c(Q.r, input[(p.m+1):(2*p.m)]))

return(re)


}
