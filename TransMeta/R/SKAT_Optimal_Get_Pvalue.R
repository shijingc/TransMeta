SKAT_Optimal_Get_Pvalue <-
function(Q.all, var, K, r.all, method){


n.r = length(r.all)
Q.all = t(as.matrix(as.numeric(Q.all)))
n.q<-dim(Q.all)[1]
p.m<-length(var)
Sigma.Hat = diag(var, p.m)

lambda.all<-list()
Sigma.neg.half = Sigma.Hat %^% (-0.5)

for(i in 1:n.r){
r.corr <- r.all[i]
R.M <- (1-r.corr) * K + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
K2 = Sigma.neg.half%*%R.M%*%Sigma.neg.half
lambda.all[[i]]<-SKAT:::Get_Lambda(K2) 

}
lambda.all1<-lambda.all

# Get Mixture param 
param.m <- SKAT_Optimal_Param(Sigma.Hat, K, r.all)
Each_Info <- SKAT_Optimal_Each_Q(param.m, Q.all, r.all, lambda.all, method=method)
pmin.q <- Each_Info$pmin.q
pmin <- Each_Info$pmin
pval <- rep(0,n.q)

if(method == "davies" || method=="optimal" || method=="optimal.mod" || method=="optimal.adj"){

for(i in 1:n.q){
pval[i]<-SKAT_Optimal_PValue_Davies(pmin.q[i,],param.m,r.all, pmin[i])

}


} else if(method =="liu" || method =="liu.mod" || method=="optimal.moment" || method=="optimal.moment.adj" ){

for(i in 1:n.q){
pval[i]<-SKAT_Optimal_PValue_Liu(pmin.q[i,],param.m,r.all, pmin[i])
}

} else {
stop("Invalid Method!")
}

# Check the pval 
# Since SKAT-O is between burden and SKAT, SKAT-O p-value should be <= min(p-values) * 2
# To correct conservatively, we use min(p-values) * 3

multi<-3
if(length(r.all) < 3){
multi<-2
}

for(i in 1:n.q){
pval.each<-Each_Info$pval[i,]
IDX<-which(pval.each > 0)

pval1<-min(pval.each) * multi
if(pval[i] <= 0 || length(IDX) < length(r.all)){
pval[i]<-pval1
}

# if pval==0, use nonzero min each.pval as p-value
if(pval[i] == 0){
if(length(IDX) > 0){
pval[i] = min(pval.each[IDX])
}
}

}

 
 return(list(p.value=pval,p.val.each=Each_Info$pval, opt.rho = r.all[which.min(Each_Info$pval)]))

}
