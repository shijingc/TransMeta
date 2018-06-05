SKAT_Optimal_PValue_Liu <-
function(pmin.q,param.m,r.all, pmin=NULL){

 re<-integrate(SKAT_Optimal_Integrate_Func_Liu, lower=0, upper=40, subdivisions=2000
,pmin.q=pmin.q,param.m=param.m,r.all=r.all,abs.tol = 10^-25)

pvalue<-1-re[[1]]

if(!is.null(pmin)){
if(pmin *length(r.all) < pvalue){
pvalue = pmin *length(r.all)
}
}

return(pvalue)

}
