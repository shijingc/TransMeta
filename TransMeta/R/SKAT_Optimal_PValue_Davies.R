SKAT_Optimal_PValue_Davies <-
function(pmin.q,param.m,r.all, pmin=NULL){

#re<-try(integrate(SKAT_Optimal_Integrate_Func_Davies, lower=0, upper=30, subdivisions=500, pmin.q=pmin.q,param.m=param.m,r.all=r.all,abs.tol = 10^-15), silent = TRUE)

re<-try(integrate(SKAT_Optimal_Integrate_Func_Davies, lower=0, upper=40, subdivisions=1000, pmin.q=pmin.q,param.m=param.m,r.all=r.all,abs.tol = 10^-25), silent = TRUE)

if(class(re) == "try-error"){
re<-SKAT_Optimal_PValue_Liu(pmin.q,param.m,r.all, pmin)
return(re)
} 

pvalue<-1-re[[1]]
if(!is.null(pmin)){
if(pmin *length(r.all) < pvalue){
pvalue = pmin *length(r.all)
}
}


return(pvalue)

}
