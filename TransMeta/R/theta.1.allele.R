theta.1.allele <-
function(maf.pool, n, pair){

maf = maf.pool[,pair]
n1 = n[,pair[1]]; n2 = n[,pair[2]]  ##? is n.all, the study size of all studies changing from snp to snp, if not, just make n.all 1*27 vector, instead of 2000*27 matrix
p = (n1*maf[,1] + n2*maf[,2])/(n1 + n2)
theta = ((maf[,1]-p)^2+(maf[,2]-p)^2)/(p*(1-p))

theta = sum(theta, na.rm = TRUE)/sum( !is.na( theta) ) 
return(theta)
}
