Get_TransMeta <-
function(FileName, type = "effect size", K.type = "Fst", K = NULL, r.all = c(0, 0.09, 0.25, 1), method = "davies"){

n.study = length(FileName) - 1



maf.all = n.all  = wgt.maf = wgt.n = NULL  ### maf.all is the 2000*27 study_size*maf's for all 27 studies, 
   ### n.all is the 2000*27 study sizes, may change it to 1*27 vector,
   ### wgt.maf is the 2000*9 weighted maf's for all 2000 snps from each of the 9 populations
Beta = Var = forK = NULL   ### beta is the 2000*27 matrix which stores all the beta.hat (missing value is denoted as NA),
  ### var is the 2000*27 matrix which stores all the se^2,
  ###forK 2000*27 matrix, each row records all the studie # that observe the correspongding snp, 
  ### if a study does not observe the snp, then the value is 0 in the forK matrix

### read in Meta Datasets
objects <- list()
for (i in 1:n.study) {
Check_File_Exists(FileName[i])
objects[[i]] <- read.table(FileName[i],header = TRUE)
colnames(objects[[i]]) = c("SNP", "N", "Test_Allele", "Alt_Allele", "MAF", "Beta", "SE", "P_Val")
if(i == 1){
obj = data.frame(objects[[i]]$SNP, objects[[i]]$Test_Allele )
colnames(obj) = c("SNP","Test_Allele") 
}else{
obj.mrg = data.frame(objects[[i]]$SNP, objects[[i]]$Test_Allele  )
colnames(obj.mrg) = c("SNP","Test_Allele")
obj = merge(obj, obj.mrg, by=c("SNP","Test_Allele"), all=TRUE)
obj = obj[!duplicated(obj$SNP),]

}
}  ### obj stores all observed snps and test_allele across all studies
all.snp = as.data.frame(obj$SNP);colnames(all.snp) = "SNP"

t.snp = dim(obj)[1] ## total number of snps across all studies
for (i in 1:n.study) {
objects[[i]] = merge(all.snp, objects[[i]],by ="SNP", all=TRUE)
objects[[i]] = objects[[i]][order(objects[[i]]$SNP),]
Missing_Ind = !is.na(objects[[i]]$Beta)
objects[[i]] = cbind(objects[[i]], Missing_Ind)
maf.all = cbind(maf.all, objects[[i]]$Missing_Ind * objects[[i]]$N * objects[[i]]$MAF)
n.all = cbind(n.all, objects[[i]]$Missing_Ind * objects[[i]]$N ) 
forK = cbind(forK, i*Missing_Ind)

dir = 2*(obj$Test_Allele == objects[[i]]$Test_Allele)-1

if(type == "effect size"){
Beta = cbind(Beta, dir*objects[[i]]$Beta)
Var = cbind(Var, objects[[i]]$SE^2)
}else{
Check_MAF_Exists(i, objects[[i]], FileName, t.snp)
z = abs(qnorm(objects[[i]]$P_Val/2))*sign(dir*objects[[i]]$Beta)
tranf = objects[[i]]$N * objects[[i]]$MAF *(1-objects[[i]]$MAF)
Beta = cbind(Beta, z/sqrt(tranf))
Var = cbind(Var, 1/tranf)
}
}


### read in Ancestry Group data
Check_File_Exists(FileName[n.study+1])
pop_id = read.table(FileName[n.study+1])
colnames(pop_id)="Ancestry"

if(is.null(K)){
if(K.type == "Fst"){
Check_All_MAF(objects, n.study, t.snp)

### 'n.pop.rep' stores the number of studies for each ancestry group
pop = as.vector(unique(pop_id));n.pop = dim(pop)[1] 
n.pop.rep = NULL
for(t in 1:n.pop){
n.pop.rep = c(n.pop.rep, length(which(pop_id[,1] == pop[t,1])))  
}

##### 'position' stores the starting and ending positions for each ancestry group
start = c(1, cumsum(n.pop.rep[-length(n.pop.rep)])+1)
end = cumsum(n.pop.rep)
position = cbind(start, end)

for(t in 1:n.pop){ ### wgt.maf is the 2000*9 weighted maf's for all 2000 snps from each of the 9 populations
wgt = Get_wgt_maf( pos = position[t,], maf = maf.all, n = n.all)
wgt.maf = cbind(wgt.maf, wgt[,2])
wgt.n = cbind(wgt.n, wgt[,1])
}

F.st = matrix(0, nrow = n.pop, ncol = n.pop) ##F.st is the 9*9 population K matrix
for(j in 1:(n.pop-1)){
for(k in (j+1):n.pop){ 
F.st[j,k] = F.st[k,j] =  theta.1.allele(wgt.maf, wgt.n, c(j,k))
}
}
F.st = 1 - F.st/max(F.st)

K.st = NULL ## K.st is the 27*27 study K matrix
for(t in 1:n.pop){
fst = rep(F.st[t,], n.pop.rep)  ## one row of the 27*27 study pair-wise F.st values
m = n.pop.rep[t]
n = 1
while( n <= m){
K.st = rbind(K.st, fst)
n = n + 1
}
}
}else{
K.st = get_K_indep(FileName)
}
}else{
K.st = K
}

study = seq(1:n.study)
result = t(apply(as.matrix(seq(1:t.snp)), 1, 'Calculate_TransMeta', study = study, n.study = n.study, forK = forK, K.st = K.st, Beta = Beta, Var = Var, r.all = r.all, method = method))
result = data.frame(SNP = objects[[1]]$SNP, result)

s = c("SNP", "Pval")
for(i in 1:length(r.all)){
s.rho = paste("Pval_rho",i, sep="")
s = c(s, s.rho)
}
s = c(s, "Opt_rho")
colnames(result) = s

return(result)
}
