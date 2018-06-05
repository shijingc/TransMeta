get_K_Fst <-
function(FileName){ 
data(F_ST)
Check_File_Exists(FileName[length(FileName)])
pop_id = read.table(FileName[length(FileName)])
colnames(pop_id)="Ancestry"
n.pop = dim(pop_id)[1]

HMP3.pop = c("CEU", "CHD", "GIH", "JPT_CHB", "LWK", "MEX", "MKK", "TSI", "YRI")
K = group = NULL
for ( i in 1:dim(pop_id)[1]){
if( !(pop_id$Ancestry[i]  %in% HMP3.pop) ) {
stop('One or more population identifier is not identifiable!')
}
else{
for ( j in 1:9 ){
if ( as.character(pop_id$Ancestry[i]) == HMP3.pop[j] ) group = c(group, j)
}
}
}

D = 0
for ( p in 1:(n.pop-1) ){
for( q in (p+1):n.pop ){
D = max(D, F_ST[ group[p], group[q] ])
}
}


for ( m in 1:dim(pop_id)[1] ){
K.st = NULL
for( n in 1:dim(pop_id)[1] ){
K.st = c(K.st, 1- F_ST[ group[m], group[n] ]/D)
}
K = rbind(K, K.st)
}
K = as.matrix(K/max(K))
dimnames(K) <-list(rep("", dim(K)[1]), rep("", dim(K)[2]))

return(K)
}
