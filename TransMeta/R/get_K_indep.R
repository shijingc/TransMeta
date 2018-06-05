get_K_indep <-
function(FileName){
data(G_ind)
Check_File_Exists(FileName[length(FileName)])
pop_id = read.table(FileName[length(FileName)])
colnames(pop_id)="Ancestry"

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

for ( m in 1:dim(pop_id)[1] ){
K.indep = NULL
for( n in 1:dim(pop_id)[1] ){
K.indep = c(K.indep, G_ind[ group[m], group[n] ])
}
K = rbind(K, K.indep)
}
K = as.matrix(K)
dimnames(K) <-list(rep("", dim(K)[1]), rep("", dim(K)[2]))

return(K)
}
