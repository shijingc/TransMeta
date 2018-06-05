Get_wgt_maf <-
function(pos, maf, n){
if(pos[1] != pos[2]){
wgt.n = rowSums(n[ ,(pos[1]:pos[2])], na.rm = TRUE)
wgt.maf = rowSums(maf[ ,(pos[1]:pos[2])], na.rm = TRUE)/ rowSums(n[ ,(pos[1]:pos[2])], na.rm = TRUE)
}else{
wgt.n = n[ ,pos[1]]
wgt.maf = maf[ ,pos[1]]/n[ ,pos[1]]
}
res = cbind(wgt.n, wgt.maf)
return(res)

}
