Check_MAF_Exists <-
function(i, object,FileName,t.snp){

if( sum(is.na(object$MAF)) == t.snp ){
Msg<-sprintf("File %s does not have valid MAF\n",FileName[i])
stop(Msg)
}

}
