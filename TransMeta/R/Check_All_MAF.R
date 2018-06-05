Check_All_MAF <-
function(object, n.study, t.snp){

count = 0
for(i in 1:n.study){
if( sum(is.na(object[[i]]$MAF)) == t.snp ){
count = count + 1
}
}
if(count == n.study){
Msg<-sprintf("No MAF data, cannot estimate F.st")
stop(Msg)
}
}
