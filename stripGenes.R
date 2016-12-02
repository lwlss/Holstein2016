library(qfa)

stripFile=function(fname,fout,strip,skip=13){
  factual=fname
  header=paste(readLines(fname,n=skip))
  body=read.delim(factual,skip=skip,stringsAsFactors=FALSE)
  body=body[!body$ORF%in%strip,]
  write.table(header,fout,quote=FALSE,row.names=FALSE,col.names=FALSE)
  write.table(body,"tmp.out",quote=FALSE,row.names=FALSE,sep="\t")
  file.append(fout,"tmp.out")
  file.remove("tmp.out")
}

flist=list.files("OriginalFitnessReports",pattern="*.txt",full.names=TRUE)
getQ=function(fname) return(strsplit(strsplit(fname,"/")[[1]][2],"_")[[1]][1])
qnums=sapply(flist,getQ)

commonStrip=c("YDR173C","YER069W","YHR018C","YJL071W","YJL088W","YML099C","YMR042W","YMR062C","YOL058W","YOL140W","YBR248C","YCL030C","YFR025C","YER055C","YIL020C","YIL116W","YCL018W","YGL009C","YHR002W","YLR451W","YNL104C","YOR108W","YBR115C","YDL131W","YDL182W","YDR034C","YDR234W","YGL154C",
"YIL094C","YIR034C","YNR050C","YMR038C")
names(commonStrip)=rep("Common",length(commonStrip))

sgd=readSGD()
neighbs=getNeighbours("exo1",kb=20,sgd,geneCoord="Start")

bgrnds=list()
bgrnds[["QFA0051"]]=c("cdc13","exo1")
bgrnds[["QFA0131"]]=c("rfa3")
bgrnds[["QFA0132"]]=c("lyp1")
bgrnds[["QFA0136"]]=c("stn1")
bgrnds[["QFA0139"]]=c("yku70")
bgrnds[["QFA0140"]]=c("cdc13")
bgrnds[["QFA0141"]]=c("ura3")
bgrnds[["QFA0142"]]=c("cdc13","rad9")

neighbs=list()
for(qnum in names(bgrnds)) neighbs[[qnum]]=getNeighbours(bgrnds[[qnum]],kb=20,sgd,geneCoord="Start")$FName

tostrip=list()
for(qnum in names(bgrnds)) tostrip[[qnum]]=c(neighbs[[qnum]],commonStrip)

newdir="GeneStripped"
tostrip=list()
for(qnum in names(bgrnds)) tostrip[[qnum]]=c(neighbs[[qnum]],commonStrip)
dir.create(newdir, showWarnings = FALSE)
for(f in flist){
  basef=strsplit(f,"/")[[1]][2]
  fout=file.path(newdir,basef)
  qnum=qnums[[f]]
  stripFile(f,fout,tostrip[[qnum]])  
}

newdir="CommonStripped"
tostrip=list()
for(qnum in names(bgrnds)) tostrip[[qnum]]=c(commonStrip)
dir.create(newdir, showWarnings = FALSE)
for(f in flist){
  basef=strsplit(f,"/")[[1]][2]
  fout=file.path(newdir,basef)
  qnum=qnums[[f]]
  stripFile(f,fout,tostrip[[qnum]])  
}

newdir="AllStripped"
tostrip=list()
for(qnum in names(bgrnds)) tostrip[[qnum]]=c(commonStrip,unlist(neighbs))
dir.create(newdir, showWarnings = FALSE)
for(f in flist){
  basef=strsplit(f,"/")[[1]][2]
  fout=file.path(newdir,basef)
  qnum=qnums[[f]]
  stripFile(f,fout,tostrip[[qnum]])  
}
