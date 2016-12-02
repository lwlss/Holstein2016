source("trimComplexes.R",local=TRUE)
fc=read.delim("FunctionalComplexesTrimmed.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)
fc=rbind(c("None","None","None",""),fc)
choiceList=as.list(1:length(fc$Notes))
names(choiceList)=fc$GroupName

# Read in all data once
dfs=list()
for(f in flist){
	dfnew=parseFitReport(f)
	dfs[[dfnew$ScrLabel[1]]]=dfnew
}

scrNms=names(dfs)
starting=1:length(scrNms)
#starting=c(23,26,27,28,30,31,32)
screens=scrNms[starting]

# Merge all datasets by ORF
fits=Reduce(function(...) {merge(...,by=c("ORF","Gene"),all=TRUE)}, dfs)
colroots=colnames(dfs[[1]])[3:length(colnames(dfs[[1]]))] 
cnames=c("ORF","Gene",c(t(sapply(colroots,paste,names(dfs),sep="."))))
colnames(fits)=cnames

# Eliminate some duplicated genes
fits=fits[!duplicated(fits$Gene),]

# Load fitness data consistent with default user input values
vals=list()
input=list()
input$dType="mahalanobis"
vals$prof=makeProfiles(fits,scrNms,"MeanFit")
vals$nExp=dim(vals$prof)[2]-2
vals$sims=similarities(vals$prof,"mahalanobis")
gdict=rownames(vals$prof)
names(gdict)=vals$prof$ORF

fc=read.delim("FunctionalComplexes.txt",sep="\t",header=TRUE,stringsAsFactors=FALSE)

# Discard rows where ORFs are not to be found in prof...
fcout=fc
reject=c()
for(i in 1:dim(fc)[1]){
	ORFs=strsplit(fc$GroupORFs[i]," ")[[1]]
	if(length(intersect(ORFs,names(gdict)))==0) reject=c(reject,i)
}
fc=fc[-1*reject,]
write.table(fc,"FunctionalComplexesTrimmed.txt",sep="\t",quote=FALSE,row.names=FALSE)

scrids=sapply(scrNms,function(x) {vec=fits[[paste("Screen_ID",x,sep=".")]]; noNAind=min(which(!is.na(vec))); vec[noNAind]})

#labs=colnames(vals$prof)[1:(length(colnames(vals$prof))-2)]
labs=paste(scrNms," (",scrids[scrNms],")",sep="")
chooses=1:length(labs)
names(chooses)=labs
