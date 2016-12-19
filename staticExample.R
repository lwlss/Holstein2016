#source("profilingFunctions.R")
#source("global.R")

# Plot fitness profiles for most similar and most highly correlated genes
plotSimilarFit=function(targ,prof,nearest=c(),farthest=c(),ylim=c(),legendval=NULL,cexlabs=1.5,mainadd="",showgenes=TRUE){
	# Number of experiments
	nExp=dim(prof)[2]-2

	# Point characteristics
	df=data.frame()
	for(trg in names(targ)){
		N=length(targ[[trg]]$genes)

		dftrg=data.frame(
					genes=targ[[trg]]$genes,
					cols=targ[[trg]]$cols[((1:N)-1)%%length(targ[[trg]]$cols)+1],
					pchs=targ[[trg]]$pchs[((1:N)-1)%%length(targ[[trg]]$pchs)+1],
					ltys=targ[[trg]]$ltys[((1:N)-1)%%length(targ[[trg]]$ltys)+1],
					cexs=targ[[trg]]$cexs[((1:N)-1)%%length(targ[[trg]]$cexs)+1],
					types=targ[[trg]]$types[((1:N)-1)%%length(targ[[trg]]$types)+1],
					lwds=targ[[trg]]$lwds[((1:N)-1)%%length(targ[[trg]]$lwds)+1],
					group=trg,
					stringsAsFactors=FALSE
		)
		df=rbind(df,dftrg)
	}
	
	if(length(ylim)==0) ylim=quantile(as.numeric(unlist(prof[,1:nExp])),c(0,1),na.rm=TRUE)
	#ylim=pmin(ylim,150)
	# Visualise similarity of similar genes
	nearprof=data.frame(prof[rownames(prof)%in%names(nearest),1:nExp])
	farprof=data.frame(prof[rownames(prof)%in%names(farthest),1:nExp])
	# For some reason, if prof only contains one experiment, rownames are obliterated, so replace them:
	rownames(nearprof)=rownames(prof)[rownames(prof)%in%names(nearest)]
	rownames(farprof)=rownames(prof)[rownames(prof)%in%names(farthest)]
	if((length(targ)==1)&(length(targ[1]$genes)==1)) {
		mlab=paste(mainadd,tolower(targ[1]$genes),"profile")
		cols=c("red")
	}else{
		mlab=paste(mainadd,"Fitness profiling")
		cols=rainbow(length(targ))
	}
	
	# Similarity
	#if(is.null(legendval)){
	  op=par(mai=c(2,1,1,1))
	#}else{
	#  op=par(mai=c(4,2,1,1),mar=par()$mar+c(0,0,0,3))
	#}
	profdat=prof[,1:nExp]
	delta=0.5
	bcol="gray92"
	pcol=rgb(1,0.65,0.65)
	#plot(NULL,xlab="",ylab="Fitness",main=mlab,ylim=ylim,cex.lab=2,xlim=c(1-delta,nExp+delta),cex.main=2,xaxt="n",yaxt="n", mgp=c(5, 2, 0))
	maxvals=apply(profdat,2,max,na.rm=TRUE); minvals=apply(profdat,2,min,na.rm=TRUE)
	#polygon(c(1-delta,1:nExp,nExp+delta,nExp+delta,rev(1:nExp),1-delta),c(maxvals[1],maxvals,maxvals[length(maxvals)],minvals[length(minvals)],rev(minvals),minvals[1]),col=bcol,border=bcol)
	#boxplot(profdat,col="lightblue",notch=FALSE,outline=FALSE,border="black",las=2,range=1.5,cex.axis=2,add=TRUE,lwd=1.5,pars = list(boxwex = 0.45, staplewex = 0.5, outwex = 0.5))

	stripchart(profdat,vertical=TRUE,method="jitter",jitter=0.25,pch=16,col=rgb(0,0,0,0.1),cex=0.75,las=2,main=mlab,ylim=ylim,cex.lab=2,xlim=c(1-delta,nExp+delta),cex.main=2,cex.axis=2,ylab="Fitness")
	
	#points(maxvals,pch="-",col="darkgrey",cex=3)
	#points(minvals,pch="-",col="darkgrey",cex=3)
	
	pnk="pink"
	pnk=rgb(1.0,0.5,0.5)
	pnk=rgb(0.1176471,0.5647059,1.0000000,0.5)
	
	if(length(nearest)>0) for(i in 1:dim(nearprof)[1]) points(1:nExp,nearprof[i,],type="b",col=pnk,lwd=2)
	if(length(farthest)>0) for(i in 1:dim(farprof)[1]) points(1:nExp,farprof[i,],type="b",col="grey",lwd=2)
	for(trg in names(targ)) {
		for(i in seq_along(targ[[trg]]$genes)) {
			typeorig=targ[[trg]]$types[(i-1)%%length(targ[[trg]]$types)+1]
			types=typeorig
			if(types=="p") {xvals=1:nExp+runif(nExp,-0.2,0.2)}else{xvals=1:nExp}
			types[types=="b"]="p"
			types[types=="c"]="n"
			
			points(
				x=xvals,
				y=prof[targ[[trg]]$genes[i],1:nExp],
				col=targ[[trg]]$cols[(i-1)%%length(targ[[trg]]$cols)+1],
				type=types,
				lwd=targ[[trg]]$lwds[(i-1)%%length(targ[[trg]]$lwds)+1],
				cex=targ[[trg]]$cexs[(i-1)%%length(targ[[trg]]$cexs)+1],
				pch=targ[[trg]]$pchs[(i-1)%%length(targ[[trg]]$pchs)+1],
				lty=targ[[trg]]$ltys[(i-1)%%length(targ[[trg]]$ltys)+1]
			)
			types=typeorig
			types[types=="p"]="n"
			types[types=="b"]="c"
			points(
				x=1:nExp,
				y=prof[targ[[trg]]$genes[i],1:nExp],
				col=adjustcolor(targ[[trg]]$cols[(i-1)%%length(targ[[trg]]$cols)+1],alpha.f=0.25),
				type=types,
				lwd=targ[[trg]]$lwds[(i-1)%%length(targ[[trg]]$lwds)+1],
				cex=targ[[trg]]$cexs[(i-1)%%length(targ[[trg]]$cexs)+1],
				pch=targ[[trg]]$pchs[(i-1)%%length(targ[[trg]]$pchs)+1],
				lty=targ[[trg]]$ltys[(i-1)%%length(targ[[trg]]$ltys)+1]
			)
		}
	}
	if(showgenes){
	  if(dim(nearprof)[1]>0){
		text(nExp,nearprof[,nExp],tolower(rownames(nearprof)),cex=cexlabs,pos=4,col=pnk)
		text(1,nearprof[,1],tolower(rownames(nearprof)),cex=cexlabs,pos=2,col=pnk)
	  }
	  if(dim(farprof)[1]>0){
		text(nExp,farprof[,nExp],tolower(rownames(farprof)),cex=cexlabs,pos=4,col="grey")
		text(1,farprof[,1],tolower(rownames(farprof)),cex=cexlabs,pos=2,col="grey")
	  }
	  for(trg in names(targ)){
		for(i in seq_along(targ[[trg]]$genes)){
			text(nExp,prof[targ[[trg]]$genes[i],nExp],tolower(targ[[trg]]$genes[i]),cex=cexlabs,pos=4,col=targ[[trg]]$cols[(i-1)%%length(targ[[trg]]$cols)+1])
			text(1,prof[targ[[trg]]$genes[i],1],tolower(targ[[trg]]$genes[i]),cex=cexlabs,pos=2,col=targ[[trg]]$cols[(i-1)%%length(targ[[trg]]$cols)+1])
		}
	  }
	}
	if(!is.null(legendval)){
		plot.new()
		legend(x=legendval,legend=tolower(df$genes),col=df$cols,cex=mean(df$cexs),pch=df$pchs,inset=c(-0.3,0),xpd=TRUE)
		#legend(x=legendval,legend=tolower(df$genes),col=df$cols,cex=mean(df$cexs),pch=df$pchs)
	}
	par(op)
}


format=list(
NARROW=list(width=4.75,height=10,cexlabs=1.6,fname="teloplots%02d.pdf",pfname="teloplotsLEGEND%02d.pdf",pcexlabs=0.8),
WIDE=list(width=10,height=7,cexlabs=1.6,fname="Wideteloplots%02d.pdf",pfname="WideteloplotsLEGEND%02d.pdf",pcexlabs=0.8)
)

for(f in names(format)){

profALL=vals$prof
colnames(profALL)[1:length(prettyNames)]=prettyNames
prof=profALL[,c(1,5,8,11,12,13)]
lwd=1.25
cex=2.0
type="b"
lty=1
pch=1

targ_all=list(
CHK=list(genes=c("DDC1","RAD24","RAD17"),cols="darkgreen",pchs=pch,ltys=lty,lwds=lwd,cexs=cex,types=type),
NMD=list(genes=c("NAM7","NMD2","UPF3"),cols="dodgerblue",pchs=pch,ltys=lty,lwds=lwd,cexs=cex,types=type),
KU=list(genes=c("YKU70","YKU80"),cols="black",pchs=pch,ltys=lty,lwds=lwd,cexs=cex,types=type),
MRX=list(genes=c("MRE11","XRS2","RAD50"),cols="red",pchs=pch,ltys=lty,lwds=lwd,cexs=cex,types=type),
TELO=list(genes=c("RAD9","CHK1","TEL1","EXO1"),cols=c("purple","orange","darkcyan","magenta"),pchs=pch,ltys=lty,lwds=lwd,cexs=cex,types=type)
)


cairo_pdf(file=format[[f]]$fname,height=format[[f]]$height,width=format[[f]]$width,pointsize=6)
plotSimilarFit(targ_all,prof,cexlabs=format[[f]]$cexlabs)
targ=list(CHK=targ_all$CHK)
plotSimilarFit(targ,prof,cexlabs=format[[f]]$cexlabs,mainadd="CHK")
targ=list(NMD=targ_all$NMD)
plotSimilarFit(targ,prof,cexlabs=format[[f]]$cexlabs,mainadd="NMD")
targ=list(KU=targ_all$KU)
plotSimilarFit(targ,prof,cexlabs=format[[f]]$cexlabs,mainadd="KU")
targ=list(MRX=targ_all$MRX)
plotSimilarFit(targ,prof,cexlabs=format[[f]]$cexlabs,mainadd="MRX")
targ=list(TELO=targ_all$TELO)
plotSimilarFit(targ,prof,cexlabs=format[[f]]$cexlabs,mainadd="TELO")
dev.off()

cairo_pdf(file=format[[f]]$pfname,height=format[[f]]$height,width=format[[f]]$width,pointsize=6)

pch=1:9
pch[pch==3]=0
type="p"
cex=4.0
targ=list(
CHK=list(genes=c("DDC1","RAD24","RAD17"),cols="darkgreen",pchs=pch,ltys=lty,lwds=lwd,cexs=cex,types=type),
NMD=list(genes=c("NAM7","NMD2","UPF3"),cols="dodgerblue",pchs=pch,ltys=lty,lwds=lwd,cexs=cex,types=type),
KU=list(genes=c("YKU70","YKU80"),cols="black",pchs=pch,ltys=lty,lwds=lwd,cexs=cex,types=type),
MRX=list(genes=c("MRE11","XRS2","RAD50"),cols="red",pchs=pch,ltys=lty,lwds=lwd,cexs=cex,types=type),
TELO=list(genes=c("RAD9","CHK1","TEL1","EXO1"),cols=c("purple","orange","darkcyan","magenta"),pchs=pch,ltys=lty,lwds=lwd,cexs=cex,types=type)
)
plotSimilarFit(targ,prof,cexlabs=format[[f]]$pcexlabs,showgenes=FALSE,legendval="topright")
dev.off()

}
