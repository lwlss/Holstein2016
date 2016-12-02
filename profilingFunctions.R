parseFitReport=function(fname){
	con=file(fname,open="r")
	lskip=1
	linedat=""
	header=list()
	while(!grepl("##########",linedat)) {
		linedat=readLines(con,n=1,warn=FALSE)
		vvec=strsplit(linedat,": ")[[1]]
		if(!grepl("##########",vvec[1])) header[[vvec[1]]]=vvec[2]
		lskip=lskip+1
	}
	close(con)
	names(header)=gsub(" ","_",names(header),fixed=TRUE)
	df=read.delim(fname,skip=lskip-1,stringsAsFactors=FALSE,sep="\t",header=TRUE)
	for (h in names(header)) df[[h]]=header[[h]]
	df$ScrLabel=paste(df$Screen_name,df$Treatment)
	return(df)
}

makeProfiles=function(fits,scrNms,genChar="MedianFit"){
	# First attempt, ignoring uncertainty on means
	df=list()
	for(f in scrNms){
		df[[f]]=fits[[paste(genChar,f,sep=".")]]
	}
	df=data.frame(df)
	rownames(df)=fits$Gene
	# Cannot calculate correlation between rows where one row is filled with identical values 
	# (e.g. all strains have zero fitness)
	df$SD=apply(df,1,sd,na.rm=TRUE)
	df$ORF=fits$ORF
	if(length(scrNms)>1) df=df[df$SD>0,]
	return(df)
}

parseGeneList=function(glist,gdict){
	# Parse comma, space, tab, newline delimited list of genes
	# Replacing systematic gene names with standard where appropriate
	# and stripping any gene names not found in library
	glist=gsub(" ",",",glist)
	glist=gsub("\t",",",glist)
	glist=gsub("\n",",",glist)
	gene=toupper(strsplit(glist,",")[[1]])
	gene[!gene%in%gdict]=gdict[gene[!gene%in%gdict]]
	gene=gene[!is.na(gene)]
	return(gene)
}

fastPwMahal=function(x1_raw,invCovMat) {
	x1=as.matrix(x1_raw)
	SQRT=with(svd(invCovMat), u %*% diag(d^0.5) %*% t(v))
	as.matrix(dist(x1 %*% SQRT))
}

# Generate correlation, distance and Mahalanobis matrices
similarities=function(prof,type="mahalanobis"){
	nExp=dim(prof)[2]-2
	# Allow us to use complete.cases where SD=NA (i.e. one experiment only)
	prof$SD=0
	prof=prof[complete.cases(prof),]
	res=matrix()
	if(type=="correlation") res=1-cor(t(as.matrix(prof[,1:nExp])),use="complete.obs")
	if(type=="distance") res=as.matrix(dist(prof[,1:nExp]))
	if(type=="mahalanobis") {
		invCovMat = solve(cov(as.matrix(prof[,1:nExp])))
		res=fastPwMahal(prof[,1:nExp],invCovMat)
	}
	# For some reason, row names are lost if nExp == 1
	if(nExp==1) {
		rownames(res)=rownames(prof)
		colnames(res)=rownames(prof)
	}
	return(res)
}

# Find genes most similar to a target
findSimilar=function(target,sims){
	res=sims[,target][order(sims[,target],decreasing=FALSE)]
	return(res)
}

# Plot fitness profiles for most similar and most highly correlated genes
plotSimilarFit=function(targ,prof,nearest=c(),farthest=c(),ylim=c()){
	# Number of experiments
	nExp=dim(prof)[2]-2
	
	if(length(ylim)==0) ylim=quantile(as.numeric(unlist(prof[,1:nExp])),c(0,1),na.rm=TRUE)
	#ylim=pmin(ylim,150)
	# Visualise similarity of similar genes
	nearprof=data.frame(prof[rownames(prof)%in%names(nearest),1:nExp])
	farprof=data.frame(prof[rownames(prof)%in%names(farthest),1:nExp])
	# For some reason, if prof only contains one experiment, rownames are obliterated, so replace them:
	rownames(nearprof)=rownames(prof)[rownames(prof)%in%names(nearest)]
	rownames(farprof)=rownames(prof)[rownames(prof)%in%names(farthest)]
	if(length(targ)==1) {
		mlab=paste(tolower(targ),"profile")
		cols=c("red")
	}else{
		mlab="Fitness profiles"
		cols=rainbow(length(targ))
	}
	cexlabs=1.5
	if(length(targ)>1) cexlabs=1.5
	
	# Similarity
	op=par(mai=c(4,2,1,1))
	profdat=prof[,1:nExp]
	delta=0.5
	bcol="gray92"
	pcol=rgb(1,0.65,0.65)
	plot(NULL,xlab="",ylab="Fitness",main=mlab,ylim=ylim,cex.lab=2,xlim=c(1-delta,nExp+delta),cex.main=2,xaxt="n",yaxt="n", mgp=c(5, 2, 0))
	maxvals=apply(profdat,2,max,na.rm=TRUE); minvals=apply(profdat,2,min,na.rm=TRUE)
	polygon(c(1-delta,1:nExp,nExp+delta,nExp+delta,rev(1:nExp),1-delta),c(maxvals[1],maxvals,maxvals[length(maxvals)],minvals[length(minvals)],rev(minvals),minvals[1]),col=bcol,border=bcol)
	
	boxplot(profdat,col="lightblue",notch=FALSE,outline=FALSE,border="black",las=2,range=1.5,cex.axis=2,add=TRUE,lwd=1.5,pars = list(boxwex = 0.45, staplewex = 0.5, outwex = 0.5))
	#points(maxvals,pch="-",col="darkgrey",cex=3)
	#points(minvals,pch="-",col="darkgrey",cex=3)
	
	if(length(nearest)>0) for(i in 1:dim(nearprof)[1]) points(1:nExp,nearprof[i,],type="b",col="pink",lwd=2)
	if(length(farthest)>0) for(i in 1:dim(farprof)[1]) points(1:nExp,farprof[i,],type="b",col="grey",lwd=2)
	for(i in seq_along(targ)) points(1:nExp,prof[targ[i],1:nExp],col=cols[i],type="b",lwd=3)
	if(dim(nearprof)[1]>0){
		text(nExp,nearprof[,nExp],tolower(rownames(nearprof)),cex=cexlabs,pos=4,col="pink")
		text(1,nearprof[,1],tolower(rownames(nearprof)),cex=cexlabs,pos=2,col="pink")
	}
	if(dim(farprof)[1]>0){
		text(nExp,farprof[,nExp],tolower(rownames(farprof)),cex=cexlabs,pos=4,col="grey")
		text(1,farprof[,1],tolower(rownames(farprof)),cex=cexlabs,pos=2,col="grey")
	}
	for(i in seq_along(targ)){
		text(nExp,prof[targ[i],nExp],tolower(targ[i]),cex=cexlabs,pos=4,col=cols[i])
		text(1,prof[targ[i],1],tolower(targ[i]),cex=cexlabs,pos=2,col=cols[i])
	}
	par(op)
}

# Plot similarity rank
plotGenomewideSimilarity=function(targ,sims){
	ylim=range(sims)
	if(length(targ)==1){
		ylab=paste("Difference from",tolower(targ))
		xlab=paste("ORFs ranked by difference from",tolower(targ))
		cols=c("red")
	}else{
		ylab="Difference from target"
		xlab="ORFs ranked by difference from target"
		cols=rainbow(length(targ))
	}
	op=par(mai=c(1,1,1,1))
	plot(NULL,xlab=xlab,ylab=ylab,main="Genome-wide distribution of differences",type="l",ylim=ylim,xlim=c(0,1.1*dim(sims)[2]),cex.axis=2,cex.lab=2,cex.main=2)
	for(i in seq_along(targ)){
		sim=findSimilar(targ[i],sims)
		points(sim,col=cols[i],pch=16,cex=0.2)
		#text(length(sim),sim[length(sim)],tolower(targ[i]),cex=1.5,pos=4,col=cols[i])
	}
	if(length(targ)>0) legend("topleft",tolower(targ),col=cols,lwd=1)
	par(op)
}

