library(shiny)
library(xtable)
source("trimComplexes.R",local=TRUE)
source("profilingFunctions.R",local=TRUE)

shinyServer(function(input, output) {

	output$ui <- renderUI({
    if (is.null(input$checkGroup))
      return()

    # Depending on length of input$checkGroup, we'll generate a different
    # UI component and send it to the client.
	lengrp=length(input$checkGroup)
	if (lengrp==0) return()
	if (lengrp>2) lengrp=3
	
	lenplt=vals$nExp
	print(lengrp)
	print(lenplt)
	default="mahalanobis"
	if (lenplt<2) default="distance"  
	
	switch(lengrp,{
			radioButtons("dType",
			label=h4("Distance Measure"),
            choices=c("Euclidean" = "distance"),
			selected="distance"
			   )		
		},{
			radioButtons("dType",
			label=h4("Distance Measure"),
            choices=c("Euclidean" = "distance","Mahalanobis" = "mahalanobis"),
			selected=default
			   )		
		},{
			radioButtons("dType",
			label=h4("Distance Measure"),
            choices=c("Euclidean" = "distance","Correlation" = "correlation","Mahalanobis" = "mahalanobis"),
			selected=default
			   )		
	})
  })
		
		updateScreenSelection <- reactive({
			tmp = input$drawplot
			return(isolate({scrNms[as.numeric(input$checkGroup)]}))
		})
		
		updateScreens <- reactive({
			vals = list()
			screens=updateScreenSelection()
			# Read in all data once
			dfs=list()
			#if(input$fDef=="nAUC") flist=fnames_nAUC
			#if(input$fDef=="MDRMDP") flist=fnames_MDRMDP
			
			for(f in flist){
				dfnew=parseFitReport(f)
				dfs[[dfnew$ScrLabel[1]]]=dfnew
			}

			# Merge all datasets by ORF
			fits=Reduce(function(...) {merge(...,by=c("ORF","Gene"),all=TRUE)}, dfs)
			colroots=colnames(dfs[[1]])[3:length(colnames(dfs[[1]]))] 
			cnames=c("ORF","Gene",c(t(sapply(colroots,paste,names(dfs),sep="."))))
			colnames(fits)=cnames

			# Eliminate some duplicated genes
			fits=fits[!duplicated(fits$Gene),]
			vals$prof=makeProfiles(fits,screens,input$pType)
			vals$nExp=dim(vals$prof)[2]-2
			if (is.null(input$dType)){
				dT="mahalanobis"
			}else{
				dT=input$dType
			}
			vals$sims=similarities(vals$prof,dT)
			colnames(vals$prof) = beautifulNames[colnames(vals$prof)]
			vals$simtarg=c()
			return(vals)
		})
		
		updateTarget <- reactive({
			vals <<- updateScreens()
			if(input$ggroup==1){
				gene=parseGeneList(input$glist,gdict)
			}else{
				gene=parseGeneList(fc$GroupORFs[as.numeric(input$ggroup)],gdict)
			}
			if((length(gene)==1)&(gene%in%rownames(vals$sims))) {
				simtarg = findSimilar(gene,vals$sims)
			}else{
				simtarg=c()
			}
			return(simtarg)
		})

		output$downloadPlot <- downloadHandler(
			filename = function() { "FitnessProfile.pdf" },
			content = function(file) {
			vals <<- updateScreens()
			if(input$ggroup==1){
				gene=parseGeneList(input$glist,gdict)
			}else{
				gene=parseGeneList(fc$GroupORFs[as.numeric(input$ggroup)],gdict)
			}
			simtarg = updateTarget()
			near = head(simtarg,n=input$nearNum)
			far = tail(simtarg,n=input$farNum)
			cairo_pdf(file,width=14,height=10,onefile=TRUE)
				print(plotSimilarFit(gene,vals$prof,nearest=near,farthest=far))
			dev.off()
		})
		
		output$profiles <- renderPlot({
			vals <<- updateScreens()
			if(input$ggroup==1){
				gene=parseGeneList(input$glist,gdict)
			}else{
				gene=parseGeneList(fc$GroupORFs[as.numeric(input$ggroup)],gdict)
			}
			simtarg = updateTarget()
			near = head(simtarg,n=input$nearNum)
			far = tail(simtarg,n=input$farNum)
			plotSimilarFit(gene,vals$prof,nearest=near,farthest=far) 
		})
		
		output$ranks <- renderPlot({
			vals <<- updateScreens()
			if(input$ggroup==1){
				gene=parseGeneList(input$glist,gdict)
			}else{
				gene=parseGeneList(fc$GroupORFs[as.numeric(input$ggroup)],gdict)
			}
			gene=intersect(gene,rownames(vals$sims))
			plotGenomewideSimilarity(gene,vals$sims)
		})
			
		output$nearest <- renderText({
			vals <<- updateScreens()
			simtarg = updateTarget()
			if((length(simtarg)>0)&(input$nearNum>0)){
				near = head(simtarg,n=input$nearNum)
				dftab=data.frame('Gene'=paste('<a href="http://www.yeastgenome.org/locus/',tolower(names(near)),'/overview" target="_blank">',tolower(names(near)),'</a>',sep=""),Distance=near,row.names=NULL)
				colnames(dftab)=c("Gene Name","Distance")
				return(as.character(print(xtable(dftab,align=c("r","c","c"),digits=4),type="html",sanitize.text.function = function(x){x})))
			}else{return("")}
		})
			
		output$farthest <- renderText({
			vals <<- updateScreens()
			simtarg = updateTarget()
			if((length(simtarg)>0)&(input$farNum>0)){
				far = rev(tail(simtarg,n=input$farNum))
				dftab=data.frame(Gene=paste('<a href="http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=',tolower(names(far)),'" target="_blank">',tolower(names(far)),'</a>',sep=""),Distance=far,row.names=NULL)
				colnames(dftab)=c("Gene Name","Distance")
				return(as.character(print(xtable(dftab,align=c("r","c","c"),digits=4),type="html",sanitize.text.function = function(x){x})))
			}else{return("")}
		})
})
