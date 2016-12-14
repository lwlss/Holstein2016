library(shiny)

shinyUI(fluidPage(title="Profilyzer",
	titlePanel(h1("Fitness Profiling: Responses to defective telomeres in budding yeast")),
	HTML('<p>To explore evidence for genetic interaction in detail, use our fitness plot visualisation tool: <a href="http://bsu-srv.ncl.ac.uk/dixy-telo">DIXY</a></p>'),
	p("Press 'Draw Plot' button after changing input selection in grey box.  Plots and tables are updated in response to all other input in real time."),
	
	sidebarLayout(
		sidebarPanel(
		  checkboxGroupInput("checkGroup", 
			label = h4("Screens to profile"), 
			choices = chooses,
			selected = starting),
			HTML('<br><p>QFA numbers refer to entries in the internal <a href="http://minch-moor.ncl.ac.uk/fmi/iwp/cgi?-db=QFAScreenDatabase&-startsession">Lydall lab QFA database</a>, which contains more detailed information about screens.</p>'),
			actionButton("drawplot", label = "Draw Plot"),
			width=2
		),
		
  mainPanel(
  fluidRow(
		
	column(2,
		radioButtons("pType",
			label=h4("Fitness summary"),
            choices=c("Mean" = "MeanFit","Median" = "MedianFit"),
			selected="MeanFit"
			   ),
		h6("Method of summarising replicate fitness observations for individual genotypes.")#,
#		radioButtons("fDef", 
#			label = h4("Fitness definition"), 
#			choices = c("nAUC"="nAUC","MDRMDP"="MDRMDP"),
#			selected = "nAUC")
		),
	column(2,
		uiOutput("ui"),
		HTML("<h6><a href='https://en.wikipedia.org/wiki/Euclidean_distance'>Euclidean distance</a> discovers profiles that are close or overlap.  <a href='https://en.wikipedia.org/wiki/Correlation_and_dependence'>Correlation</a> discovers profiles with the same pattern, but potentially offset from each other.  <a href='https://en.wikipedia.org/wiki/Mahalanobis_distance'>Mahalanobis distance</a> can be considered as a combination of the two.</h6>")
			   ),
		
	column(3,
		textInput("glist", label = h4("Target gene(s)"),value = "rad24"),
		h6("Nearest & farthest genes are only highlighted when just a single gene is specified in the box above.  Note that nearest and farthest genes can only be identified when target gene is present in all selected screens."),
		selectInput("ggroup", label = h4("Gene group"),choices = choiceList,selected = 1),
		h6("Optionally select groups of related genes from the drop-down list above instead of specifying them manually.  Select 'None' to return to highlighting 'Target gene(s)'")
		), 
    
    column(2,
		numericInput("nearNum",label = h4("No. genes nearest to target"),value = 12,min=0),
		numericInput("farNum",label = h4("No. genes farthest from target"),value = 0,min=0),
		h6("Select number of nearest (or farthest) profiles to plot & tabulate alongside a single 'Target gene'.  Nearest genes may have a similar function to the target.  Farthest genes may have an opposite function, but such observations are probably more difficult to interpret.")
		) 
    ),
	  
   fluidRow(
    plotOutput("profiles", height="1000px", width = "100%"),
	h6("Boxplot comparing average fitness distributions across multiple QFA screens, overlaid with average fitness profiles for individual library deletion strains.  Horizontal black bar is median fitness, edges of blue box correspond to the 1st & 3rd quartiles of each average fitness distribution.  Whiskers extend to most extreme average fitness which is no more than 1.5 times the length of the box away from the box.  Grey background corresponds to the range of average fitnesses in each QFA screen.")
	),
	downloadButton('downloadPlot', 'Download Plot'),
	fluidRow(
	column(2,
	h4("Nearest"),
	htmlOutput("nearest")
	),
	column(2,
	h4("Farthest"),
	htmlOutput("farthest")
	),
	column(8,
	plotOutput("ranks", height="500px", width = "100%"),
	h6("Difference distribution: plot of the magnitude of profile difference from target against difference rank.  Examining this plot helps us to decide whether library deletions identified as least different to target are significant outliers or simply part of the average behaviour of deletions across all screens.  The flat part of this curve represents the average profile difference from the target deletion. Deletions with differences from the target around this value have essentially indistinguishable profiles.  A flat curve immediately after the origin suggests that the target profile is not very different from that of many other deletions (e.g. when target is set to his3): ranking by difference will not provide much information and the ranked order of nearest deletions is not likely to be informative.  On the other hand, if the curve has a steep slope near the origin, then deletions on the steep section have profiles that are unusually similar to the target and should be candidates for further investigation.")
	)	
   ),
   HTML('<p>Data, source code & documentation for this instance of profilyzer are hosted on <a href="https://github.com/lwlss/profilyzer">GitHub</a></p>')
   )
)))
