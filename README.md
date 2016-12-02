# Data & source code Holstein et al. 2016

## Profilyzer

This directory contains a [Shiny](https://shiny.rstudio.com/) app for interactive exploration of the [QFA](http://research.ncl.ac.uk) fitness profiles found in the GeneStripped directory below.

![A static version of a profilyzer plot](Demo.png?raw=true)

Above is a static version of a plot that can be produced by profilyzer.

To see profilyzer in action, a live instance of this particular dataset can be found at [this page](http://research.ncl.ac.uk/Holstein2016).  Alternatively, you can download the code & data from this repository and run a local session.

To run profilyzer in a local session: 
* Check out this repository
* Open R session and set the current working directory to the one just above the profilyzer directory
* Ensure that you have shiny installed.  If you do not, then execute the following in the R terminal: `install.packages("shiny")`
* Load the shiny library by executing the following in the R terminal: `library("shiny")`
* Finally, use profilyzer to browse the included LydallLab dataset: `runApp("profilyzer")`

## Data

This repository contains the following data directories:

* AllStripped: QFA fitness estimate reports for each of the independent strains examined in screens carried out for this manuscript.  Fitnesses of strains with deleted genes found within 20kb of screen background mutation and genes on a list of "slow-growing" strains that do not grow well during SGA are stripped.
* CommonStripped: QFA fitness estimate reports for each of the independent strains examined in screens carried out for this manuscript.  Fitnesses of strains with deleted genes on a list of "slow-growing" strains that do not grow well during SGA are stripped.
* GeneStripped: QFA fitness estimate reports for each of the independent strains examined in screens carried out for this manuscript.  Fitnesses of strains with deleted genes found within 20kb of screen background mutation and genes on a list of "slow-growing" strains that do not grow well during SGA and all impossible "double deletions" (e.g. rad9D rad9D)  are stripped.
* GIS: Genetic interaction strength estimates generated by comparing fitness in matched pairs of QFA screens, as used in our online interactive plotting tool [DIXY](http://bsu-srv.ncl.ac.uk/dixy-telo/viz/), for example.  

For detailed description of how to interpret columns and how fitnesses and genetic interaction strengths are calculated, please see the documentation for the [R package](http://qfa.r-forge.r-project.org/) which was used to generate them.

