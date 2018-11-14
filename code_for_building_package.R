##### this script contains the code used to build the package.


### strategy for shiny-in-a-package deployment: taken from https://www.r-bloggers.com/packaging-shiny-applications-a-deep-dive/
### code written to emulate: https://github.com/MangoTheCat/shinyAppDemo


#### package building code: from https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/



library("devtools")
library(roxygen2)


create_package("DEplotter")

setwd("DEplotter")
document()
usethis::use_testthat()
devtools::test()

setwd("..")
install("DEplotter")
library("DEplotter")




# to run the app:
ReadLossPlotter::launchApp()

#### idea: in gene set plot, show a background bar for the mean up-pval to the mean dn-pval, wiht the bar colored by whichever directoin is more impressive?

#** ---> try it out!

#**- geneset plotter: format genesets object the same way as msigdb downloads
#**- geneset plotter: provide args for controlling which gene sets go in there


##### code for exploring EGSEAdata datasets
# note: pathways available in library("EGSEAdata"); DOI:10.18129/B9.bioc.EGSEAdata 
library(EGSEAdata)

??`EGSEAdata-package`
??EGSEAdata
egsea.data()
?gsetdb.human
?kegg.pathways
?msigdb

objs = c("gsetdb.human", "kegg.pathways", "msigdb")
for(o in objs){
  print((get(o)[[1]][[1]]))
}

## structure we atually want: just a longform df: columns are: species, database name, gene set name, gene name
# (or each geneset could be a list, and a separate lookup table could have species and database-relevant info (and maybe even num. genes))
# basic user goals:
# - try out different thematic genesets
# - 

### simpler alternative: just provide KEGG with an easy format; let user add own if they want
