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


