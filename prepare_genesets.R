##### get a few genesets by hand from the EGSEAdata package:


## get KEGG human:
names(kegg.pathways)

names(kegg.pathways[["human"]][["kg.sets"]])
str(kegg.pathways[["human"]][["kg.sets"]])
str(kegg.pathways[["human"]][["kg.sets"]][[1]])
