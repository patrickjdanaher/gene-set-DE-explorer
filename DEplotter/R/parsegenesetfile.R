#' Parse a gene set file into a gene set object.
#'
#' This function takes a gene set file like those output by msigdb and modifies it into a 2 column mapping of 
#' pathways to their component genes
#'
#' @param filename The filename of a gene set file ala those output by MSigDB
#' @return A 2-column data frame of gene sets and their member genes 
#' @export
parsegenesetfile = function(filename){
  ## load in file
  # get number of columns
  x = scan(filename, what = "", sep = "\n")
  x = strsplit(x, "[ \t]+") # split string by white space
  max.col = max(sapply(x, length))
  # read the file as a table:
  temp = read.table(filename, fill = TRUE, stringsAsFactors = F, header = F, col.names = paste0("x", 1:max.col))
  # remove the second column if it's web addresses, as expected from msigdb:
  if (grepl("http", temp[1,2])){
    temp = temp[, -2, drop = F]
  }
  temp = t(temp)
  
  ## assemble the object for use in plotting, a list of genesets:
  genesets = list()
  for (i in 1:ncol(temp)){
    gsname = temp[1, i]
    genesets[[gsname]] = unique(setdiff(as.vector(temp[-1, i]), ""))
  }
  
  return(genesets)
}
                    

#infile = "c2.all.v6.2.symbols.gmt"



intersectgenesetswithdata = function(genesets, genenames){
  
}