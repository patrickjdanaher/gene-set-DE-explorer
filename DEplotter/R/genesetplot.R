#' Draw a geneset plot given p-values and gene set membership.
#'
#' This function draws a gene set plot given differential expression results. 
#' Options are given to control how gene sets are included. 
#'
#' @param ests A vector of parameter estimates for each gene
#' @param pvals A vector of the p-values from each gene
#' @param names A vector 
#' @param genesets A list of genesets. Each entry's name is a gene set name, and each entry is a character vector of gene names.
#'                 Defaults to KEGG pathways, which is loaded with the package.
#' @param n.genesets How many genesets to show
#' @param min.geneset.size Minimum number of genes for a geneset to have in order to be considered
#' @param min.geneset.coverage A number in (0,1) - excluded gene sets without at least this many of their genes present in the data
#' @param mandatory.genesets Character vector of gene set names to include
#' @param geneset.ranking.method One of "most.significant", "most.unidirectional" (extremely up- or down-regulated genesets),
#'        "most.sig.in.each.direction" (get the most significant up- and down-regulated genesets)
#' @param fdr.lines A vector of FDR values at which to draw lines
#' @param color.genes.up The color with which to show the names of the up-regulated genes
#' @param color.genes.dn The color with which to show the names of the down-regulated genes
#' @param color.background.up The color of background bars showing mean -log10 p-values of up-regulated genesets
#' @param color.background.dn The color of background bars showing mean -log10 p-values of down-regulated genesets
#' @param cex.genenames The size of the genenames in the plot
#' 
#' @return A plot of gene sets. 
#' @export
genesetplot = function(ests, pvals, names, genesets,
                       n.genesets = 10, min.geneset.size = 5, min.geneset.coverage = 0.3, 
                       mandatory.genesets = NULL, geneset.ranking.method = "most.significant",
                       fdr.lines = c(0.05, 0.5),
                       color.genes.up = "firebrick", color.genes.dn = "darkblue", 
                       color.background.up = rgb(1,0,0,0.2), color.background.dn = rgb(0,0,1,0.2),
                       xlab = "Estimate",
                       cex.points = 0.5, cex.genenames = 0.5, cex.legend = 0.5, ...){
  
  ## take intersection of genesets with data:
  for (gsname in setdiff(names(genesets), mandatory.genesets)){
    # exclude genesets with too little intersection with data:
    if (length(intersect(genesets[[gsname]], names)) / length(genesets[[gsname]]) < min.geneset.coverage){
      genesets[[gsname]] = NULL
    }
    # only keep the intersection:
    genesets[[gsname]] = intersect(genesets[[gsname]], names)
  }
  
  ## exclude genesets with insufficient genes
  for (gsname in setdiff(names(genesets), mandatory.genesets)){
    # exclude genesets with too little intersection with data:
    if (length(genesets[[gsname]] < min.geneset.size){
      genesets[[gsname]] = NULL
    }
  }
  
  ## choose which genesets to plot:
  scores = c()
  if (geneset.ranking.method == "most.significant"){
    for (gsname in names(genesets)){
      tempgenes = genesets[[gsname]]
      scores[gsname] = mean(-log10(pvals[tempgenes]))
    }
  }
  if (geneset.ranking.method == "most.unidirectional"){
    for (gsname in names(genesets)){
      tempgenes = genesets[[gsname]]
      uppvals = -log10(replace(pvals[tempgenes], ests[tempgenes] < 0, 1))
      dnpvals = -log10(replace(pvals[tempgenes], ests[tempgenes] > 0, 1))
      scores[gsname] = max(mean(uppvals), mean(dnpvals))
    }
  }
  if (geneset.ranking.method == "most.sig.in.each.direction"){
    upscores = c()
    dnscores = c()
    for (gsname in names(genesets)){
      tempgenes = genesets[[gsname]]
      uppvals = -log10(replace(pvals[tempgenes], ests[tempgenes] < 0, 1))
      dnpvals = -log10(replace(pvals[tempgenes], ests[tempgenes] > 0, 1))
      upscores[gsname] = mean(uppvals)
      dnscores[name] = mean(dnpvals)
    }
    # score each geneset by whichever is more impressive: its rank in the up-regulated scores or among the dn-regulated scores
    scores = apply(cbind(rank(upscores), rank(dnscores)), 1, max)
  }
  
  ## which genesets to plot: any mandatory, plus any additional in decreasing order of score up to n.genesets
  show = c(mandatory.genesets, names(genesets)[order(scores, decreasing = T)])
  show = show[max(length(mandatory.genesets), n.genesets)]
  
  ## draw the plot:
  
  
  
}