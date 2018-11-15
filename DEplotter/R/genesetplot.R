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
#' @param geneset.ranking.method One of "most.significant", "most.unidirectional".
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
  
  ## choose which genesets to plot:
  
  
  ## draw the plot:
  
}