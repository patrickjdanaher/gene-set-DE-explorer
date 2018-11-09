#' Draw a single volcano plot given vectors of p-values and estimates
#'
#' This function draws a basic volcano plot given differential expression results. Options are given
#' to control the threshold at which gene names are shown, as well as basic plotting options. 
#'
#' @param ests A vector of parameter estimates for each gene
#' @param pvals A vector of the p-values from each gene
#' @param names A vector 
#' @param genesets A character matrix of gene set/ gene name key value pairs, 
#'                 with 2 columns, the first giving gene set name and the second giving gene name.
#'                 Defaults to KEGG pathways, which is loaded with the package.
#' @param fdr.lines A vector of FDR values at which to draw lines
#' @param color.genes.up The color with which to show the names of the up-regulated genes
#' @param color.genes.dn The color with which to show the names of the down-regulated genes
#' @param color.background.up 
#' @param color.background.dn 
#' @param cex.genenames The size of the genenames in the plot
#' 
#' @return A plot of gene sets. 
#' @export
genesetplot = function(ests, pvals, names, fdr.lines = c(0.05, 0.5),
                    show.names.top.N = 100, show.names.pval.thresh = NULL, show.names.fdr.thresh = NULL,
                    color.top.up = "firebrick", color.top.dn = "darkblue", xlab = "Estimate",
                    cex.points = 0.5, cex.genenames = 0.5, cex.legend = 0.5, ...){
  
  
  
}