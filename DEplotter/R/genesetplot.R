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
#' @param fdr.lines A vector of FDR values at which to draw lines
#' @param show.names.top.N An integer: show the gene names for the top N p-values
#' @param show.names.pval.thresh A value in (0,1): show names for all genes below this p-value
#' @param show.names.fdr.thresh A value in (0,1): show names for all genes below this FDR
#' @param color.top.up The color with which to show the names of the top up-regulated genes
#' @param color.top.dn The color with which to show the names of the top down-regulated genes
#' @param xlab The horizontal axis label
#' @param cex.points The size of the points in the plot
#' @param cex.genenames The size of the genenames in the plot
#' @param cex.legend The size of the legend detailing the FDR cutoffs
#' 
#' @param ... Arguments passed to plot()
#' @return A volcano plot. 
#' @export
genesetplot = function(ests, pvals, names, fdr.lines = c(0.05, 0.5),
                    show.names.top.N = 100, show.names.pval.thresh = NULL, show.names.fdr.thresh = NULL,
                    color.top.up = "firebrick", color.top.dn = "darkblue", xlab = "Estimate",
                    cex.points = 0.5, cex.genenames = 0.5, cex.legend = 0.5, ...){
  
  
  
}