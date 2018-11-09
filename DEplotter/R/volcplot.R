#' Draw a single volcano plot given vectors of p-values and estimates
#'
#' This function draws a basic volcano plot given differential expression results. Options are given
#' to control the threshold at which gene names are shown, as well as basic plotting options. 
#'
#' @param ests A vector of parameter estimates for each gene
#' @param pvals A vector of the p-values from each gene
#' @param fdrs A vector of False Discovery Rate values. If not give, p.adjust is used with method = "BH"
#' @param names A vector 
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
volcplot = function(ests, pvals, fdrs, names, fdr.lines = c(0.05, 0.5),
                    show.names.top.N = 100, show.names.pval.thresh = NULL, show.names.fdr.thresh = NULL,
                    color.top.up = "firebrick", color.top.dn = "darkblue", xlab = "Estimate",
                    cex.points = 0.5, cex.genenames = 0.5, cex.legend = 0.5, ...){
  # QC the input:
  if(length(ests) != length(pvals)){stop("ests and pvals don't have the same length.")}
  if(length(ests) != length(names)){stop("ests and names don't have the same length.")}
  # make sure every gene has a name:
  no.valid.name = union(which(is.na(names)), which(gsub(" ","",names) == ""))
  names[no.valid.name] = paste0("gene", no.valid.name)
  
  # calc FDRs if not provided:
  if (length(fdrs) == 0){
    fdrs = p.adjust(pvals, "BH")
  }
  # get p-values corresponding to those FDR cutoffs:
  fdr.line.positions = c()
  if (length(fdr.lines) > 0){
    fdr.lines = fdr.lines[order(fdr.lines)]
    for (i in 1:length(fdr.lines)){
      fdr.line.positions[i] = suppressWarnings(-log10(max(pvals[fdrs < fdr.lines[i]])))
    }
    
  }
  
  # determine which gene names to show using the most permissive criterion given by the user. 
  show.text = rep(FALSE, length(ests))
  # if the show.names.top.N arg is valid, apply it:
  if (length(show.names.top.N) == 1){
    show.text[order(pvals)[1:show.names.top.N]] = TRUE
  }
  # if the show.names.pval.thresh arg is valid, apply it:
  if (length(show.names.pval.thresh) == 1){
    show.text[which(pvals < show.names.pval.thresh)] = TRUE
  }
  # if the show.names.fdr.thresh arg is valid, apply it:
  if (length(show.names.fdr.thresh) == 1){
    show.text[which(fdrs < show.names.fdr.thresh)] = TRUE
  }
  
  # draw the volcano plot
  plot(ests, -log10(pvals), pch=16, cex = cex.points, col=c(rgb(0,0,0,.5),"white")[1 + show.text],
       xlab = xlab,
       ylab = "-log10(p-value)", ...)
  
  # add the FDR cutoff lines
  for (i in 1:length(fdr.line.positions)){
    abline(h = fdr.line.positions[i], lty = 1+i, col = "grey20")
  }
  
  # add the gene names:
  if (sum(show.text) > 0){
    text(ests[show.text], -log10(pvals[show.text]), names[show.text],
       cex=cex.genenames, col=c(color.top.dn, color.top.up)[1+(ests[show.text]>0)])  
  }
  
  # draw a legend for the FDR cutoffs:
  if (length(fdr.lines) > 0){
    legend("bottomleft", lty = 1 + (1:length(fdr.lines)), legend = paste0("FDR = ", fdr.lines))
  }
}
  