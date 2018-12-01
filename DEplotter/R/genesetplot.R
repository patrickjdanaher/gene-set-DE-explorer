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
#' @param draw.volcano Logical, whether to draw a volcano plot alongside it. If this is selected,
#'        then various formatting choices will be made with layout() and par()$mar.
#' @param volcano.to.geneset.width.ratio A number in (0,1), how much of the total width to devote
#'       to the volcano plot. (Default 0.2.)
#' @param fdr.lines A vector of FDR values at which to draw lines
#' @param fdr.legend Logical, whether to show a legend for the FDR lines
#' @param color.genes.up The color with which to show the names of the up-regulated genes
#' @param color.genes.dn The color with which to show the names of the down-regulated genes
#' @param color.background.up The color of background bars showing mean -log10 p-values of up-regulated genesets
#' @param color.background.dn The color of background bars showing mean -log10 p-values of down-regulated genesets
#' @param cex.genenames The size of the genenames in the plot
#' @param cex.genesetnames The size of the genesetnames in the horizontal axis labels
#' @param cex.legend The size of the legend showing FDR values
#' @param bottom.margin If not NULL, par()$mar will be reset to this.
#' @param show.names.top.N Only used in the volcano plot. An integer: show the gene names for the
#'        top N p-values.
#' @param show.names.pval.thresh Only used in the volcano plot. A value in (0,1): show names
#'        for all genes below this p-value.
#' @param show.names.fdr.thresh Only used in the volcano plot. A value in (0,1): show names for all
#'        genes below this FDR.
#' @param xlab Only used in the volcano plot. The horizontal axis label for the volcano plot.
#' @return A plot of gene sets.
#' @export
genesetplot = function(ests, pvals, fdrs = NULL, names, genesets,
                       n.genesets = 10, min.geneset.size = 5, min.geneset.coverage = 0.3,
                       mandatory.genesets = NULL,
                       geneset.ranking.method = "most.significant",
                       draw.volcano = TRUE, volcano.to.geneset.width.ratio = 0.2,
                       fdr.lines = c(0.05, 0.5), fdr.legend = TRUE,
                       color.genes.up = "firebrick", color.genes.dn = "darkblue",
                       color.background.up = rgb(1,0,0,0.2), color.background.dn = rgb(0,0,1,0.2),
                       ylim = NULL, cex.points = 0.5, cex.genenames = 0.7, cex.genesetnames = 0.6,
                       cex.legend = 0.5, bottom.margin = 12,
                       show.names.top.N = 100, show.names.pval.thresh = NULL,
                       show.names.fdr.thresh = NULL,xlab = "Estimate"){
  ## name ests and pvals with gene names:
  if (length(names) > 0){
    names(ests) = names(pvals) = names
  }
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
    if (length(genesets[[gsname]]) < min.geneset.size){
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
  show = show[1:max(length(mandatory.genesets), n.genesets)]

  # parse the ylim: take what's been given, or if NULL, get the range of the available ones
  if (length(ylim) == 0){
    # get all relevant pvals to find their range:
    temp = c()
    for (name in show){
      temp = c(temp, pvals[match(genesets[[name]], names)])
      ylim = c(0, -log10(min(temp, na.rm = T)))
    }
  }
  # for each geneset, calculate the mean -log10(pval) of the negative and positive genes
  means = means.neg = means.pos = c()
  for (name in show){
    tempgenes = genesets[[name]]
    tempgenes.neg = tempgenes[ests[tempgenes] < 0]
    tempgenes.pos = tempgenes[ests[tempgenes] > 0]
    means[name] = mean(-log10(pvals[tempgenes]))
    means.neg[name] = mean(-log10(pvals[tempgenes.neg]))
    means.pos[name] = mean(-log10(pvals[tempgenes.pos]))
  }

  ## draw the plot:

  # save the original margins in order to revert to them:
  mars0 = par()$mar

  # if draw.volcano has been selected:
  if (draw.volcano){
    # lay out the plotting window:
    layout(matrix(c(1,2),1), widths = c(volcano.to.geneset.width.ratio, 1 - volcano.to.geneset.width.ratio))

    # define margins for the volcplot:
    mars = mars0
    if (length(bottom.margin)>0){
      mars[1] = bottom.margin
    }
    mars[4] = 0.5
    par(mar = mars)
    # draw the volcano plot:
    volcplot(ests = ests, pvals = pvals, fdrs = fdrs, names = names, fdr.lines = fdr.lines,
            show.names.top.N = show.names.top.N, show.names.pval.thresh = show.names.pval.thresh,
            show.names.fdr.thresh = show.names.fdr.thresh, xlab = xlab,
            color.top.up = color.genes.up, color.top.dn = color.genes.dn,
            cex.points = 0.5, cex.genenames = cex.genenames, cex.legend = 0.5, ylim = ylim)

    # define margins for the geneset plot:
    mars[2] = 0.5
    par(mar = mars)

    # disable the FDR legend for the genesetplot:
    fdr.legend = FALSE
  }

  # draw the geneset plot:
  #bp = barplot(means, xaxt = "n", ylab = "-log10(p-value)", main = "", col = 0,
  #             border = F, ylim = ylim, axes = !draw.volcano
  plot(means, xaxt = "n", ylab = "-log10(p-value)", main = "", col = 0,
         ylim = ylim, yaxt = 'n', xlab = "")
  bp = 1:length(means)
  # add FDR lines:
  if (length(fdr.lines) > 0){
    # get FDRs if not provided:
    if (length(fdrs) == 0){
      fdrs = p.adjust(pvals, "BH")
    }
    # calculate the positions of the FDR lines:
    fdr.line.positions = c()
    fdr.lines = fdr.lines[order(fdr.lines)]
    for (i in 1:length(fdr.lines)){
      fdr.line.positions[i] = suppressWarnings(-log10(max(pvals[fdrs < fdr.lines[i]])))
    }
    # add the FDR cutoff lines
    for (i in 1:length(fdr.line.positions)){
      abline(h = fdr.line.positions[i], lty = 1+i, col = "grey20")
    }
  }

  # draw boxes for mean -log10pvals in up vs. down genes:
  boxwidth = (bp[2] - bp[1]) / 2
  # draw alternating background boxes:
  for (i in 1:length(show)){
    if (i%%2 == 0){
        rect(bp[i] - boxwidth, ylim[1], bp[i] + boxwidth, ylim[2], col = rgb(0,0,0,0.1), border = F)
    }
  }
  #for (i in 1:length(show)){
  #  if(means.pos[i] > means.neg[i]){
  #    rect(bp[i] - boxwidth, means.neg[i], bp[i] + boxwidth, means.pos[i], col = color.background.up, border = F)
  #  }
  #  if(means.pos[i] < means.neg[i]){
  #    rect(bp[i] - boxwidth, means.pos[i], bp[i] + boxwidth, means.neg[i], col = color.background.dn, border = F)
  #  }
  #}
  # draw gene names:
  for (i in 1:length(show)){
    name = show[i]
    tempgenes = genesets[[name]]
    text(rep(bp[i], length(tempgenes)), -log10(pvals[tempgenes]), tempgenes, cex = cex.genenames,
         col = c(color.genes.dn, color.genes.up)[1 + (ests[tempgenes] > 0)])
  }
  # draw a y-axis if the volcano plot hasn't been drawn:
  if (!draw.volcano){
    axis(2)
  }
  # draw the x-axis, the geneset names:
  axis(1, at = bp, show, las = 2, cex.axis = cex.genesetnames)
  # add fdr lines legend:
  if (fdr.legend & (length(fdr.lines)>0)){
    legend("bottomleft", lty = 1 + (1:length(fdr.lines)), legend = paste0("FDR = ", fdr.lines))
  }

  # revert to old margins:
  par(mar = mars0)
}
