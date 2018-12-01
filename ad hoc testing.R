##### ad-hoc testing of package functions:


## load permuted results from a real study:
temp = read.csv("example data.csv", stringsAsFactors = F, row.names = 1)
ests = temp[,"log2.fold.change.from.baseline"]
pvals = temp[,"p.value"]
fdrs = p.adjust(pvals,"BH")
names = rownames(temp)
volcplot(ests, pvals, names = names,fdrs = fdrs,
            show.names.top.N = 100)#,   show.names.fdr.thresh = 0.05)


## load gene sets:
gsets = parsegenesetfile("MSigDB_datasets/c2.all.v6.2.symbols.gmt")
gs2 = parsegenesetfile("MSigDB_datasets/c2.cp.biocarta.v6.2.symbols.gmt.txt")




## run gene set plot:
genesetplot(ests = ests, pvals = pvals, fdrs = NULL, names = names, genesets = gsets,
            n.genesets = 10, min.geneset.size = 5, min.geneset.coverage = 0.1, 
            mandatory.genesets = NULL, geneset.ranking.method = "most.significant",
            fdr.lines = c(0.05, 0.5),
            color.genes.up = "firebrick", color.genes.dn = "darkblue", 
            color.background.up = rgb(1,0,0,0.2), color.background.dn = rgb(0,0,1,0.2),
            xlab = "Estimate", ylim = NULL, 
            cex.points = 0.5, cex.genenames = 0.7, cex.legend = 0.5)
  