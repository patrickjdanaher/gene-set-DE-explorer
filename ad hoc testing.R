##### ad-hoc testing of package functions:


ests = -100 : 100
pvals = runif(201, 0, 1)/pmax(pmin(abs(ests),50),1)

pvals[abs(ests)>20] = pvals[abs(ests)>20]


volcplot(ests, pvals, names = paste0("gene",1:length(ests)),
      #   show.names.top.N = 20,
         show.names.fdr.thresh = 0.1)
