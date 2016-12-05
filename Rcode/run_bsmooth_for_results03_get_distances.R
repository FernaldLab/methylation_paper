








########################

# compute distances between dmrs and genes
dmrsDistToGenes = lapply(as.list(dmrsGR), function(f) .distance.named(geneGR, f));
# dmrsDistToGenes_summary = as.data.frame(cbind(do.call('rbind', 
#                                                       lapply(dmrsDistToGenes, summary)), 
#                                               num=sapply(dmrsDistToGenes, length)));

dmrsDistToGenes_thresh = list();
test_dists = c(0, 1000, 2000, 3000, 4000, 
               seq(5000,1e5,5000),
               seq(1e5,1e6,1e5)
              # seq(5000, 1e6, 5000) 
              # ,seq(1e6, 7036168, 100000)
               );
#test_dists = seq(0,1e6,1000)
test_dists = unique(test_dists);
for (d in test_dists) {
  cat(as.character(d),'...')
  dmrsDistToGenes_thresh[[as.character(d)]] = lapply(dmrsDistToGenes, function(f) f[f <= d]);
}; rm(d);

# make tables where dmrs are rows and distances are columns

dmrsDistToGenes_thresh_countDf = as.data.frame(sapply(dmrsDistToGenes_thresh[[1]], length));
for (d in 2:length(dmrsDistToGenes_thresh)) {
  dmrsDistToGenes_thresh_countDf = as.data.frame(cbind(dmrsDistToGenes_thresh_countDf,
                                                       sapply(dmrsDistToGenes_thresh[[d]], length) ))
}; rm(d);
names(dmrsDistToGenes_thresh_countDf) = names(dmrsDistToGenes_thresh);

dmrs2=as.data.frame(cbind(dmrs[,7:15], fc=log2(dmrs$group2.mean/dmrs$group1.mean)));
dmrs2 = as.data.frame(cbind(dmrs2[,1:4], 
                            areaStat.abs=abs(dmrs2$areaStat), 
                            meanDiff=dmrs2[,6], 
                            meanDiff.abs=abs(dmrs2$meanDiff), 
                            log2fc=dmrs2[,10], 
                            log2fc.abs=abs(dmrs2$fc)));
dmrsDistToGenes_thresh_countDf_testStats = as.data.frame(apply(dmrs2, 2, 
                                                              # function(ff) apply(dmrsDistToGenes_thresh_countDf[,1:27], 2, 
                                                              function(ff) apply(dmrsDistToGenes_thresh_countDf, 2, 
                                                                                function(f) kruskal.test(ff, 
                                                                                                           f >= mean(f))$p.value  )  ))



tmphigh = apply(apply(dmrsDistToGenes_thresh_countDf, 2, 
                      function(f) f >= mean(f)), 2, 
                function(ff) table(dmrs$direction[ff]));
tmplow = apply(apply(dmrsDistToGenes_thresh_countDf, 2, 
                      function(f) f < mean(f)), 2, 
                function(ff) table(dmrs$direction[ff]));
dodds = c();
for (i in 1:ncol(tmphigh)) {
  dodds = c(dodds, fisher.test(cbind(tmphigh[,i],tmplow[,i]))$p.value)
}; rm(i);
plot(as.numeric(names(dmrsDistToGenes_thresh_countDf)),
     -log10(dodds))


# get expression and fc for genes around dmrs

par(mfrow=c(3,7))
for (wincol in 5:ncol(dmrsDistToGenes_thresh_countDf)) {
  thismean = mean(dmrsDistToGenes_thresh_countDf[,wincol])
  highcountdmrs = dmrsDistToGenes_thresh_countDf[,wincol] >= thismean
  highcountgenes = lapply(dmrsDistToGenes_thresh[[wincol]][highcountdmrs], 
                          function(f) names(f));
  highcountgenesmeans = sapply(highcountgenes, function(f) mean(resexpr_in$log2FoldChange[rownames(resexpr_in) %in% f]));
  lowcountgenes = lapply(dmrsDistToGenes_thresh[[wincol]][!highcountdmrs], 
                          function(f) names(f));
  lowcountgenesmeans = sapply(lowcountgenes, function(f) mean(resexpr_in$log2FoldChange[rownames(resexpr_in) %in% f]));
  vals = c(highcountgenesmeans, lowcountgenesmeans);
  facts = c(rep('high',length(highcountgenesmeans)), rep('low',length(lowcountgenesmeans)))
  # highcountgenes = unique(unlist(sapply(dmrsDistToGenes_thresh[[wincol]][highcountdmrs], 
  #                                       function(f) names(f))));
  # lowcountgenes = unique(unlist(sapply(dmrsDistToGenes_thresh[[wincol]][!highcountdmrs], 
  #                                      function(f) names(f))));
  # vals = c(resexpr_in$baseMean[rownames(resexpr_in) %in% highcountgenes],
  #          resexpr_in$baseMean[rownames(resexpr_in) %in% lowcountgenes]);
  # vals = c(resexpr_in$log2FoldChange[rownames(resexpr_in) %in% highcountgenes],
  #          resexpr_in$log2FoldChange[rownames(resexpr_in) %in% lowcountgenes]);
  # facts = c(rep('high', sum(rownames(resexpr_in) %in% highcountgenes)),
  #           rep('low', sum(rownames(resexpr_in) %in% lowcountgenes)))
  verboseBoxplot(vals, as.factor(facts), 
                 #ylim=c(0,10000),
                 xlab='',ylab='',
                 main=paste0(names(dmrsDistToGenes_thresh_countDf)[wincol],'\n'))
}


dmrs3 = as.data.frame(cbind(dmrs2, 
                            scaffold_length=slGeneStats$length[match(dmrs$chr, rownames(slGeneStats))],
                            dmrsDistToGenes_thresh_countDf,
                            overlapTE=rownames(dmrs3) %in% dmrsOV$te$queryHits));
for (i in 1:length(dmrsDistToGenes_thresh)) {
  cat(names(dmrsDistToGenes_thresh)[i],'...')
  dmrs3 = as.data.frame(cbind(dmrs3, t(sapply(lapply(dmrsDistToGenes_thresh[[i]], 
                                                     function(f) resexpr_in[match(names(f), 
                                                                                  rownames(resexpr_in)), ]), 
                                              function(ff) c(mean(ff$baseMean,na.rm=T), 
                                                             mean(ff$log2FoldChange,na.rm=T), 
                                                             mean(abs(ff$log2FoldChange),na.rm=T))))));
  names(dmrs3)[(ncol(dmrs3)-2):ncol(dmrs3)] = c(paste0(names(dmrsDistToGenes_thresh)[i],'.baseMean'),
                                                paste0(names(dmrsDistToGenes_thresh)[i],'.log2FoldChange'),
                                                paste0(names(dmrsDistToGenes_thresh)[i],'.abs.log2FoldChange'))
}; rm(i);


# 
# plot(as.numeric(rownames(dmrsDistToGenes_thresh_countDf_testStats)), 
#      -log10(dmrsDistToGenes_thresh_countDf_testStats$meanDiff), 
#      type='l', ylim=c(0,5.5));
# par(new=T); 
# plot(as.numeric(rownames(dmrsDistToGenes_thresh_countDf_testStats)), 
#      -log10(dmrsDistToGenes_thresh_countDf_testStats$fc), 
#      type='l', col='green', ylim=c(0,5.5))
dev.off();
colors=c('black', 'midnightblue','blue','skyblue','darkgreen','green','grey','purple','orange','pink', 'cyan','lightgreen')
for (i in c(1:3,5,7,9)) {
  if(names(dmrsDistToGenes_thresh_countDf_testStats)[i] %in% c('n','width','invdensity','areaStat.abs')) {
    thiscol = 'darkgrey'
  } else if (names(dmrsDistToGenes_thresh_countDf_testStats)[i]=='meanDiff.abs') {
    thiscol = 'purple'
  } else {
    thiscol = 'green'
  }
  #thiscol = colors[i]
  plot(as.numeric(rownames(dmrsDistToGenes_thresh_countDf_testStats)), 
       -log10(dmrsDistToGenes_thresh_countDf_testStats[,i]), 
       type='b', pch=19, col=thiscol, ylim=c(0,5.5)); cat(paste0(names(dmrsDistToGenes_thresh_countDf_testStats)[i], ':', thiscol,'\n'));par(new=T)
}; rm(i); par(new=T);
plot(as.numeric(names(dmrsDistToGenes_thresh_countDf)),
     -log10(dodds),
     type='b', pch=19, col='darkgrey', ylim=c(0,5.5), lty='dashed');
abline(h=-log10(c(.05,.01,.001)), col='red');
abline(v=c(20000,50000,100000), col='darkgrey');
par(new=F)

dmrsDistToGenes_thresh_countMeans = apply(dmrsDistToGenes_thresh_countDf, 2, mean);


dmrsDistToGenes_thresh_distDf = as.data.frame(sapply(dmrsDistToGenes_thresh[[1]], mean));
for (d in 2:length(dmrsDistToGenes_thresh)) {
  dmrsDistToGenes_thresh_distDf = as.data.frame(cbind(dmrsDistToGenes_thresh_distDf,
                                                      sapply(dmrsDistToGenes_thresh[[d]], mean) ))
}; rm(d);
names(dmrsDistToGenes_thresh_distDf) = names(dmrsDistToGenes_thresh);

dmrsDistToGenes_thresh_summaryDf = as.data.frame(cbind(count=apply(dmrsDistToGenes_thresh_countDf, 2,
                                                                   mean, na.rm=T), 
                                                       dist=apply(dmrsDistToGenes_thresh_distDf, 2,
                                                                  mean, na.rm=T)));
dmrs_scaffold_lengths = slGeneStats_in70$length[match(dmrs$chr, rownames(slGeneStats_in70))];



test_scaffold_lengths = c(10000, 50000, 1e5, 2e5, 5e5, 1e6, 2e6, 3e6, 4e6, 5e6);
#test_scaffold_lengths = seq(10000, 2e5, 10000)
df = dmrsDistToGenes_thresh_countDf;

dfpvals = as.data.frame(matrix(ncol=length(test_scaffold_lengths),
                               nrow=ncol(df),
                               dimnames=list(names(df), as.character(test_scaffold_lengths))));

for (thiswin in 1:length(as.numeric(names(df)))) {
  #cat(thiswin,'...')
  for (thislen in 1:length(test_scaffold_lengths)) {
    #cat(thislen,'...')
    rows = dmrs_scaffold_lengths > test_scaffold_lengths[thislen];
    valvec = c(df[,thiswin], 
               df[rows, thiswin]);
    factvec = c(rep('all',709), 
                rep(test_scaffold_lengths[thislen], sum(rows))); 
    dfpvals[thiswin, thislen] = kruskal.test(valvec, as.factor(factvec))$p.value
  }
}
rm(thiswin, thislen, rows, valvec, factvec)

tmp = dmrsDistToGenes_thresh_summaryDf[seq(1,nrow(dmrsDistToGenes_thresh_summaryDf),5),  ]
WGCNA::verboseScatterplot(tmp$count, 
                          tmp$dist, 
                          frame.plot=F, abline=F, ylab='', pch=21 ,
                          bg=WGCNA::numbers2colors(tmp$dist / as.numeric(rownames(tmp))                      ))



############

#tmp = dmrsDistToGenes_thresh_summaryDf[seq(1,nrow(dmrsDistToGenes_thresh_summaryDf),5),  ]
tmp = dmrsDistToGenes_thresh_summaryDf[1:201,  ]
tmpdf = as.data.frame(matrix(nrow=nrow(tmp),ncol=3,dimnames=list(rownames(tmp),c('actual','shuff','pval'))))
for (i in 1:nrow(tmp)) {
  cat(rownames(tmp)[i],'...')
  thisactual = dmrsDistToGenes_thresh_summaryDf$count[rownames(dmrsDistToGenes_thresh_summaryDf)==rownames(tmp)[i]];
  thisshuff = sapply(shuffDistToGenesAcrossRuns, 
                     function(ff) mean(sapply(ff, 
                                              function(f) sum(f <= as.numeric(rownames(tmp)[i])))));
  tmpdf[i,] = c(thisactual, mean(thisshuff), sum(thisshuff >= thisactual)/length(shuffDistToGenesAcrossRuns))
}
rm(i,thisactual,thisshuff)

shuffDistToGenesAcrossRuns50kbcounts = t(sapply(shuffDistToGenesAcrossRuns, 
                                                function(ff) summary(sapply(ff, 
                                                                            function(f) sum(f <= 50000)))));  
shuffDistToGenesAcrossRuns50kb = lapply(shuffDistToGenesAcrossRuns, function(f) lapply(f, function(ff) ff[ff <= 50000]));
shuffDistToGenesAcrossRuns20kb = lapply(shuffDistToGenesAcrossRuns, function(f) lapply(f, function(ff) ff[ff <= 20000]));
shuffDistToGenesAcrossRuns5kb = lapply(shuffDistToGenesAcrossRuns, function(f) lapply(f, function(ff) ff[ff <= 5000]))

shuffDistToGenesAcrossRuns20kbcounts = t(sapply(shuffDistToGenesAcrossRuns, 
                                                function(ff) summary(sapply(ff, 
                                                                            function(f) sum(f <= 20000))))); 

shuffDistToGenesAcrossRuns5kbcounts = t(sapply(shuffDistToGenesAcrossRuns, 
                                               function(ff) summary(sapply(ff, 
                                                                           function(f) sum(f <= 5000))))); 

shuffDistToGenesAcrossRuns1mbcounts = t(sapply(shuffDistToGenesAcrossRuns, 
                                               function(ff) summary(sapply(ff, 
                                                                           function(f) sum(f <= 1e6))))); 

shuffDistToGenesAcrossRuns50kbdists= t(sapply(shuffDistToGenesAcrossRuns, 
                                                function(ff) summary(sapply(ff, 
                                                                            function(f) summary(f[f <= 50000]))))); 

shuffDistToGenesAcrossRuns5kbdists= t(sapply(shuffDistToGenesAcrossRuns, 
                                              function(ff) summary(sapply(ff, 
                                                                          function(f) summary(f[f <= 5000]))))); 

shuffDistToGenesAcrossRuns1mbdists= t(sapply(shuffDistToGenesAcrossRuns, 
                                             function(ff) summary(sapply(ff, 
                                                                         function(f) summary(f[f <= 1e6]))))); 



# dmrsDistToGenes_thresh_num_dmrs = sapply(dmrsDistToGenes_thresh, 
#                                          function(f) sum(sapply(f, length) > 0));
# 
# dmrsDistToGenes_thresh_avg_genes = sapply(dmrsDistToGenes_thresh, 
#                                           function(f) mean(sapply(f, length)));
# 
# dmrsDistToGenes_thresh_avg_dist = sapply(dmrsDistToGenes_thresh, 
#                                          function(f) mean(sapply(f, mean), na.rm=T));
# 
# dmrsDistToGenes_thresh_num_genes = as.data.frame(sapply(dmrsDistToGenes_thresh, 
#                                                         function(f) sapply(f, length)));
# 
# #
# dmrsDistToGenes_statDf = as.data.frame(cbind(num.dmrs=dmrsDistToGenes_thresh_num_dmrs, 
#                                              avg.dist.to.gene=dmrsDistToGenes_thresh_avg_dist, 
#                                              avg.num.genes=dmrsDistToGenes_thresh_avg_genes));

#
par(mfrow=c(1,3));     
plot(as.numeric(colnames(dmrsDistToGenes_thresh_num_genes)), 
     apply(dmrsDistToGenes_thresh_num_genes[dmrs$direction=='hypo',],
           2,mean), 
     type='l', xlab='distance',ylab='avg.num.genes');         
plot(as.numeric(colnames(dmrsDistToGenes_thresh_num_genes)), 
     apply(dmrsDistToGenes_thresh_num_genes[dmrs$direction=='hyper',],
           2,mean), 
     type='l', xlab='distance',ylab='avg.num.genes');
    
plot(as.numeric(colnames(dmrsDistToGenes_thresh_num_genes)), 
     apply(dmrsDistToGenes_thresh_num_genes[dmrs$chr %in% dmr_scaffolds_1mb_long,],
           2,mean), 
     type='l', xlab='distance',ylab='avg.num.genes');         





#
par(mfrow=c(3,2));
cex=1; cex.axis=1.5; cex.lab=1.5; cex.main=1.5;
maxind = which(dmrsDistToGenes_thresh_num_dmrs == max(dmrsDistToGenes_thresh_num_dmrs))[1]
plot(as.numeric(names(dmrsDistToGenes_thresh_num_dmrs)),
     dmrsDistToGenes_thresh_num_dmrs,
     type='b', frame.plot=F, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main,
     ylab='number of dmrs', xlab='distance with at least 1 gene');
abline(h=max(dmrsDistToGenes_thresh_num_dmrs), 
       v=as.numeric(names(dmrsDistToGenes_thresh_num_dmrs)[maxind]), col='red');
plot(as.numeric(names(dmrsDistToGenes_thresh_num_dmrs))[1:(maxind+1)],
     dmrsDistToGenes_thresh_num_dmrs[1:(maxind+1)],
     type='b', frame.plot=F, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main,
     ylab='number of dmrs', xlab='distance with at least 1 gene');
abline(h=max(dmrsDistToGenes_thresh_num_dmrs), 
       v=as.numeric(names(dmrsDistToGenes_thresh_num_dmrs)[maxind]), col='red');
verboseScatterplot(as.numeric(names(dmrsDistToGenes_thresh_avg_genes)), 
                   dmrsDistToGenes_thresh_avg_genes, 
                   abline=T, abline.col='red', abline.lty='dashed', 
                   type='b', frame.plot=F, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main,
                   xlab='distance', ylab='avg num genes');
abline(v=as.numeric(names(dmrsDistToGenes_thresh_num_dmrs)[maxind]), col='red');
verboseScatterplot(as.numeric(names(dmrsDistToGenes_thresh_avg_genes))[1:(maxind+1)], 
                   dmrsDistToGenes_thresh_avg_genes[1:(maxind+1)], 
                   abline=T, abline.col='red', abline.lty='dashed', 
                   type='b', frame.plot=F, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main,
                   xlab='distance', ylab='avg num genes');
abline(v=as.numeric(names(dmrsDistToGenes_thresh_num_dmrs)[maxind]), col='red');
verboseScatterplot(as.numeric(names(dmrsDistToGenes_thresh_avg_dist)), 
                   dmrsDistToGenes_thresh_avg_dist, 
                   abline=T, abline.col='red', abline.lty='dashed', 
                   type='b', frame.plot=F, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main,
                   xlab='distance', ylab='avg distance to genes');
abline(v=as.numeric(names(dmrsDistToGenes_thresh_num_dmrs)[maxind]), col='red');
verboseScatterplot(as.numeric(names(dmrsDistToGenes_thresh_avg_dist))[1:(maxind+1)], 
                   dmrsDistToGenes_thresh_avg_dist[1:(maxind+1)], 
                   abline=T, abline.col='red', abline.lty='dashed', 
                   type='b', frame.plot=F, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main,
                   xlab='distance', ylab='avg distance to genes');
abline(v=as.numeric(names(dmrsDistToGenes_thresh_num_dmrs)[maxind]), col='red');

dmrsGR_1mb_scaffoldLengths = dmrsGR[seqnames(dmrsGR) %in% rownames(subset(slGeneStatsDmrs, length >= 1e6))]

dmrsDistToGenes_1mb_scaffoldLengths = lapply(as.list(dmrsGR_1mb_scaffoldLengths), 
                                             function(f) .distance.named(geneGR, f));

dmrsDistToGenes_1mb_scaffoldLengths_summary = as.data.frame(cbind(do.call('rbind', 
                                                                          lapply(dmrsDistToGenes_1mb_scaffoldLengths, 
                                                                                 summary)),
                                                                  num=sapply(dmrsDistToGenes_1mb_scaffoldLengths, 
                                                                             length)));


.getDmrGeneDistStatsForCertainLengthScaffolds = function (dmrsGR, geneGR, slGeneStats, reqScaffoldLength=1e6) {
  filteredGR = dmrsGR[seqnames(dmrsGR) %in% rownames(subset(slGeneStats, length >= reqScaffoldLength))];
  geneDistListByDmr = lapply(as.list(filteredGR), 
                             function(f) .distance.named(geneGR, f));
  geneDistListByDmrSummary = as.data.frame(cbind(do.call('rbind', lapply(geneDistListByDmr, summary)),
                                                 num=sapply(geneDistListByDmr, length)));
  return(list(distlist=geneDistListByDmr,
              summary=geneDistListByDmrSummary));
}

reqlengths = c(1000, 10000, 50000, 100000, 250000, 500000, 1e6, 1e6+5e5, 2e6, 2e6+5e5, 3e6, 4e6, 5e6, 7e6+5e5, 10e6);
dmrsDistToGenes_different_req_scaffoldLengths = list();
for (i in reqlengths) {
  cat(i,'...')
  dmrsDistToGenes_different_req_scaffoldLengths[[as.character(i)]] = .getDmrGeneDistStatsForCertainLengthScaffolds(dmrsGR,
                                                                                                                   geneGR,
                                                                                                                   slGeneStats,
                                                                                                                   i);
}; rm(i);

# for scaffolds of a certain length, how many genes are within a certain distance?
x = rbind(sapply(dmrsDistToGenes_different_req_scaffoldLengths[1:13], 
                 function(ff) mean(sapply(ff$distlist, function(f) sum(f <= 1000)))),
          sapply(dmrsDistToGenes_different_req_scaffoldLengths[1:13], 
                 function(ff) mean(sapply(ff$distlist, function(f) sum(f <= 5000)))),
          sapply(dmrsDistToGenes_different_req_scaffoldLengths[1:13], 
                 function(ff) mean(sapply(ff$distlist, function(f) sum(f <= 10000)))),
          sapply(dmrsDistToGenes_different_req_scaffoldLengths[1:13], 
                 function(ff) mean(sapply(ff$distlist, function(f) sum(f <= 20000)))),
          sapply(dmrsDistToGenes_different_req_scaffoldLengths[1:13], 
                 function(ff) mean(sapply(ff$distlist, function(f) sum(f <= 50000)))),
          sapply(dmrsDistToGenes_different_req_scaffoldLengths[1:13], 
                 function(ff) mean(sapply(ff$distlist, function(f) sum(f <= 100000)))),
          sapply(dmrsDistToGenes_different_req_scaffoldLengths[1:13], 
                 function(ff) mean(sapply(ff$distlist, function(f) sum(f <= 500000)))),
          sapply(dmrsDistToGenes_different_req_scaffoldLengths[1:13], 
                 function(ff) mean(sapply(ff$distlist, function(f) sum(f <= 1e6)))))

thiswindow = 50000
xl=(sapply(dmrsDistToGenes_different_req_scaffoldLengths[1:13],
           function(ff) (sapply(ff$distlist, function(f) sum(f <= thiswindow)))));
par(mfrow=c(3,4), oma=c(0,0,2,0));
for(i in 2:13) {
  verboseBoxplot(c(xl$'1000', 
                   xl[[i]]), 
                 c(rep('1000',length(xl$'1000')), 
                   rep(names(xl)[i],length(xl[[i]]))), 
                 ylab='', xlab='', frame.plot=F, col='grey')
}; rm(i);
title(paste0('genes within ',thiswindow,'bp'), outer=T)

dmrsDistToGenes_1mb_scaffoldLengths_thresh = list();
test_dists = seq(0, 1e6, 5000);
for (d in test_dists) {
  dmrsDistToGenes_1mb_scaffoldLengths_thresh[[as.character(d)]] = lapply(dmrsDistToGenes_1mb_scaffoldLengths, 
                                                                         function(f) f[f <= d]);
}; rm(d);
dmrsDistToGenes_1mb_scaffoldLengths_thresh_num_dmrs = sapply(dmrsDistToGenes_1mb_scaffoldLengths_thresh, 
                                         function(f) sum(sapply(f, length) > 0));

dmrsDistToGenes_1mb_scaffoldLengths_thresh_avg_genes = sapply(dmrsDistToGenes_1mb_scaffoldLengths_thresh, 
                                          function(f) mean(sapply(f, length)));

dmrsDistToGenes_1mb_scaffoldLengths_thresh_avg_dist = sapply(dmrsDistToGenes_1mb_scaffoldLengths_thresh, 
                                         function(f) mean(sapply(f, mean), na.rm=T));

par(mfrow=c(3,2), oma=c(0,0,2,0));
cex=1; cex.axis=1.5; cex.lab=1.5; cex.main=1.5;
maxind = which(dmrsDistToGenes_1mb_scaffoldLengths_thresh_num_dmrs == max(dmrsDistToGenes_1mb_scaffoldLengths_thresh_num_dmrs))[1]
plot(as.numeric(names(dmrsDistToGenes_1mb_scaffoldLengths_thresh_num_dmrs)),
     dmrsDistToGenes_1mb_scaffoldLengths_thresh_num_dmrs,
     type='b', frame.plot=F, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main,
     ylab='number of dmrs', xlab='distance with at least 1 gene');
abline(h=max(dmrsDistToGenes_1mb_scaffoldLengths_thresh_num_dmrs), 
       v=as.numeric(names(dmrsDistToGenes_1mb_scaffoldLengths_thresh_num_dmrs)[maxind]), col='red');
plot(as.numeric(names(dmrsDistToGenes_1mb_scaffoldLengths_thresh_num_dmrs))[1:(maxind+1)],
     dmrsDistToGenes_1mb_scaffoldLengths_thresh_num_dmrs[1:(maxind+1)],
     type='b', frame.plot=F, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main,
     ylab='number of dmrs', xlab='distance with at least 1 gene');
abline(h=max(dmrsDistToGenes_1mb_scaffoldLengths_thresh_num_dmrs), 
       v=as.numeric(names(dmrsDistToGenes_1mb_scaffoldLengths_thresh_num_dmrs)[maxind]), col='red');
verboseScatterplot(as.numeric(names(dmrsDistToGenes_1mb_scaffoldLengths_thresh_avg_genes)), 
                   dmrsDistToGenes_1mb_scaffoldLengths_thresh_avg_genes, 
                   abline=T, abline.col='red', abline.lty='dashed', 
                   type='b', frame.plot=F, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main,
                   xlab='distance', ylab='avg num genes');
abline(v=as.numeric(names(dmrsDistToGenes_1mb_scaffoldLengths_thresh_num_dmrs)[maxind]), col='red');
verboseScatterplot(as.numeric(names(dmrsDistToGenes_1mb_scaffoldLengths_thresh_avg_genes))[1:(maxind+1)], 
                   dmrsDistToGenes_1mb_scaffoldLengths_thresh_avg_genes[1:(maxind+1)], 
                   abline=T, abline.col='red', abline.lty='dashed', 
                   type='b', frame.plot=F, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main,
                   xlab='distance', ylab='avg num genes');
abline(v=as.numeric(names(dmrsDistToGenes_1mb_scaffoldLengths_thresh_num_dmrs)[maxind]), col='red');
verboseScatterplot(as.numeric(names(dmrsDistToGenes_1mb_scaffoldLengths_thresh_avg_dist)), 
                   dmrsDistToGenes_1mb_scaffoldLengths_thresh_avg_dist, 
                   abline=T, abline.col='red', abline.lty='dashed', 
                   type='b', frame.plot=F, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main,
                   xlab='distance', ylab='avg distance to genes');
abline(v=as.numeric(names(dmrsDistToGenes_1mb_scaffoldLengths_thresh_num_dmrs)[maxind]), col='red');
verboseScatterplot(as.numeric(names(dmrsDistToGenes_1mb_scaffoldLengths_thresh_avg_dist))[1:(maxind+1)], 
                   dmrsDistToGenes_1mb_scaffoldLengths_thresh_avg_dist[1:(maxind+1)], 
                   abline=T, abline.col='red', abline.lty='dashed', 
                   type='b', frame.plot=F, cex=cex, cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main,
                   xlab='distance', ylab='avg distance to genes');
abline(v=as.numeric(names(dmrsDistToGenes_1mb_scaffoldLengths_thresh_num_dmrs)[maxind]), col='red');
title('only dmrs on scaffolds at least 1mb long', outer=T)

# find the nearest dmr to every gene
genesNNdmr = as.data.frame(distanceToNearest(geneGR, dmrsGR));
genesNNdmr$queryHits = mcols(geneGR)[genesNNdmr$queryHits,1];
genesNNdmr$subjectHits = names(dmrsGR)[genesNNdmr$subjectHits];

genesNNdmr=as.data.frame(cbind(genesNNdmr, an$gffGenesDF[match(genesNNdmr$queryHits, an$gffGenesDF$geneSym), ]))

genesNNdmr1mb = subset(genesNNdmr, distance <= 1e6);
genesNNdmr1mb_byDmrList = split(genesNNdmr1mb, genesNNdmr1mb$subjectHits);
genesNNdmr1mb_dmrdf = data.frame(numgenes=sapply(genesNNdmr1mb_byDmrList,nrow),
                                 meandist=sapply(genesNNdmr1mb_byDmrList, function(f) mean(f$distance)));

genesNNdmr1mb_avgDistanceByScaffold = sapply(split(genesNNdmr1mb, 
                                                   sapply(strsplit(genesNNdmr1mb$subjectHits,':'), 
                                                          function(f) f[1])), 
                                             function(ff) mean(ff$distance));
genesNNdmr_avgDistanceByScaffold = sapply(split(genesNNdmr, 
                                                sapply(strsplit(genesNNdmr$subjectHits,':'), 
                                                       function(f) f[1])), 
                                          function(ff) mean(ff$distance));

slGeneStats = as.data.frame(cbind(slGeneStats, 
                                  avgDistGeneToNearestDmr=genesNNdmr_avgDistanceByScaffold[match(rownames(slGeneStats), 
                                                                                                 names(genesNNdmr_avgDistanceByScaffold))],
                                  avgDistGeneToNearestDmr1mb=genesNNdmr1mb_avgDistanceByScaffold[match(rownames(slGeneStats), 
                                                                                                       names(genesNNdmr1mb_avgDistanceByScaffold))]));

scaffolds_with_dmrs_but_no_genes = unlist(sapply(apply(slGeneStats, 2, 
                                                       function(f) which(is.infinite(f))), 
                                                 names));
slGeneStats$dmrsPerGene[rownames(slGeneStats) %in% scaffolds_with_dmrs_but_no_genes] = NA;

x = .order2(slGeneStats, 'length')


dmrrows = match(rownames(genesNNdmr1mb_dmrdf), rownames(dmrs));
genesNNdmr1mb_dmrdf = as.data.frame(cbind(genesNNdmr1mb_dmrdf, 
                                          dmrs[dmrrows, 
                                               match(c('n','width','invdensity','areaStat','meanDiff','direction'), 
                                                     names(dmrs))],
                                          abs.invdensity=abs(dmrs$invdensity[dmrrows]),
                                          abs.areaStat=abs(dmrs$areaStat[dmrrows]),
                                          abs.meanDiff=abs(dmrs$meanDiff[dmrrows]),
                                          fc=log2(dmrs$group2.mean/dmrs$group1.mean)[dmrrows],
                                          abs.fc=abs(log2(dmrs$group2.mean/dmrs$group1.mean)[dmrrows]),
                                          dmrrank=na.omit(dmrrows))); rm(dmrrows);

dmrsdf = as.data.frame(cbind(dmrs[, match(c('n','width','invdensity','areaStat','meanDiff','direction'),names(dmrs))],
                             abs.invdensity=abs(dmrs$invdensity),
                             abs.areaStat=abs(dmrs$areaStat),
                             abs.meanDiff=abs(dmrs$meanDiff),
                             fc=log2(dmrs$group2.mean/dmrs$group1.mean),
                             abs.fc=abs(log2(dmrs$group2.mean/dmrs$group1.mean)),
                             dmrrank=1:nrow(dmrs))); 

summary(genesNNdmr1mb_dmrdf$numgenes)
# Min. 1st Qu. Median  Mean  3rd Qu. Max. 
#    1       6     15  19.2      26  102

less15 = genesNNdmr1mb_dmrdf$numgenes < 15;
hypo = genesNNdmr1mb_dmrdf$direction == 'hypo';
genesNNdmr1mb_dmrdf15numgene_stats = as.data.frame(rbind(apply(genesNNdmr1mb_dmrdf[-8], 2, median), 
                                                         apply(genesNNdmr1mb_dmrdf[hypo, -8], 2, median),
                                                         apply(genesNNdmr1mb_dmrdf[!hypo, -8], 2, median),
                                                         apply(genesNNdmr1mb_dmrdf[less15, -8], 2, median),
                                                         apply(genesNNdmr1mb_dmrdf[!less15, -8], 2, median),
                                                         apply(genesNNdmr1mb_dmrdf[hypo & less15, -8], 2, median),
                                                         apply(genesNNdmr1mb_dmrdf[!hypo & less15, -8], 2, median),
                                                         apply(genesNNdmr1mb_dmrdf[hypo & !less15, -8], 2, median),
                                                         apply(genesNNdmr1mb_dmrdf[!hypo & !less15, -8], 2, median)));
rownames(genesNNdmr1mb_dmrdf15numgene_stats) = c('all', 'hypo', 'hyper', '<15', '>=15',
                                                 'hypo_<15', 'hyper_<15',
                                                 'hypo_>=15', 'hyper_>=15');






genesNNdmr50kb = subset(genesNNdmr, distance <= 50000);
genesNNdmr50kb_byDmrList = split(genesNNdmr50kb, genesNNdmr50kb$subjectHits);
genesNNdmr50kb_dmrdf = data.frame(numgenes=sapply(genesNNdmr50kb_byDmrList,nrow),
                                 meandist=sapply(genesNNdmr50kb_byDmrList, function(f) mean(f$distance)));

dmrrows = match(rownames(genesNNdmr50kb_dmrdf), rownames(dmrs));
genesNNdmr50kb_dmrdf = as.data.frame(cbind(genesNNdmr50kb_dmrdf, 
                                          dmrs[dmrrows, 
                                               match(c('n','width','invdensity','areaStat','meanDiff','direction'), 
                                                     names(dmrs))],
                                          abs.invdensity=abs(dmrs$invdensity[dmrrows]),
                                          abs.areaStat=abs(dmrs$areaStat[dmrrows]),
                                          abs.meanDiff=abs(dmrs$meanDiff[dmrrows]),
                                          fc=log2(dmrs$group2.mean/dmrs$group1.mean)[dmrrows],
                                          abs.fc=abs(log2(dmrs$group2.mean/dmrs$group1.mean)[dmrrows]),
                                          dmrrank=na.omit(dmrrows))); rm(dmrrows);



genesNNdmr100kb = subset(genesNNdmr, distance <= 100000);
genesNNdmr100kb_byDmrList = split(genesNNdmr100kb, genesNNdmr100kb$subjectHits);
genesNNdmr100kb_dmrdf = data.frame(numgenes=sapply(genesNNdmr100kb_byDmrList,nrow),
                                   meandist=sapply(genesNNdmr100kb_byDmrList, function(f) mean(f$distance)));

dmrrows = match(rownames(genesNNdmr100kb_dmrdf), rownames(dmrs));
genesNNdmr100kb_dmrdf = as.data.frame(cbind(genesNNdmr100kb_dmrdf, 
                                            dmrs[dmrrows, 
                                                 match(c('n','width','invdensity','areaStat','meanDiff','direction'), 
                                                       names(dmrs))],
                                            abs.invdensity=abs(dmrs$invdensity[dmrrows]),
                                            abs.areaStat=abs(dmrs$areaStat[dmrrows]),
                                            abs.meanDiff=abs(dmrs$meanDiff[dmrrows]),
                                            fc=log2(dmrs$group2.mean/dmrs$group1.mean)[dmrrows],
                                            abs.fc=abs(log2(dmrs$group2.mean/dmrs$group1.mean)[dmrrows]),
                                            dmrrank=na.omit(dmrrows))); rm(dmrrows);


### plots ####
###



summary(genesNNdmr1mb_dmrdf$numgenes)
# Min. 1st Qu. Median  Mean  3rd Qu. Max. 
#    1       6     15  19.2      26  102 

par(oma=c(0,0,2,0)); toplot = genesNNdmr1mb_dmrdf[genesNNdmr1mb_dmrdf$numgenes < 15, ];
.verboseBoxplotColumns(toplot[-8], 
                       toplot$direction=='hypo',
                       names=c('ND-hyper','D-hyper'), frame.plot=F, col=c('lightblue','greenyellow'))
title('all dmrs within 1mb <=15 genes', outer=T);

par(oma=c(0,0,2,0)); toplot = genesNNdmr1mb_dmrdf[genesNNdmr1mb_dmrdf$numgenes >= 15, ];
.verboseBoxplotColumns(toplot[-8], 
                       toplot$direction=='hypo',
                       names=c('ND-hyper','D-hyper'), frame.plot=F, col=c('lightblue','greenyellow'))
title('all dmrs within 1mb >15 genes', outer=T);

par(oma=c(0,0,2,0)); toplot = genesNNdmr1mb_dmrdf;
.verboseBoxplotColumns(toplot[-8], 
                       toplot$numgenes < 15,
                       names=c('>=15 genes','<15 genes'), 
                       frame.plot=F, col=c('lightblue','greenyellow'))
title('dmrs within 1mb', outer=T);

par(oma=c(0,0,2,0)); toplot = genesNNdmr1mb_dmrdf[genesNNdmr1mb_dmrdf$direction == 'hypo', ];
.verboseBoxplotColumns(toplot[-8], 
                       toplot$numgenes < 15,
                       names=c('>=15 genes','<15 genes'), 
                       frame.plot=F, col=c('lightblue','greenyellow'))
title('hypo dmrs within 1mb', outer=T);

par(oma=c(0,0,2,0)); toplot = genesNNdmr1mb_dmrdf[genesNNdmr1mb_dmrdf$direction == 'hyper', ];
.verboseBoxplotColumns(toplot[-8], 
                       toplot$numgenes < 15,
                       names=c('>=15 genes','<15 genes'), 
                       frame.plot=F, col=c('lightblue','greenyellow'))
title('hyper dmrs within 1mb', outer=T);





toplotname = 'genesNNdmr1mb_dmrdf';

toplot = get(toplotname);
par(mfrow=c(3,2)); breaks=25; col='grey'; border='grey'
hist(toplot$numgenes, breaks=breaks, col=col, border=border, xlab=paste0(toplotname,'$numgenes'),
     main=paste0(.quantile.tails(toplot$numgenes),collapse=', '));
abline(v=.quantile.tails(toplot$numgenes), col='red')
hist(toplot$meandist, breaks=breaks, col=col, border=border, xlab=paste0(toplotname,'$meandist'),
     main=paste0(round(.quantile.tails(toplot$meandist)),collapse=', '));
abline(v=.quantile.tails(toplot$meandist), col='red');

toplot = subset(toplot, numgenes > 1);
hist(toplot$numgenes, breaks=breaks, col=col, border=border, xlab=paste0(toplotname,'$numgenes'),
     main=paste0(.quantile.tails(toplot$numgenes),collapse=', '));
abline(v=.quantile.tails(toplot$numgenes), col='red')
hist(toplot$meandist, breaks=breaks, col=col, border=border, xlab=paste0(toplotname,'$meandist'),
     main=paste0(round(.quantile.tails(toplot$meandist)),collapse=', '));
abline(v=.quantile.tails(toplot$meandist), col='red');

toplot = subset(toplot, numgenes > 15);
hist(toplot$numgenes, breaks=breaks, col=col, border=border, xlab=paste0(toplotname,'$numgenes'),
     main=paste0(.quantile.tails(toplot$numgenes),collapse=', '));
abline(v=.quantile.tails(toplot$numgenes), col='red')
hist(toplot$meandist, breaks=breaks, col=col, border=border, xlab=paste0(toplotname,'$meandist'),
     main=paste0(round(.quantile.tails(toplot$meandist)),collapse=', '));
abline(v=.quantile.tails(toplot$meandist), col='red');

