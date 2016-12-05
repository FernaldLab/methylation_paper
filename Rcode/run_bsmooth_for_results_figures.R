# figure 2a
req_scaffold_lengths = list(0,50000,5e5,1e6,2e6,3e6,4e6,5e6)
toplot = as.data.frame(sapply(req_scaffold_lengths, 
                              function(f) apply(dmrsDistToGenes_thresh_countDf[dmrs_scaffold_lengths >= f,], 2, 
                                                mean)));
names(toplot) = req_scaffold_lengths;

cols = c('black','powderblue','lightblue','skyblue','cornflowerblue','steelblue','royalblue','blue','blue2');
par(mfrow=c(2,1))
for (i in 1:ncol(toplot)) {
  if (i < ncol(toplot)) {
    xlab=''; ylab=''; axes=F;
  } else {
    xlab='window size'; ylab='num. genes'; axes=T;
  }
  plot(as.numeric(rownames(toplot)),
       toplot[,i],
       ylim=c(0,max(toplot)+10),
       type='b',
       pch=19,
       frame.plot=F,
       col=cols[i],
       xlab=xlab,ylab=ylab,axes=axes);
  par(new=T)
}; rm(i); par(new=F);
abline(h=10,v=1e5,col='darkgrey')
toplot = toplot[1:25, ]
for (i in 1:ncol(toplot)) {
  if (i < ncol(toplot)) {
    xlab=''; ylab=''; axes=F;
  } else {
    xlab='window size'; ylab='num. genes'; axes=T;
  }
  plot(as.numeric(rownames(toplot)),
       toplot[,i],
       ylim=c(0,10),
       type='b',
       pch=19,
       frame.plot=F,
       col=cols[i],
       xlab=xlab,ylab=ylab,axes=axes);
  par(new=T)
}; rm(i); par(new=F);
abline(v=55000,h=c(4.2581100, 4.9113924),col='red')
barplot(rep(.1,8), col=cols, space=0, border=NA, horiz=T, names.arg=req_scaffold_lengths,las=2,axes=F)



##
req_scaf_len = 1e6;
par(mfrow=c(4,9))
for (win_to_test in as.numeric(names(dmrsDistToGenes_thresh_countDf))) {
  rows = dmrs_scaffold_lengths >= req_scaf_len;
  wincol = which(names(dmrsDistToGenes_thresh_countDf) == as.character(win_to_test));
  vals = c(dmrsDistToGenes_thresh_countDf[,wincol], 
           dmrsDistToGenes_thresh_countDf[rows,wincol]);
  facts = as.factor(c(rep('all', nrow(dmrsDistToGenes_thresh_countDf)),
                      rep(paste0('>',as.character(req_scaf_len)), sum(rows))));
  verboseBoxplot(vals, facts, 
                 xlab='', ylab='num.genes', main=paste0(as.character(win_to_test),'\n'), 
                 col='grey', border='darkgrey', frame.plot=F);
}

win_to_test = 50000;
wincol = which(names(dmrsDistToGenes_thresh_countDf) == as.character(win_to_test));
par(mfrow=c(2,4));
for (req_scaf_len in req_scaffold_lengths) {
  rows = dmrs_scaffold_lengths >= req_scaf_len;
  vals = c(dmrsDistToGenes_thresh_countDf[,wincol], 
           dmrsDistToGenes_thresh_countDf[rows,wincol]);
  facts = as.factor(c(rep('all', nrow(dmrsDistToGenes_thresh_countDf)),
                      rep(paste0('>',as.character(req_scaf_len)), sum(rows))));
  verboseBoxplot(vals, facts, 
                 xlab='', ylab='num.genes', main=paste0(as.character(win_to_test),'\n'), 
                 col='grey', border='darkgrey', frame.plot=F);
}


##

# figure 2b
tmp = dmrsDistToGenes_thresh_summaryDf[1:25,];
tmpdf = as.data.frame(matrix(nrow=nrow(tmp),ncol=3,
                             dimnames=list(rownames(tmp),
                                           c('actual','shuff','pval'))))
tmppvals = c();
for (i in 1:nrow(tmp)) {
  cat(rownames(tmp)[i],'...')
  thisactual = dmrsDistToGenes_thresh_summaryDf$count[rownames(dmrsDistToGenes_thresh_summaryDf)==rownames(tmp)[i]];
  thisshuff = sapply(shuffDistToGenesAcrossRuns, 
                     function(ff) mean(sapply(ff, 
                                              function(f) sum(f <= as.numeric(rownames(tmp)[i])))));
  tmpdf[i,] = c(thisactual, mean(thisshuff), sum(thisshuff >= thisactual)/length(shuffDistToGenesAcrossRuns))
  tmppvals = c(tmppvals, sum(thisshuff >= thisactual)/length(shuffDistToGenesAcrossRuns))
}
rm(i,thisactual,thisshuff);


# figure 2c
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
dev.off();
toplot = dmrsDistToGenes_thresh_countDf_testStats[1:25,];
dodds = dodds[1:25];
colors=c('black', 'midnightblue','blue','skyblue','darkgreen','green','grey','purple','orange','pink', 'cyan','lightgreen')
for (i in c(1:3,5,7,9)) {
  if(names(toplot)[i] %in% c('n','width','invdensity','areaStat.abs')) {
    thiscol = 'darkgrey'
  } else if (names(toplot)[i]=='meanDiff.abs') {
    thiscol = 'purple'
  } else {
    thiscol = 'green'
  }
  #thiscol = colors[i]
  plot(as.numeric(rownames(toplot)), 
       -log10(toplot[,i]), 
       type='b', pch=19, col=thiscol, ylim=c(0,5), 
       ylab='', xlab='', frame.plot=F); 
  cat(paste0(names(toplot)[i], ':', thiscol,'\n'));par(new=T)
}; rm(i); par(new=T);
plot(as.numeric(rownames(toplot)),
     -log10(dodds),
     type='b', pch=19, col='darkgrey', ylim=c(0,5), lty='dashed',
     ylab='-log10(pval)',xlab='window size', frame.plot=F);
abline(h=-log10(c(.05,.01)), col='red');
abline(v=c(5000,20000,50000,100000), col='darkgrey');
par(new=F);
#
boxplot(-log10(toplot), notch=T, 
        col=c(rep('darkgrey',6),'purple','darkgrey','green'),
        border=c(rep('grey60',6),'purple3','grey60','darkgreen'), 
        frame.plot=F, ylab='-log10(pval)', ylim=c(0,5), pch=19)
abline(h=-log10(c(.05,.01)), col='red');
#
rows = dmrs_scaffold_lengths > 0e6;
.verboseScatterplotVecAgainstColumns(dmrs3[rows,1:9], 
                                     dmrsDistToGenes_thresh_countDf$'5000'[rows], 
                                     corOptions="method='s'", 
                                     pch=19, col=labels2colors(dmrs$direction[rows]), 
                                     xlab='num.genes');



# supplemental figure 1
dmrs2=as.data.frame(cbind(dmrs[,7:15], fc=log2(dmrs$group2.mean/dmrs$group1.mean)));
dmrs2 = as.data.frame(cbind(dmrs2[,1:4], 
                            areaStat.abs=abs(dmrs2$areaStat), 
                            meanDiff=dmrs2[,6], 
                            meanDiff.abs=abs(dmrs2$meanDiff), 
                            log2fc=dmrs2[,10], 
                            log2fc.abs=abs(dmrs2$fc)));
par(mfrow=c(2,5));
for (i in 1:ncol(dmrs2)) { 
  hist(dmrs2[,i], breaks=200, 
       col='grey', border='grey',
       ylab='', xlab='', 
       main=paste0(names(dmrs2)[i],
                   '\nmedian=',
                   signif(median(dmrs2[,i]), 2),
                   ', mean=', 
                   signif(mean(dmrs2[,i]), 2))); 
  abline(v=c(median(dmrs2[,i]), mean(dmrs2[,i])), col='red');     
}; rm(i)


par(mfrow=c(1,8));
for (i in 2:ncol(dmrs2)) {
  verboseScatterplot(dmrs2$n, dmrs2[,i], 
                     xlab='n', ylab=names(dmrs2)[i],
                     abline=T, frame.plot=F, pch=19, col=labels2colors(dmrs$direction))
}; rm(i);

par(mfrow=c(1,8));
for (i in c(3:ncol(dmrs2))) {
  verboseScatterplot(dmrs2$width, dmrs2[,i], 
                     xlab='width', ylab=names(dmrs2)[i],
                     abline=T, frame.plot=F, pch=19, col=labels2colors(dmrs$direction))
}; rm(i)

par(mfrow=c(1,8));
for (i in c(4:ncol(dmrs2))) {
  verboseScatterplot(dmrs2$invdensity, dmrs2[,i], 
                     xlab='invdensity', ylab=names(dmrs2)[i],
                     abline=T, frame.plot=F, pch=19, col=labels2colors(dmrs$direction))
}; rm(i)

par(mfrow=c(1,8));
for (i in c(6:ncol(dmrs2))) {
  verboseScatterplot(dmrs2$areaStat, dmrs2[,i], 
                     xlab='areaStat', ylab=names(dmrs2)[i],
                     abline=T, frame.plot=F, pch=19, col=labels2colors(dmrs$direction))
}; rm(i)

par(mfrow=c(1,8));
for (i in c(6:ncol(dmrs2))) {
  verboseScatterplot(dmrs2$areaStat.abs, dmrs2[,i], 
                     xlab='areaStat,abs', ylab=names(dmrs2)[i],
                     abline=T, frame.plot=F, pch=19, col=labels2colors(dmrs$direction))
}; rm(i)

par(mfrow=c(1,8));
for (i in c(8:ncol(dmrs2))) {
  verboseScatterplot(dmrs2$meanDiff, dmrs2[,i], 
                     xlab='meanDiff', ylab=names(dmrs2)[i],
                     abline=T, frame.plot=F, pch=19, col=labels2colors(dmrs$direction))
}; rm(i)

par(mfrow=c(1,8));
for (i in c(8:ncol(dmrs2))) {
  verboseScatterplot(dmrs2$meanDiff.abs, dmrs2[,i], 
                     xlab='meanDiff.abs', ylab=names(dmrs2)[i],
                     abline=T, frame.plot=F, pch=19, col=labels2colors(dmrs$direction))
}; rm(i)