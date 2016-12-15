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

# figure 5a, will only use one of 3 plots for main figure, others can go in supplements
par(mfrow=c(1,3));
ylim=c(0,4000);
col=c('grey','lightblue'); 
border=c('darkgrey','blue');
names = c('0 dmrs', '>0 dmrs');
subgenes = names(dmrsOVbyGene$gene);
thisres = resexpr;
thisfactor = thisres$dmr.num > 0 & rownames(thisres) %in% subgenes;
verboseBoxplot(thisres$baseMean, thisfactor, 
               ylab='DESeq2 baseMean expression', xlab='', 
               names=names, main='all expressed genes\n',
               frame.plot=F, ylim=ylim, col=col, border=border);
thisres = resexprSmoothed;
thisfactor = thisres$dmr.num > 0 & rownames(thisres) %in% subgenes;
verboseBoxplot(thisres$baseMean, thisfactor,  ylab='', xlab='',
               names=names, main='expressed genes on smoothed scaffolds\n',
               frame.plot=F, ylim=ylim, col=col, border=border);
thisres = resexprSmoothedDmrs;
thisfactor = thisres$dmr.num > 0 & rownames(thisres) %in% subgenes;
verboseBoxplot(thisres$baseMean, thisfactor, ylab='', xlab='',
               names=names, main='expressed genes on dmr scaffolds\n',
               frame.plot=F, ylim=ylim, col=col, border=border);

# figure 5b
par(mfrow=c(1,3));
ylim=c(0,.3);
verboseBoxplot(abs(resexpr$log2FoldChange), resexpr$dmr.num > 0, 
               ylab='abs(DESeq2 log2FoldChange)', xlab='',
               names=c('0 dmrs', '>0 dmrs'), 
               frame.plot=F, ylim=ylim,
               col=col, border=border,
               main='all expressed genes\n');
verboseBoxplot(abs(resexprSmoothed$log2FoldChange), resexprSmoothed$dmr.num > 0, 
               ylab='', xlab='',
               names=c('0 dmrs', '>0 dmrs'), 
               frame.plot=F, ylim=ylim,
               col=col, border=border,
               main='all expressed genes on smoothed scaffolds\n');
verboseBoxplot(abs(resexprSmoothedDmrs$log2FoldChange), resexprSmoothedDmrs$dmr.num > 0, 
               ylab='', xlab='',
               names=c('0 dmrs', '>0 dmrs'), 
               frame.plot=F, ylim=ylim,
               col=col, border=border,
               main='all expressed genes on dmr scaffolds\n');

# figure 5c
par(mfrow=c(1,3));
ylim=c(0,2.2); cex=1.8; col='black';
verboseScatterplot(resexpr$baseMean, abs(resexpr$log2FoldChange), 
                   abline=T, abline.col='red', frame.plot=F, ylim=ylim,
                   ylab='abs(DESeq2 log2FoldChange)', 
                   xlab='DESeq2 baseMean expression', 
                   bg=numbers2colors(-as.numeric(resexpr$dmr.num!=0)), cex=cex,
                   col=col,pch=21, corOptions="method='s'",
                   main='all expressed genes\n');
verboseScatterplot(resexprSmoothed$baseMean, abs(resexprSmoothed$log2FoldChange), 
                   abline=T, abline.col='red', frame.plot=F, ylim=ylim,
                   ylab='', 
                   xlab='', 
                   bg=numbers2colors(-as.numeric(resexprSmoothed$dmr.num!=0)), cex=cex,
                   col=col,pch=21, corOptions="method='s'",
                   main='all expressed genes on smoothed scaffolds\n');
verboseScatterplot(resexprSmoothedDmrs$baseMean, abs(resexprSmoothedDmrs$log2FoldChange), 
                   abline=T, abline.col='red', frame.plot=F, ylim=ylim,
                   ylab='', 
                   xlab='', 
                   bg=numbers2colors(-as.numeric(resexprSmoothedDmrs$dmr.num!=0)), cex=cex,
                   col=col,pch=21, corOptions="method='s'",
                   main='all expressed genes on dmr scaffolds\n');

# figure 5d
ylim=c(0,.3)
thisres = subset(resexprSmoothedDmrs, dmr.num > 0);
thisres = thisres[rownames(thisres) %in% names(dmrsOVbyGene$gene), ]
verboseScatterplot(thisres$baseMean, abs(thisres$log2FoldChange), 
                   abline=T, abline.col='red', frame.plot=F, ylim=ylim,
                   ylab='abs(DESeq2 log2FoldChange)', 
                   xlab='DESeq2 baseMean expression', 
                   bg=numbers2colors(-as.numeric(thisres$dmr.num)), cex=cex,
                   col=col,pch=21, corOptions="method='s'",
                   main='dmr-OV genes\n');

# figure ???
# ----------------------
dmrres = subset(resexprSmoothedDmrs, dmr.num>0)

par(mfrow=c(3,6));
col = c('grey','lightblue')
border = c('darkgrey','blue')
for (i in 25:30) {
  if (i==25) {
    ylab='DESeq2 log2FoldChange'
  } else {
    ylab=''
  }
  verboseBoxplot((dmrres$log2FoldChange), dmrres[,i], 
                 main=paste0(names(dmrres)[i],'\n'),
                 ylim=c(-.3,.3), ylab=ylab, xlab='', frame.plot=F,
                 col=col, border=border); 
  abline(h=0,col='red')
}; rm(i)
for (i in 25:30) {
  if (i==25) {
    ylab='abs(DESeq2 log2FoldChange)'
  } else {
    ylab=''
  }
  verboseBoxplot(abs(dmrres$log2FoldChange), dmrres[,i], 
                 main=paste0(names(dmrres)[i],'\n'),
                 ylim=c(0,.3), ylab=ylab, xlab='', frame.plot=F,
                 col=col, border=border); 
}; rm(i);
for (i in 25:30) {
  if (i==25) {
    ylab='DESeq2 baseMean expression'
  } else {
    ylab=''
  }
  verboseBoxplot((dmrres$baseMean), dmrres[,i], 
                 main=paste0(names(dmrres)[i],'\n'),
                 ylim=c(0,5000), ylab=ylab, xlab='', frame.plot=F,
                 col=col, border=border); 
}; rm(i)

# plot expr vs dmr fold changes for genes based on feature type and expression level
# ----------------------
par(mfrow=c(3,6));
for (i in 25:30) {
  subgenes = as.logical(dmrres[,i]);
  bgcol = numbers2colors(dmrres$dmr.num[subgenes]);
  verboseScatterplot(dmrres$dmrfc[subgenes], dmrres$log2FoldChange[subgenes],
                     abline=T, abline.col='red', frame.plot=F, 
                     xlab='dmrfc', ylab='exprfc', main=paste0('all dmr genes: ',names(dmrres)[i],'\n'),
                     bg=bgcol, pch=21, col='black');
  abline(h=0, v=0, col='grey');
}; rm(i);
thisres = subset(dmrres, baseMean >= median(baseMean));
for (i in 25:30) {
  subgenes = as.logical(thisres[,i]);
  bgcol = numbers2colors(thisres$dmr.num[subgenes]);
  verboseScatterplot(thisres$dmrfc[subgenes], thisres$log2FoldChange[subgenes],
                     abline=T, abline.col='red', frame.plot=F, 
                     xlab='dmrfc', ylab='exprfc', main=paste0('high expr: ', names(thisres)[i],'\n'),
                     bg=bgcol, pch=21, col='black');
  abline(h=0, v=0, col='grey');
}; rm(i);
thisres = subset(dmrres, baseMean < median(baseMean));
for (i in 25:30) {
  subgenes = as.logical(thisres[,i]);
  bgcol = numbers2colors(thisres$dmr.num[subgenes]);
  verboseScatterplot(thisres$dmrfc[subgenes], thisres$log2FoldChange[subgenes],
                     abline=T, abline.col='red', frame.plot=F, 
                     xlab='dmrfc', ylab='exprfc', main=paste0('low expr: ',names(thisres)[i],'\n'),
                     bg=bgcol, pch=21, col='black');
  abline(h=0, v=0, col='grey');
}; rm(i);




# ----------------------

## supplemental figure?
par(mfrow=c(1,2))
col=c('grey','lightblue'); 
border=c('darkgrey','blue');
ylim=c(0,1)
thisres = subset(resexprSmoothedDmrs, dmr.num > 0);
#thisres = thisres[rownames(thisres) %in% names(dmrsOVbyGene$gene), ]
verboseBoxplot(abs(thisres$log2FoldChange), thisres$dmrfc > 0, 
               ylab='abs(DESeq2 log2FoldChange)', xlab='',
               names=c('hyper', 'hypo'), 
               frame.plot=F, ylim=ylim,
               col=col, border=border,
               main='dmr-OV genes\n');
ylim=c(0,4000)
verboseBoxplot(thisres$baseMean, thisres$dmrfc > 0, 
               ylab='DESeq2 baseMean expression', xlab='',
               names=c('hyper', 'hypo'), 
               frame.plot=F, ylim=ylim,
               col=col, border=border,
               main='dmr-OV genes\n');

# figure 5e
par(mfrow=c(1,3));
ylim=c(-.3,.3);
verboseBoxplot((resexpr$log2FoldChange), resexpr$dmr.num > 0, 
               ylab='DESeq2 log2FoldChange', xlab='',
               names=c('0 dmrs', '>0 dmrs'), 
               frame.plot=F, ylim=ylim,
               col=col, border=border,
               main='all expressed genes\n');abline(h=0,col='red')
verboseBoxplot((resexprSmoothed$log2FoldChange), resexprSmoothed$dmr.num > 0, 
               ylab='', xlab='',
               names=c('0 dmrs', '>0 dmrs'), 
               frame.plot=F, ylim=ylim,
               col=col, border=border,
               main='all expressed genes on smoothed scaffolds\n');abline(h=0,col='red')
verboseBoxplot((resexprSmoothedDmrs$log2FoldChange), resexprSmoothedDmrs$dmr.num > 0, 
               ylab='', xlab='',
               names=c('0 dmrs', '>0 dmrs'), 
               frame.plot=F, ylim=ylim,
               col=col, border=border,
               main='all expressed genes on dmr scaffolds\n');abline(h=0,col='red')

# figure 6a



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