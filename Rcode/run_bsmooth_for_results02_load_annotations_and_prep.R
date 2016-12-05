# clear workspace, set options and working directory, load functions
rm(list=ls()); options(stringsAsFactors=F); 
setwd('~/Documents/_BS-seq_analysis/');
library('GenomicRanges');
library('Rsamtools');
library('bsseq')
source('/Volumes/fishstudies-1/_code/run_bsmooth_functions.R');

###
# load dmrs
load('run_bsmooth_for_results01_load_common_dmrs_hc_and_null_filterWORKSPACE.RData');
dmrs = dmrs_filtered_null.5;
rm(list=ls()[ls()!='dmrs']);

# load workspace from 'annotations_prep_for_dmrs_wgcna.R'
# should load:
#  handannos, sl0, slGeneStats, zcne, abHsMap, an, an2, annoCombo, sl

load('~/Documents/_annotationsDec2015/WORKSPACE_annotations_prep_for_dmrs_wgcna.RData')
rm(handannos, sl0, zcne, an2, miGR, mrnaGR, sl); gc();
burtoniGois = as.vector(na.omit(annoCombo$new$gene[match(an$gois, annoCombo$new$Dbxref)]));
names(mcols(cdsGR))[1:2] = c('mcols.geneSym','mcols.CG');
names(mcols(exonGR))[1:2] = c('mcols.geneSym','mcols.CG');
names(mcols(intronGR))[1:2] = c('mcols.geneSym','mcols.id');
names(mcols(zGR))[1:2] = c('mcols.id','mcols.score');

names(geneGR) = mcols(geneGR)[,1];
names(brGR) = mcols(brGR)[,1];
names(br3primeGR) = mcols(br3primeGR)[,1];
names(intronGR) = mcols(intronGR)[,1];
names(exonGR) = mcols(exonGR)[,1];
names(cdsGR) = mcols(cdsGR)[,1];

# get %CG for basal regulatory regions
genomeFasta = FaFile('/Volumes/fishstudies-1/_Burtoni_genome_files/H_burtoni_v1.assembly.fa');
brSeqs = getSeq(genomeFasta, brGR);
brSeqsNucs = alphabetFrequency(brSeqs);
rownames(brSeqsNucs) = mcols(brGR)[,1];
brSeqsNucsCG = apply(brSeqsNucs, 1, 
                     function(f) sum(f[match(c('C','G'), colnames(brSeqsNucs))]) / sum(f) );
mcols(brGR) <- as.data.frame(cbind(as.data.frame(mcols(brGR)), CG=brSeqsNucsCG));
rm(brSeqs, brSeqsNucs, brSeqsNucsCG);

# load names of scaffolds that were smoothed
# comes from cpgcommon object in run_bsmooth_for_methods_June2016.R
# should load: scaffolds_smoothed
load('scaffold_names_input_to_BSmooth.RData')

# add dmr information to slGeneStats
dmrs_per_scaffold = table(dmrs$chr);
slGeneStats = as.data.frame(cbind(slGeneStats, 
                                  dmrs=dmrs_per_scaffold[match(rownames(slGeneStats),
                                                               names(dmrs_per_scaffold))]));
slGeneStats = as.data.frame(cbind(slGeneStats,
                                  smoothed=as.numeric(rownames(slGeneStats) %in% scaffolds_smoothed),
                                  avgDistBetweenDmrs=slGeneStats$length/slGeneStats$dmrs,
                                  dmrsPerMb=slGeneStats$dmrs/slGeneStats$length*1e6,
                                  dmrsPerGene=slGeneStats$dmrs/slGeneStats$numGenes));
slGeneStats$dmrs[is.na(slGeneStats$dmrs) & slGeneStats$smoothed] = 0;
rm(dmrs_per_scaffold,scaffolds_smoothed);

# add CG content
# get %CG for scaffolds
genomeSeqs = getSeq(genomeFasta);
genomeSeqsNucs = alphabetFrequency(genomeSeqs);
rownames(genomeSeqsNucs) = names(genomeSeqs);
genomeSeqsNucsCG = apply(genomeSeqsNucs, 1, 
                         function(f) sum(f[match(c('C','G'), colnames(genomeSeqsNucs))]) / sum(f) );
slGeneStats = as.data.frame(cbind(slGeneStats, 
                                  scaffoldGC=genomeSeqsNucsCG[match(rownames(slGeneStats), 
                                                                    names(genomeSeqsNucsCG))]));
rm(genomeSeqs, genomeSeqsNucs, genomeSeqsNucsCG);

slGeneStats_in70 = subset(slGeneStats, smoothed==1);

##########
#genes_pctCDS = sapply(geneGR_dmrScaffolds, function(f) sum(width(reduce(cdsGR[names(cdsGR)==names(f)]))) / width(f))
########
 
# build GRanges for dmrs
mcolCols = which(names(dmrs)=='n'):ncol(dmrs);
dmrsGR = GRanges(seqnames=dmrs$chr, 
                 ranges=IRanges(start=dmrs$start, end=dmrs$end), 
                 mcols=dmrs[,mcolCols]);
# get %CG for dmrs
dmrsSeqs = getSeq(genomeFasta, dmrsGR);
dmrsSeqsNucs = alphabetFrequency(dmrsSeqs);
rownames(dmrsSeqsNucs) = names(dmrsGR);
dmrsSeqsNucsCG = apply(dmrsSeqsNucs, 1, 
                       function(f) sum(f[match(c('C','G'), colnames(dmrsSeqsNucs))]) / sum(f) );
mcols(dmrsGR) <- as.data.frame(cbind(name=names(dmrsGR), as.data.frame(mcols(dmrsGR)), GC=dmrsSeqsNucsCG));
dmrs = as.data.frame(cbind(dmrs, GC=dmrsSeqsNucsCG))
rm(genomeFasta, dmrsSeqs, dmrsSeqsNucs, dmrsSeqsNucsCG, mcolCols);

save(list=ls(), file='run_bsmooth_for_results02_load_annotations_and_prepWORKSPACE.RData');

#########################
### plots

# histograms of slGeneStats columns
par(mfrow=c(2,7), oma=c(0,0,2,0)); breaks=50;
for (i in 1:ncol(slGeneStats)) {
  hist(slGeneStats[,i], breaks=breaks, col='grey', border='darkgrey', main=names(slGeneStats)[i], xlab='');
  abline(v=quantile(slGeneStats[,i], c(.05,.5,.95), na.rm=T), col='red')
}; rm(i,breaks);
title('all scaffolds\nred lines = 5%, 50%, 95%', outer=T);

# histograms of slGeneStats_in70 columns
par(mfrow=c(2,7), oma=c(0,0,2,0)); breaks=50;
for (i in 1:ncol(slGeneStats_in70)) {
  hist(slGeneStats_in70[,i], breaks=breaks, col='grey', border='darkgrey', main=names(slGeneStats_in70)[i], xlab='');
  abline(v=quantile(slGeneStats_in70[,i], c(.05,.5,.95), na.rm=T), col='red')
}; rm(i);
title('all scaffolds input to BSmooth with >=70 CpGs \nred lines = 5%, 50%, 95%', outer=T);

# compare dmr-containing scaffolds to the rest
par(oma=c(0,0,2,0));
.verboseBoxplotColumns(slGeneStats_in70[,-grep('dmr|smoothed',names(slGeneStats_in70),ignore.case=T)], 
                       slGeneStats_in70$dmrs > 0,
                       frame.plot=F, col='grey', border='darkgrey',
                       names=c('0 dmrs', '>0 dmrs')
                       )
title('all scaffolds input to BSmooth with >=70 CpGs', outer=T);

# rbind(apply(slGeneStats_in70[slGeneStats_in70$dmrs > 0,], 2, median, na.rm=T), 
#       apply(slGeneStats_in70[slGeneStats_in70$dmrs == 0,], 2, median, na.rm=T));
#        length numGenes avgGeneWidth avgGeneGC avgNumIsoforms pctBpGene genesPerMb dmrs smoothed avgDistBetweenDmrs dmrsPerMb dmrsPerGene scaffoldGC
# [1,] 965616.0       31      19962.6   0.39685          1.901  0.020745   31.45046    1        1           546059.8  1.831301  0.05882353  0.3599692
# [2,]  72252.5        2      11345.0   0.38660          1.250  0.083870   33.71146    0        1                 NA        NA          NA  0.3202028

# rbind(apply(slGeneStats_in70[slGeneStats_in70$dmrs > 0,], 2, mean, na.rm=T), 
#       apply(slGeneStats_in70[slGeneStats_in70$dmrs == 0,], 2, mean, na.rm=T));
#         length numGenes avgGeneWidth avgGeneGC avgNumIsoforms  pctBpGene genesPerMb     dmrs smoothed avgDistBetweenDmrs dmrsPerMb dmrsPerGene scaffoldGC
# [1,] 1325566.2 43.20299      24077.2 0.3957894       1.922188 0.04689187   34.01965 2.116418        1             741025  4.360121         Inf  0.3466554
# [2,]  202384.1  5.85023      17331.8 0.3827374       1.496006 0.17246235   48.28995 0.000000        1                NaN       NaN         NaN  0.3042358


# compare dmr-containing scaffolds to the rest, required scaffold length >1mb
par(oma=c(0,0,2,0));
.verboseBoxplotColumns(slGeneStats_in70[slGeneStats_in70$length > 1e6,
                                        -grep('dmr|smoothed',names(slGeneStats_in70),ignore.case=T)], 
                       (slGeneStats_in70$dmrs > 0)[slGeneStats_in70$length > 1e6],
                       frame.plot=F, col='grey', border='darkgrey',
                       names=c('>0 dmrs', '0 dmrs'))
title('all scaffolds at least 1mb long input to BSmooth with >=70 CpGs', outer=T);


# histograms of stats for dmr-containing scaffolds
slGeneStatsDmrs = slGeneStats[slGeneStats$dmrs > 0, ];
slGeneStatsDmrs = as.data.frame(cbind(slGeneStatsDmrs, t(sapply(as.list(rownames(slGeneStatsDmrs)), 
                                                                function(f) apply(subset(dmrs, chr == f)[,7:15], 2,
                                                                                  mean)))));
par(mfrow=c(2,11), oma=c(0,0,2,0)); breaks=50;
for (i in 1:ncol(slGeneStatsDmrs)) {
  hist(slGeneStatsDmrs[,i], breaks=breaks, 
       col='grey', border='darkgrey', main=names(slGeneStatsDmrs)[i], xlab='');
  abline(v=quantile(slGeneStatsDmrs[,i], c(.05,.5,.95), na.rm=T), col='red')
}; rm(i);
title('all scaffolds with dmrs \nred lines = 5%, 50%, 95%', outer=T);

# histograms of stats for dmr-containing scaffolds with dmrsPerGene <= .2
goodrows = which(slGeneStats_in70$dmrsPerGene <= .2);
#
par(mfrow=c(2,7), oma=c(0,0,2,0)); breaks=50;
for (i in 1:ncol(slGeneStats_in70)) {
  hist(slGeneStats_in70[goodrows,i], breaks=breaks, col='grey', border='darkgrey', main=names(slGeneStats_in70)[i], xlab='');
  abline(v=quantile(slGeneStats_in70[goodrows,i], c(.05,.5,.95), na.rm=T), col='red')
}; rm(i);
title('scaffolds with >=70 CpGs and dmrsPerGene < .2\nred lines = 5%, 50%, 95%', outer=T);


#
par(oma=c(0,0,2,0))
.verboseScatterplotAllColumnPairs(slGeneStats[, c(8,10:13)], 
                                  col='grey', mfrow=c(3,5), 
                                  corOptions="method='p',use='p'");
title('all scaffolds, pearson', outer=T);

#
par(oma=c(0,0,2,0))
.verboseScatterplotAllColumnPairs(slGeneStats_in[, 8:13], 
                                  col='grey', mfrow=c(3,5), 
                                  corOptions="method='p',use='p'");
title('all scaffolds input to BSmooth, pearson', outer=T);

#
par(oma=c(0,0,2,0))
.verboseScatterplotAllColumnPairs(slGeneStats_in70[, 8:13], 
                                  col='grey', mfrow=c(3,5), 
                                  corOptions="method='p',use='p'");
title('all scaffolds input to BSmooth with >= 70 CpGs, pearson', outer=T);

#
par(oma=c(0,0,2,0))
.verboseScatterplotAllColumnPairs(slGeneStats[, 8:13], 
                                  col='grey', mfrow=c(3,5), 
                                  corOptions="method='s',use='p'");
title('all scaffolds, spearman', outer=T);

#
par(oma=c(0,0,2,0))
.verboseScatterplotAllColumnPairs(slGeneStats[goodrows, 8:13], 
                                  col='grey', mfrow=c(3,5), 
                                  corOptions="method='p',use='p'");
title('scaffolds with dmrsPerGene < .175, pearson', outer=T);

#
par(oma=c(0,0,2,0))
.verboseScatterplotAllColumnPairs(slGeneStats_in[goodrows, 8:13], 
                                  col='grey', mfrow=c(3,5), 
                                  corOptions="method='p',use='p'");
title('scaffolds input to BSMooth with dmrsPerGene < .175, pearson', outer=T);

par(oma=c(0,0,2,0))
.verboseScatterplotAllColumnPairs(slGeneStats_in70[goodrows, 8:13], 
                                  col='grey', mfrow=c(3,5), 
                                  corOptions="method='p',use='p'");
title('scaffolds input to BSMooth with >= 70 CpGs and dmrsPerGene < .2, pearson', outer=T);

#
par(oma=c(0,0,2,0))
.verboseScatterplotAllColumnPairs(slGeneStats[goodrows, 8:13], 
                                  col='grey', mfrow=c(3,5), 
                                  corOptions="method='s',use='p'");
title('scaffolds with dmrsPerGene < .175, spearman', outer=T);

#
par(oma=c(0,0,2,0));
.verboseBoxplotColumns(slGeneStats[-11], slGeneStats[11] < .175, 
                       col='grey', border='darkgrey', frame.plot=F); 
title(' scaffolds with < .175 dmrsPerGene',outer=T)

#
par(oma=c(0,0,2,0));
.verboseBoxplotColumns(slGeneStats_in[-11], slGeneStats_in[11] < .175, 
                       col='grey', border='darkgrey', frame.plot=F); 
title(' scaffolds input to BSmooth with < .175 dmrsPerGene',outer=T)



#
sgstoplot = slGeneStats;
par(mfrow=c(6,7), oma=c(0,0,2,0))
for (i in 8:ncol(sgstoplot)) {
  for (j in 1:7) {
    verboseScatterplot(sgstoplot[,i], sgstoplot[,j], xlab=names(sgstoplot)[i], ylab=names(sgstoplot)[j],
                       frame.plot=F, abline=T, abline.col='red',
                       pch=21, bg=numbers2colors(sgstoplot$dmrsPerGene), col='grey', cex=1.2)
  }
}
rm(i,j)
title('all scaffolds, dots colored by dmrsPerGene', outer=T)

#
sgstoplot = slGeneStats_in;
par(mfrow=c(6,7), oma=c(0,0,2,0))
for (i in 8:ncol(sgstoplot)) {
  for (j in 1:7) {
    verboseScatterplot(sgstoplot[,i], sgstoplot[,j], xlab=names(sgstoplot)[i], ylab=names(sgstoplot)[j],
                       frame.plot=F, abline=T, abline.col='red',
                       pch=21, bg=numbers2colors(sgstoplot$dmrsPerGene), col='grey', cex=1.2)
  }
}
rm(i,j)
title('all scaffolds input to BSmooth, dots colored by dmrsPerGene', outer=T)

#
sgstoplot = slGeneStats_in70;
par(mfrow=c(6,7), oma=c(0,0,2,0))
for (i in 8:ncol(sgstoplot)) {
  for (j in 1:7) {
    verboseScatterplot(sgstoplot[,i], sgstoplot[,j], xlab=names(sgstoplot)[i], ylab=names(sgstoplot)[j],
                       frame.plot=F, abline=T, abline.col='red',
                       pch=21, bg=numbers2colors(sgstoplot$dmrsPerGene), col='grey', cex=1.2)
  }
}
rm(i,j)
title('all scaffolds input to BSmooth with >=70 CpGs, dots colored by dmrsPerGene', outer=T)

#
sgstoplot = slGeneStats[goodrows, ];
par(mfrow=c(6,7), oma=c(0,0,2,0))
for (i in 8:ncol(sgstoplot)) {
  for (j in 1:7) {
    verboseScatterplot(sgstoplot[,i], sgstoplot[,j], xlab=names(sgstoplot)[i], ylab=names(sgstoplot)[j],
                       frame.plot=F, abline=T, abline.col='red',
                       pch=21, bg=numbers2colors(sgstoplot$dmrsPerGene), col='grey', cex=1.2)
  }
}
rm(i,j)
title('scaffolds with dmrsPerGene < .175, dots colored by dmrsPerGene', outer=T)


goodrows = which(slGeneStats_in$dmrsPerGene < .175);

#
sgstoplot = slGeneStats_in[goodrows, ];
par(mfrow=c(6,7), oma=c(0,0,2,0))
for (i in 8:ncol(sgstoplot)) {
  for (j in 1:7) {
    verboseScatterplot(sgstoplot[,i], sgstoplot[,j], xlab=names(sgstoplot)[i], ylab=names(sgstoplot)[j],
                       frame.plot=F, abline=T, abline.col='red',
                       pch=21, bg=numbers2colors(sgstoplot$dmrsPerGene), col='grey', cex=1.2)
  }
}
rm(i,j)
title('scaffolds input to BSmooth with dmrsPerGene < .175, dots colored by dmrsPerGene', outer=T)

#
sgstoplot = slGeneStats_in70[goodrows, ];
par(mfrow=c(7,7), oma=c(0,0,2,0))
for (i in 8:ncol(sgstoplot)) {
  for (j in 1:7) {
    verboseScatterplot(sgstoplot[,i], sgstoplot[,j], xlab=names(sgstoplot)[i], ylab=names(sgstoplot)[j],
                       frame.plot=F, abline=T, abline.col='red',
                       pch=21, bg=numbers2colors(sgstoplot$dmrsPerGene), col='grey', cex=1.2)
  }
}
rm(i,j)
title('scaffolds input to BSmooth with >=70 CpGs and dmrsPerGene < .2, dots colored by dmrsPerGene', outer=T)
