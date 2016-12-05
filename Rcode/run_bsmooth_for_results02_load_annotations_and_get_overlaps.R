# clear workspace, set options and working directory, load functions
rm(list=ls()); options(stringsAsFactors=F); 
setwd('~/Documents/_BS-seq_analysis/');
#setwd('/Volumes/FISHSTUDIES/_methylation/');
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


# check dmrs hypo/hyper ratio by scaffold
hypo_hyper_byScaffold = sapply(split(dmrs, dmrs$chr), function(f) table(f$direction));
hypo_hyper_byScaffold_df = data.frame(hypo=rep(NA,length(hypo_hyper_byScaffold)), 
                                      hyper=rep(NA,length(hypo_hyper_byScaffold)));
rownames(hypo_hyper_byScaffold_df) = names(hypo_hyper_byScaffold);
for (i in 1:length(hypo_hyper_byScaffold)) {
  hypo_hyper_byScaffold_df[i,] = hypo_hyper_byScaffold[[i]][match(names(hypo_hyper_byScaffold_df), names(hypo_hyper_byScaffold[[i]]))]
}; rm(i);
hypo_hyper_byScaffold_df[is.na(hypo_hyper_byScaffold_df)] = 0;
hypo_hyper_byScaffold_df = as.data.frame(cbind(hypo_hyper_byScaffold_df,
                                               length=slGeneStats$length[match(rownames(hypo_hyper_byScaffold_df), 
                                                                               rownames(slGeneStats))]));
# # 
# x=c(); xnum=c(); 
# for (d in test_dists) { 
#   xnum=c(xnum, nrow(subset(dmrs, chr %in% rownames(subset(slGeneStats, length <= d))))) ; 
#   x=c(x, sum(subset(dmrs, chr %in% rownames(subset(slGeneStats, length <= d)))$direction == 'hypo') / nrow(subset(dmrs, chr %in% rownames(subset(slGeneStats, length <= d)))))    }
# verboseScatterplot(test_dists, x, abline=T, abline.col='red', type='p', pch=19, col='lightgrey', frame.plot=F); 
# text(test_dists, x+.004, labels=as.character(xnum))


# dmrs_1mb_scaffoldLengths = subset(dmrs, chr %in% rownames(subset(slGeneStats, length >= 1e6)));

# add dmr information to slGeneStats
dmrs_per_scaffold = table(dmrs$chr);
slGeneStats = as.data.frame(cbind(slGeneStats, 
                                  dmrs=dmrs_per_scaffold[match(rownames(slGeneStats),
                                                               names(dmrs_per_scaffold))]));
slGeneStats = as.data.frame(cbind(slGeneStats,
                                  smoothed=rownames(slGeneStats) %in% scaffolds_smoothed,
                                  avgDistBetweenDmrs=slGeneStats$length/slGeneStats$dmrs,
                                  dmrsPerMb=slGeneStats$dmrs/slGeneStats$length*1e6,
                                  dmrsPerGene=slGeneStats$dmrs/slGeneStats$numGenes));
slGeneStats$dmrs[is.na(slGeneStats$dmrs) & slGeneStats$smoothed] = 0;


load('aligned_trimmed4-98.adapters.q30.m0_bsmap2.9.bam_methratio_samtools0.1.19-m4-CpGcombined.CG_bsdAll.RData')
cpgcommon = .filterLowCoverage2(bsdAll, reqCov=1, fix.seqlevels=T);
slGeneStats_in = slGeneStats[rownames(slGeneStats) %in% unique(as.character(seqnames(cpgcommon))), ];
numlociPerScaffold = table(as.character(seqnames(cpgcommon)));
slGeneStats_in = as.data.frame(cbind(slGeneStats_in,
                                     numCpGs=numlociPerScaffold[match(rownames(slGeneStats_in), 
                                                                      names(numlociPerScaffold))]));
slGeneStats_in70 = subset(slGeneStats_in, numCpGs >= 70);

##########
#genes_pctCDS = sapply(geneGR_dmrScaffolds, function(f) sum(width(reduce(cdsGR[names(cdsGR)==names(f)]))) / width(f))
########

# get %CG for scaffolds
genomeFasta = FaFile('/Volumes/fishstudies-1/_Burtoni_genome_files/H_burtoni_v1.assembly.fa');
genomeSeqs = getSeq(genomeFasta);
genomeSeqsNucs = alphabetFrequency(genomeSeqs);
rownames(genomeSeqsNucs) = names(genomeSeqs);
genomeSeqsNucsCG = apply(genomeSeqsNucs, 1, 
                     function(f) sum(f[match(c('C','G'), colnames(genomeSeqsNucs))]) / sum(f) );


# get %CG for basal regulatory regions
brSeqs = getSeq(genomeFasta, brGR);
brSeqsNucs = alphabetFrequency(brSeqs);
rownames(brSeqsNucs) = mcols(brGR)[,1];
brSeqsNucsCG = apply(brSeqsNucs, 1, 
                     function(f) sum(f[match(c('C','G'), colnames(brSeqsNucs))]) / sum(f) );
mcols(brGR) <- as.data.frame(cbind(as.data.frame(mcols(brGR)), CG=brSeqsNucsCG));
rm(brSeqs, brSeqsNucs, brSeqsNucsCG);

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
mcols(dmrsGR) <- as.data.frame(cbind(name=names(dmrsGR), as.data.frame(mcols(dmrsGR)), CG=dmrsSeqsNucsCG));
rm(genomeFasta, dmrsSeqs, dmrsSeqsNucs, dmrsSeqsNucsCG, mcolCols);

# get overlaps between dmrs and genome features
#
# > grep('GR$',ls(),value=T)
# [1] "br3primeGR" "brGR" "cdsGR" "dmrsGR" "exonGR" "geneGR" "intronGR" "teGR" "zGR"

# all GRanges should have name of feature in first mcols column

featureGRnames = grep('GR$', ls(), value=T);
featureGRnames = featureGRnames[featureGRnames!='dmrsGR'];

dmrsOV = list();
for (gr in featureGRnames) {
  dmrsOV[[gsub('GR','',gr)]] = .findOverlapsWithSubjectMcols(dmrsGR, get(gr));
}
rm(gr);

genebrmaps = do.call('rbind',lapply(dmrsOV[c('gene','br','br3prime')], function(f) f[,1:6]));
genebrmaps1 = subset(genebrmaps, 
                     queryHits %in% names(table(genebrmaps$queryHits))[table(genebrmaps$queryHits) == 1])
genebrmaps1 = subset(genebrmaps1,
                     mcols.geneSym %in% names(table(genebrmaps1$mcols.geneSym))[table(genebrmaps1$mcols.geneSym) == 1])
# 
multiGeneDmrs = unique(dmrsOV$gene$queryHits[which(duplicated(dmrsOV$gene$queryHits))]);
multiDmrGenes = unique(dmrsOV$gene$mcols.geneSym[which(duplicated(dmrsOV$gene$mcols.geneSym))]);

# make version of overlap results where dmrs are organized by the genes they overlapped
dmrsOVbyGene = dmrsOV[names(dmrsOV) %in% c('br3prime','br','cds','exon','gene','intron')];
for (i in 1:length(dmrsOVbyGene)) {
  dmrsOVbyGene[[i]] = lapply(split(dmrsOVbyGene[[i]], 
                                   dmrsOVbyGene[[i]]$mcols.geneSym), 
                             function(f) dmrs[rownames(dmrs) %in% f$queryHits, ]);
}
rm(i);

#
dmrsOVsimple = lapply(dmrsOV[-which(names(dmrsOV) %in% c('te','z'))], function(f) split(f$mcols.geneSym, f$queryHits))
dmrsOVsimple = lapply(dmrsOVsimple, function(f) lapply(f, unique));
alldmrgenes = unique(unlist(dmrsOVsimple));
multigenes = unique(c(multiDmrGenes, subset(genebrmaps, queryHits %in% multiGeneDmrs)$mcols.geneSym))



# make table of overlaps where rows are dmrs and columns are features
numOV = as.data.frame(matrix(0, nrow=nrow(dmrs), ncol=length(dmrsOV), 
                             dimnames=list(rownames(dmrs), names(dmrsOV))));
numOVrows = lapply(sapply(dmrsOV, 
                          function(f) f$queryHits), 
                   function(ff) match(unique(ff), rownames(numOV)));
for (i in 1:length(numOVrows)) {
  numOV[numOVrows[[i]], i] = 1;
}; rm(i);

zs = rep(0,nrow(numOV))
numOV = as.data.frame(cbind(numOV, 
                            br.only=zs, br.same=zs, br.diff=zs,
                            br3.only=zs, br3.same=zs, br3.diff=zs));
rm(zs);

# break down br and br3prime hits
brows = which(numOV$br > 0);
for (i in brows) {
  x = lapply(dmrsOV, function(f) subset(f, queryHits==rownames(numOV)[i]));
  x = x[match(c('br','gene'), names(x))];
  if (nrow(x$gene) == 0) { numOV$br.only[i] = 1; next }
  if (nrow(x$br) == 1) {
    if (x$br$mcols.geneSym %in% x$gene$mcols.geneSym) {
      numOV$br.same[i] = 1;
    } else {
      numOV$br.diff[i] = 1;
    }
  } else {
    check = x$br$mcols.geneSym %in% x$gene$mcols.geneSym
    if (any(check)) { numOV$br.same[i] = 1 }
    if (any(!check)) { numOV$br.diff[i] = 1 }
  }
}; rm(i,x,check);

b3rows = which(numOV$br3prime > 0);
for (i in b3rows) {
  x = lapply(dmrsOV, function(f) subset(f, queryHits==rownames(numOV)[i]));
  x = x[match(c('br3prime','gene'), names(x))];
  if (nrow(x$gene) == 0) { numOV$br3.only[i] = 1; next }
  if (nrow(x$br3prime) == 1) {
    if (x$br3prime$mcols.geneSym %in% x$gene$mcols.geneSym) {
      numOV$br3.same[i] = 1;
    } else {
      numOV$br3.diff[i] = 1;
    }
  } else {
    check = x$br3prime$mcols.geneSym %in% x$gene$mcols.geneSym
    if (any(check)) { numOV$br3.same[i] = 1 }
    if (any(!check)) { numOV$br3.diff[i] = 1 }
  }
}; rm(i,x,check);

for (i in 1:ncol(numOV)) {
  cat(names(numOV)[i],'...')
  print(.checkGeneListEnrichment(rownames(numOV)[numOV[,i]==1], 
                                 rownames(numOV)[dmrs$direction=='hypo'], 
                                 rownames(numOV)))
}; rm(i);


fg = as.character(apply(numOV[,-c(1,2,5)], 1, 
                        function(f) paste0(names(f)[f>0],collapse=':')));
names(fg) = rownames(numOV)

numOVtoplot = numOV[,-which(names(numOV) %in% c('gene','br','br3prime'))];
numOVtoplot = numOVtoplot[apply(numOVtoplot,1,sum) > 0, ];
numOVtoplot = as.data.frame(cbind(numOVtoplot, direction=rep(0,nrow(numOVtoplot))));
#numOVtoplot$direction[dmrs$direction=='hypo'] = 1;







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

#
par(mfrow=c(2,7), oma=c(0,0,2,0)); breaks=50;
for (i in 1:ncol(slGeneStats)) {
  hist(slGeneStats[,i], breaks=breaks, col='grey', border='darkgrey', main=names(slGeneStats)[i], xlab='');
  abline(v=quantile(slGeneStats[,i], c(.05,.5,.95), na.rm=T), col='red')
}; rm(i);
title('all scaffolds\nred lines = 5%, 50%, 95%', outer=T);

par(mfrow=c(2,7), oma=c(0,0,2,0)); breaks=50;
for (i in 1:ncol(slGeneStats_in)) {
  hist(slGeneStats_in[,i], breaks=breaks, col='grey', border='darkgrey', main=names(slGeneStats_in)[i], xlab='');
  abline(v=quantile(slGeneStats_in[,i], c(.05,.5,.95), na.rm=T), col='red')
}; rm(i);
title('all scaffolds input to BSmooth\nred lines = 5%, 50%, 95%', outer=T);

par(mfrow=c(2,7), oma=c(0,0,2,0)); breaks=50;
for (i in 1:ncol(slGeneStats_in70)) {
  hist(slGeneStats_in70[,i], breaks=breaks, col='grey', border='darkgrey', main=names(slGeneStats_in70)[i], xlab='');
  abline(v=quantile(slGeneStats_in70[,i], c(.05,.5,.95), na.rm=T), col='red')
}; rm(i);
title('all scaffolds input to BSmooth with >=70 CpGs \nred lines = 5%, 50%, 95%', outer=T);


#
par(oma=c(0,0,2,0));
.verboseBoxplotColumns(slGeneStats[,1:7], is.na(slGeneStats$dmrs),
                       frame.plot=F, col='grey', border='darkgrey',
                       names=c('>0 dmrs', '0 dmrs'))
title('all scaffolds', outer=T);

# rbind(apply(slGeneStats[!is.na(slGeneStats$dmrs),1:7], 2, median, na.rm=T), 
#       apply(slGeneStats[is.na(slGeneStats$dmrs),1:7], 2, median, na.rm=T));
#        length numGenes avgGeneWidth avgGeneGC avgNumIsoforms pctBpGene genesPerMb
# [1,] 965616.0       31      19962.6   0.39685          1.901  0.020745   31.45046
# [2,]   4433.5        0       4421.0   0.38860          1.000  0.194200   70.12769

# rbind(apply(slGeneStats[!is.na(slGeneStats$dmrs),1:7], 2, mean, na.rm=T), 
#       apply(slGeneStats[is.na(slGeneStats$dmrs),1:7], 2, mean, na.rm=T));
#         length  numGenes avgGeneWidth avgGeneGC avgNumIsoforms  pctBpGene genesPerMb
# [1,] 1325566.2 43.202985     24077.20 0.3957894       1.922188 0.04689187   34.01965
# [2,]   50527.9  1.511479     10539.43 0.3833135       1.299743 0.30438105  152.94848

par(oma=c(0,0,2,0));
.verboseBoxplotColumns(slGeneStats_in[,1:7], is.na(slGeneStats_in$dmrs),
                       frame.plot=F, col='grey', border='darkgrey',
                       names=c('>0 dmrs', '0 dmrs'))
title('all scaffolds input to BSmooth', outer=T);

# rbind(apply(slGeneStats_in[!is.na(slGeneStats_in$dmrs),1:7], 2, median, na.rm=T), 
#       apply(slGeneStats_in[is.na(slGeneStats_in$dmrs),1:7], 2, median, na.rm=T));
#      length numGenes avgGeneWidth avgGeneGC avgNumIsoforms pctBpGene genesPerMb
# [1,] 965616       31    19962.600   0.39685          1.901  0.020745   31.45046
# [2,]   7044        0     5219.333   0.38830          1.000  0.162100   58.28968

# rbind(apply(slGeneStats_in[!is.na(slGeneStats_in$dmrs),1:7], 2, mean, na.rm=T), 
#       apply(slGeneStats_in[is.na(slGeneStats_in$dmrs),1:7], 2, mean, na.rm=T));
#          length  numGenes avgGeneWidth avgGeneGC avgNumIsoforms  pctBpGene genesPerMb
# [1,] 1325566.24 43.202985     24077.20 0.3957894       1.922188 0.04689187   34.01965
# [2,]   68222.03  2.022349     11530.95 0.3829383       1.327031 0.27338514  119.75714

par(oma=c(0,0,2,0));
.verboseBoxplotColumns(slGeneStats_in70[,1:7], is.na(slGeneStats_in70$dmrs),
                       frame.plot=F, col='grey', border='darkgrey',
                       names=c('>0 dmrs', '0 dmrs'))
title('all scaffolds input to BSmooth with >=70 CpGs', outer=T);

par(oma=c(0,0,2,0));
.verboseBoxplotColumns(slGeneStats_in70[slGeneStats_in70$length > 1e6,1:7], 
                       is.na(slGeneStats_in70$dmrs)[slGeneStats_in70$length > 1e6],
                       frame.plot=F, col='grey', border='darkgrey',
                       names=c('>0 dmrs', '0 dmrs'))
title('all scaffolds at least 1mb long input to BSmooth with >=70 CpGs', outer=T);


#
slGeneStatsDmrs = slGeneStats[!is.na(slGeneStats$dmrs), ];
slGeneStatsDmrs = as.data.frame(cbind(slGeneStatsDmrs, t(sapply(as.list(rownames(slGeneStatsDmrs)), 
                                                                function(f) apply(subset(dmrs, chr == f)[,7:15], 2,
                                                                                  mean)))))

par(mfrow=c(2,7), oma=c(0,0,2,0)); breaks=50;
for (i in 1:ncol(slGeneStatsDmrs)) {
  hist(slGeneStatsDmrs[,i], breaks=breaks, col='grey', border='darkgrey', main=names(slGeneStatsDmrs)[i], xlab='');
  abline(v=quantile(slGeneStatsDmrs[,i], c(.05,.5,.95), na.rm=T), col='red')
}; rm(i);
title('all scaffolds with dmrs \nred lines = 5%, 50%, 95%', outer=T);

#
goodrows = which(slGeneStats_in70$dmrsPerGene <= .2);
#
par(mfrow=c(2,7), oma=c(0,0,2,0)); breaks=50;
for (i in 1:ncol(slGeneStats)) {
  hist(slGeneStats_in70[goodrows,i], breaks=breaks, col='grey', border='darkgrey', main=names(slGeneStats_in70)[i], xlab='');
  abline(v=quantile(slGeneStats_in70[goodrows,i], c(.05,.5,.95), na.rm=T), col='red')
}; rm(i);
title('scaffolds with >=70 CpGs and dmrsPerGene < .2\nred lines = 5%, 50%, 95%', outer=T);


#
goodrows = which(slGeneStats$dmrsPerGene < .175);

#
par(mfrow=c(2,7), oma=c(0,0,2,0)); breaks=50;
for (i in 1:ncol(slGeneStats)) {
  hist(slGeneStats[goodrows,i], breaks=breaks, col='grey', border='darkgrey', main=names(slGeneStats)[i], xlab='');
  abline(v=quantile(slGeneStats[goodrows,i], c(.05,.5,.95), na.rm=T), col='red')
}; rm(i);
title('scaffolds with dmrsPerGene < .175\nred lines = 5%, 50%, 95%', outer=T);

#
goodrows = which(slGeneStats_in$dmrsPerGene < .175);

#
par(mfrow=c(2,7), oma=c(0,0,2,0)); breaks=50;
for (i in 1:ncol(slGeneStats_in)) {
  hist(slGeneStats_in[goodrows,i], breaks=breaks, col='grey', border='darkgrey', main=names(slGeneStats_in)[i], xlab='');
  abline(v=quantile(slGeneStats_in[goodrows,i], c(.05,.5,.95), na.rm=T), col='red')
}; rm(i);
title('scaffolds input to BSmooth with dmrsPerGene < .175\nred lines = 5%, 50%, 95%', outer=T);

#
par(oma=c(0,0,2,0))
.verboseScatterplotAllColumnPairs(slGeneStats[, 8:13], 
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

