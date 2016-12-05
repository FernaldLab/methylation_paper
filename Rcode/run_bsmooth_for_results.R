# clear workspace, set options and working directory, load functions
rm(list=ls()); options(stringsAsFactors=F); setwd('~/Documents/_BS-seq_analysis/');
library('GenomicRanges');
source('/Volumes/fishstudies-1/_code/run_bsmooth_functions.R');


# ==========================================================================================================
# ==========================================================================================================
# 
# ==========================================================================================================
# ==========================================================================================================
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# 
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------

# # load common dmrs
newprefix = 'CpGcombined.CG_4xCov_ns70_50_25hAll';
pattern = 'mg1e+08_fits_4xCov_Ttests_DMRs_cut0.95_mg300_tstat.corrected_n3_meanDiff0.1.RData';
# filebase = gsub('.RData', '', paste0(newprefix, pattern), fixed=T);
# load(paste0(filebase, 'exactAnd100bpOverlaps.RData'));
# 
# # combine overlapping dmrs and exact matches
# dmrs0 = as.data.frame(rbind(ovdmrs, mdmrs));
load('WORKSPACE_CpGcombined.CG_4xCov_ns70_50_25hAllmg1e+08_fits_4xCov_Ttests_DMRs_cut0.95_mg300_tstat.corrected_n3_meanDiff0.1.RData');
dmrs0 = DMRs25_50_70; rm(list=grep('70$',ls(),value=TRUE));
dmrs0ints = paste0(dmrs0$chr, ':', dmrs0$start, '-', dmrs0$end);

# load dataframe containing hand curation of dmrs
hc2 = read.table(gsub('RData','txt',paste0('curatedALL_',newprefix,pattern)), header=T, sep='\t');
hc2GR = GRanges(seqnames=hc2$chr, ranges=IRanges(start=hc2$start, end=hc2$end), mcols=hc2$status);
hc2ints = paste0(hc2$chr, ':', hc2$start, '-', hc2$end);
if (any(hc2ints[match(dmrs0ints,hc2ints)] != dmrs0ints)) { stop('hc2 and dmrs0 don\'t match') }

# get ids of dmrs that passed curation and filter
goodstatus = c('ok','ok_s','ok_d');
goodhc2ints = hc2ints[hc2$status %in% goodstatus];

dmrs = dmrs0[dmrs0ints %in% goodhc2ints, ];
dmrsints = paste0(dmrs$chr, ':', dmrs$start, '-', dmrs$end);
mcolCols = which(names(dmrs)=='n'):ncol(dmrs);
dmrsGR = GRanges(seqnames=dmrs$chr, ranges=IRanges(start=dmrs$start, end=dmrs$end), mcols=dmrs[,mcolCols]);

# rm(list=grep('dmrs',ls(),invert=T,value=T));
# save(list=ls(), file='run_bsmooth_for_results_')

# load null dmrs
load('WORKSPACE_run_bsmooth_generate_nulls.RData');

# find and filter out dmrs that are overlapped more than 50% by a null dmr
allnullints = unique(unlist(nullints));
allnullintsGR = .buildGR_for_intervals(allnullints)

nov = list()
for (d in 1:nrow(dmrs)) {
  thisw = dmrs$width[d];
  nov[[dmrsints[d]]] = as.data.frame(findOverlaps(dmrsGR[d], allnullintsGR, minoverlap=(thisw/2)))
}; rm(d,thisw);

dmrs = dmrs[!(dmrsints %in% names(nov)[sapply(nov, nrow)!=0]), ];

rm(list=ls()[ls()!='dmrs']);gc();
#save(dmrs,file='dmrs.RData')

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# load burtoni genome information and features with GRanges objects
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# load workspace from 'annotations_prep_for_dmrs_wgcna.R'
# should load:
#  handannos, sl0, slGenesStats, zcne, abHsMap, an, an2, annoCombo, sl

load('~/Documents/_annotationsDec2015/WORKSPACE_annotations_prep_for_dmrs_wgcna.RData')
rm(sl, sl0, handannos, an2, ); gc();

# get nucleotide counts for basal regulatory regions
genomeFasta = FaFile('/Volumes/fishstudies-1/_Burtoni_genome_files/H_burtoni_v1.assembly.fa');
brSeqs = getSeq(genomeFasta, brGR);
# could also use alphabetFrequency(brSeqs) to get brSeqsNucs
brSeqsNucs = lapply(brSeqs, function(f) table(unlist(strsplit(as.character(f),''))));
names(brSeqsNucs) = mcols(brGR)[,1]
brSeqsNucsCG = sapply(brSeqsNucs, function(f) (f[names(f)=='C'] + f[names(f)=='G'])/sum(f));
names(brSeqsNucsCG) = gsub('.C','',names(brSeqsNucsCG),fixed=T);

brSeqsNucsCGdmr = brSeqsNucsCG[match(names(dmrOVbyGene), names(brSeqsNucsCG))];


burtoniGois = as.vector(na.omit(annoCombo$new$gene[match(an$gois, annoCombo$new$Dbxref)]))

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# get/make GRanges objects for dmrs and find overlaps between genome features of interest
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# dmrs
dmrsints = paste0(dmrs$chr, ':', dmrs$start, '-', dmrs$end);
rownames(dmrs) = dmrsints;
mcolCols = which(names(dmrs)=='n'):ncol(dmrs);
dmrsGR = GRanges(seqnames=dmrs$chr, 
                 ranges=IRanges(start=dmrs$start, end=dmrs$end), 
                 mcols=dmrs[,mcolCols],
                 names=rownames(dmrs));
# dmrsGR20kb = .getSymmetricWindowsAroundFeatures(gr=dmrsGR, bpToAdd=20000, chrLengthVec=sl, sortOutput=F);
# mcols(dmrsGR20kb) <- as.data.frame(cbind(as.data.frame(mcols(dmrsGR20kb)), dmrsints));

# get DMR CG content
genomeFasta = FaFile('/Volumes/fishstudies-1/_Burtoni_genome_files/H_burtoni_v1.assembly.fa');
dmrSeqs = getSeq(genomeFasta, dmrsGR);
names(dmrSeqs) = rownames(dmrs);
dmrSeqsNucs = as.data.frame(alphabetFrequency(dmrSeqs));
dmrsCG = apply(dmrSeqsNucs, 1, function(f) sum(f[2:3]) / sum(f));
names(dmrsCG) = rownames(dmrs)

# gene_z_OV = .getOverlapsWithAmountsUnstranded(zGR, geneGR)
# gene_z_OVhits = gene_z_OV[!is.na(gene_z_OV)]
# gene_z_OVhitsGenes = lapply(gene_z_OVhits, function(f) f$mcols.geneSym)
# gene_z_OVhitsGenes_numz = table(unlist(gene_z_OVhitsGenes));

scaffolds_with_dmrs = unique(dmrs$chr);


# assumes gene name is in first column of mcols for gr1
.getOVgenes = function(gr1,gr2) {
  return(suppressWarnings(mcols(gr1)[unique(as.data.frame(findOverlaps(gr1,gr2))$queryHits), 1]));
}

with_te = list(br=.getOVgenes(brGR,teGR),
               br3prime=.getOVgenes(br3primeGR,teGR),
               gene=.getOVgenes(geneGR,teGR),
               intron=unique(.getOVgenes(intronGR,teGR)),
               cds=unique(.getOVgenes(cdsGR,teGR)),
               exon=unique(.getOVgenes(exonGR,teGR)));

geneGR_dmr_scaffolds = .fixSeqLevels(subset(geneGR, seqnames %in% scaffolds_with_dmrs));
with_te_dmr_scaffolds = lapply(with_te, function(x) x[x %in% names(geneGR_dmr_scaffolds)])


gene_gene_OV = as.data.frame(findOverlaps(geneGR,ignoreSelf=T,ignore.strand=T,ignoreRedundant=T));
gene_gene_OVlist = apply(gene_gene_OV, 1, function(f) geneGR[unlist(f)]);
gene_gene_OVlistGR = Reduce(c, gene_gene_OVlist);
names(gene_gene_OVlistGR) = mcols(gene_gene_OVlistGR)[,1];
gene_gene_OVlistGR_dmr_scaffolds = .fixSeqLevels(subset(gene_gene_OVlistGR, seqnames %in% scaffolds_with_dmrs));

gene_gene_OVlist_biotypes = lapply(gene_gene_OVlist, function(f) mcols(f)[,3]);
rm(gene_gene_OV); gc();

lncrows = mcols(geneGR)[,3]=='lncRNA';
gene_lnc_OV = as.data.frame(findOverlaps(geneGR[!lncrows], geneGR[lncrows], ignore.strand=T));
gene_lnc_OV$queryHits = mcols(geneGR[!lncrows])[gene_lnc_OV$queryHits,1];
gene_lnc_OV$subjectHits = mcols(geneGR[lncrows])[gene_lnc_OV$subjectHits,1];

gene_lnc_OVcodingGR = geneGR[names(geneGR) %in% unique(gene_lnc_OV$queryHits)];
gene_lnc_OVcodingGR_dmr_scaffolds = .fixSeqLevels(subset(gene_lnc_OVcodingGR, seqnames %in% scaffolds_with_dmrs));

gene_lnc_OVdmrsSyms = names(dmrOVbyGene)[names(dmrOVbyGene) %in% gene_lnc_OV$queryHits]

gffgenesyms = .parseGffMetaCol(subset(an$gff, V3=='mRNA'), pattern='gene=');
gffgenedesc = .parseGffMetaCol(subset(an$gff, V3=='mRNA'), pattern='product=');

#
geneGR_dmrs_distToNearest = as.data.frame(distanceToNearest(geneGR, dmrsGR))
geneGR_dmrs_distToNearest$queryHits = names(geneGR[geneGR_dmrs_distToNearest$queryHits]);
geneGR_dmrs_distToNearest$subjectHits = names(dmrsGR[geneGR_dmrs_distToNearest$subjectHits]);

geneGR_dmrs_distToNearest1mb = subset(geneGR_dmrs_distToNearest, distance <= 1e6);
geneGR_dmrs_distToNearest1mb_dmrcounts = table(geneGR_dmrs_distToNearest1mb$subjectHits);
geneGR_dmrs_distToNearest1mb_dmrcounts.95 = subset(geneGR_dmrs_distToNearest1mb, 
                                                   subjectHits %in% names(geneGR_dmrs_distToNearest1mb_dmrcounts[geneGR_dmrs_distToNearest1mb_dmrcounts > quantile(geneGR_dmrs_distToNearest1mb_dmrcounts, .95)]))
geneGR_dmrs_distToNearest1mb_dmrcounts.95byDmr = split(geneGR_dmrs_distToNearest1mb_dmrcounts.95, 
                                                       geneGR_dmrs_distToNearest1mb_dmrcounts.95$subjectHits)

geneGR_dmrs_distToNearest20kb = subset(geneGR_dmrs_distToNearest, distance <= 20000);
geneGR_dmrs_distToNearest20kb_dmrcounts = table(geneGR_dmrs_distToNearest20kb$subjectHits);
geneGR_dmrs_distToNearest20kb_dmrcounts.95 = subset(geneGR_dmrs_distToNearest20kb, 
                                                   subjectHits %in% names(geneGR_dmrs_distToNearest20kb_dmrcounts[geneGR_dmrs_distToNearest20kb_dmrcounts > quantile(geneGR_dmrs_distToNearest20kb_dmrcounts, .95)]))
geneGR_dmrs_distToNearest20kb_dmrcounts.95byDmr = split(geneGR_dmrs_distToNearest20kb_dmrcounts.95, 
                                                       geneGR_dmrs_distToNearest20kb_dmrcounts.95$subjectHits)

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# get dmr overlaps and nearest neighbors for all features
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# loop through different feature types to create lists of nearest neighbors and overlaps
features = c('gene', 'br', 'br3prime', 'mrna', 'exon', 'cds', 'intron', 'te', 'mi', 'z');
groups = list(hypo=dmrsints[dmrs$direction == 'hypo'], 
              hyper=dmrsints[dmrs$direction == 'hyper']);

# find overlaps and nearest neighbors for dmrs across features
.getDmrFeatureOverlapsAcrossGroups(dmrsGR, features, groups);
# .getDmrFeatureOverlapsAcrossGroupsNoNN(dmrsGR20kb, features, groups);

# get lists of genes for different features
gene = .process_dmrsGR_featureOV(dmrsGR_geneOV, dmrs, 6);
intron = .process_dmrsGR_featureOV(dmrsGR_intronOV, dmrs, 6);
cds = .process_dmrsGR_featureOV(dmrsGR_cdsOV, dmrs, 6);
te = .process_dmrsGR_featureOV(dmrsGR_teOV, dmrs, 6);
br = .process_dmrsGR_featureOV(dmrsGR_brOV, dmrs, 6);
br3prime = .process_dmrsGR_featureOV(dmrsGR_br3primeOV, dmrs, 6);

allOVgenes = unique(c(unlist(gene$byDmrSyms),
                      unlist(br$byDmrSyms),
                      unlist(br3prime$byDmrSyms)));

par(mfrow=c(2,2));
inds = rownames(res) %in% names(intron$byFeature)
verboseScatterplot(res$width[inds], res$isoforms[inds], abline=T);
inds = rownames(res) %in% names(cds$byFeature)
verboseScatterplot(res$width[inds], res$isoforms[inds], abline=T);
inds = rownames(res) %in% names(br$byFeature)
verboseScatterplot(res$width[inds], res$isoforms[inds], abline=T);
inds = rownames(res) %in% names(br3prime$byFeature)
verboseScatterplot(res$width[inds], res$isoforms[inds], abline=T);

# make lists that combine results from all features
dmrsGRhits = .combineOvAndNnByDmr(paste0('dmrsGR_',features));
dmrsGRhitsOV = lapply(dmrsGRhits, function(f) f$ov);
dmrsGRhitsNN = lapply(dmrsGRhits, function(f) f$nn);

# get num overlap hits for each feature type
numOV = .make_numOV(dmrsGRhitsOV, features);
#numOV$lncRNA = rep(0,nrow(numOV));
#numOV$lncRNA[rownames(numOV) %in% names(geneOVbyDmr[sapply(geneOVbyDmr, function(f) any(grepl('lncRNA',f)))])] = 1;

# add columns to delineate br hits for same gene that's overlapped, different gene, or no overlapping gene
numOV$br.only = rep(0,nrow(numOV));
numOV$br.same = rep(0,nrow(numOV));
numOV$br.diff = rep(0,nrow(numOV));
numOV$br3prime.only = rep(0,nrow(numOV));
numOV$br3prime.same = rep(0,nrow(numOV));
numOV$br3prime.diff = rep(0,nrow(numOV));

brDMRsCheckGene = .checkOVmatchGenes(dmrsGRhitsOV, numOV, 'br');
br3primeDMRsCheckGene = .checkOVmatchGenes(dmrsGRhitsOV, numOV, 'br3prime');

numOV$br.only[rownames(numOV) %in% brDMRsCheckGene$only] = 1;
numOV$br.same[rownames(numOV) %in% brDMRsCheckGene$same] = 1;
numOV$br.diff[rownames(numOV) %in% brDMRsCheckGene$diff] = 1;
numOV$br3prime.only[rownames(numOV) %in% br3primeDMRsCheckGene$only] = 1;
numOV$br3prime.same[rownames(numOV) %in% br3primeDMRsCheckGene$same] = 1;
numOV$br3prime.diff[rownames(numOV) %in% br3primeDMRsCheckGene$diff] = 1;



numOVyn = as.data.frame(apply(apply(numOV, 2, as.logical), 2, as.numeric));
numOVfeat=numOVyn[,-c(1:5)]
for(i in 1:ncol(numOVfeat)) { numOVfeat[(numOVfeat[,i] > 0) & !is.na(numOVfeat[,i]), i] = names(numOVfeat)[i] }; rm(i);
#numOVfeat$hypo[numOVfeat$hypo=='0'] = 'hyper';
fg = apply(numOVfeat, 1, function(f) paste0(f[f!='0' & !is.na(f)], collapse=':'))
fg[fg==''] = 'none';
names(fg) = rownames(numOV)
numOVynt = as.data.frame(t(numOVyn))
names(numOVynt) = fg;
tmp=t(numOVynt);


tmp=as.matrix(numOVyn[ apply(numOV[,-c(1:5,9,11)], 1, sum)   >0,-c(1:5,9,11)])
xxx=heatmap(tmp[,order(-apply(tmp, 2, sum))],col=heat.colors(2));
fgdendro = fg[apply(numOV[,-c(1:5,9,11)], 1, sum)   >0][xxx$rowInd];
fgdendro = fgdendro[length(fgdendro):1]

nnForNoOV = do.call('rbind',
        lapply(dmrsGRhitsNN[rownames(numOV[apply(numOV,1,sum)==0,])], 
               function(f) subset(do.call('rbind',f), 
                                  dist==min(dist,na.rm=T))))

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# load gene expression data
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
library(DESeq2);

# load estimated count data from kallisto
exprdir = '/Volumes/fishstudies-1/_LYNLEY_RNAseq/data/_kallisto_092915/';
nd1exp = read.table(paste0(exprdir, 'ATCACG/abundance.tsv'), header=T, sep='\t');
rownames(nd1exp) = sapply(strsplit(nd1exp$target_id, '|', fixed=T), function(f) f[length(f)]);
nd2exp = read.table(paste0(exprdir, 'TGACCA/abundance.tsv'), header=T, sep='\t');
rownames(nd2exp) = sapply(strsplit(nd2exp$target_id, '|', fixed=T), function(f) f[length(f)]);
d1exp = read.table(paste0(exprdir, 'CGATGT/abundance.tsv'), header=T, sep='\t');
rownames(d1exp) = sapply(strsplit(d1exp$target_id, '|', fixed=T), function(f) f[length(f)]);
d2exp = read.table(paste0(exprdir, 'TTAGGC/abundance.tsv'), header=T, sep='\t');
rownames(d2exp) = sapply(strsplit(d2exp$target_id, '|', fixed=T), function(f) f[length(f)]);

# make count matrix, columns are subjects, rows are transcripts
counts0 = as.data.frame(cbind('3157_TENNISON'=nd1exp$est_counts, 
                              '3165_BRISCOE'=d1exp$est_counts, 
                              '3581_LYNLEY'=d2exp$est_counts, 
                              '3677_MONK'=nd2exp$est_counts));

# sum across transcripts to get gene-level expression values
rownames(counts0) = rownames(nd1exp);
counts0 = counts0[match(an$lookup$transcript_id,rownames(counts0)), ];
counts0l = split(counts0,an$lookup$gene);
counts = as.data.frame(t(sapply(counts0l, function(f) apply(f, 2, sum))));
rm(counts0, counts0l); gc();

# make design matrix
coldata = as.data.frame(matrix(c('ND','D','D','ND'), ncol=1));
rownames(coldata) = names(counts);
names(coldata) = 'group';

dmrScaffoldsGenes = names(geneGR)[as.character(seqnames(geneGR)) %in% dmrScaffolds];
#counts_dmrScaffolds = counts[rownames(counts) %in% dmrScaffoldsGenes, ];

# run DESeq2
dds = DESeqDataSetFromMatrix(countData=round(counts), colData=coldata, design = ~ group);
dds = DESeq(dds, parallel=TRUE);
# dds_dmrScaffolds = dds[rownames(dds) %in% dmrScaffoldsGenes, ];
# res_raw_dmrScaffolds = results(dds_dmrScaffolds, contrast=c('group','D','ND'));
res_raw = results(dds, contrast=c('group','D','ND'));
res = as.data.frame(res_raw);
res = as.data.frame(cbind(res[,c(1,2,6)], 
                          an$gffGenesDF[match(rownames(res), an$gffGenesDF$geneSym), ],
                          brCG=brSeqsNucsCG[match(rownames(res), names(brSeqsNucsCG))]));


allGeneNearestDmr = as.data.frame(distanceToNearest(geneGR, dmrsGR));
allGeneNearestDmr$queryHits = mcols(geneGR)[allGeneNearestDmr$queryHits,1];
allGeneNearestDmr$subjectHits = names(dmrsGR)[allGeneNearestDmr$subjectHits];
allGeneNearestDmr1mb = subset(allGeneNearestDmr, distance <= 1e6);

allGeneNearestDmr1mb_dmrcounts = table(allGeneNearestDmr1mb$subjectHits);
allGeneNearestDmr1mb_dmrlist = split(allGeneNearestDmr1mb, allGeneNearestDmr1mb$subjectHits);
allGeneNearestDmr1mb_dmrdf = data.frame(numgenes=sapply(allGeneNearestDmr1mb_dmrlist,nrow),
                                        meandist=sapply(allGeneNearestDmr1mb_dmrlist, function(f) mean(f$distance)));

allGeneNearestDmr1mb_dmrdf = as.data.frame(cbind(allGeneNearestDmr1mb_dmrdf, 
                                                 n=dmrs$n[match(rownames(allGeneNearestDmr1mb_dmrdf), rownames(dmrs))],
                                                 width=dmrs$width[match(rownames(allGeneNearestDmr1mb_dmrdf), rownames(dmrs))],
                                                 invdensity=dmrs$invdensity[match(rownames(allGeneNearestDmr1mb_dmrdf), rownames(dmrs))],
                                                 meanDiff=dmrs$meanDiff[match(rownames(allGeneNearestDmr1mb_dmrdf), rownames(dmrs))],
                                                 areaStat=dmrs$areaStat[match(rownames(allGeneNearestDmr1mb_dmrdf), rownames(dmrs))],
                                                 direction=dmrs$direction[match(rownames(allGeneNearestDmr1mb_dmrdf), rownames(dmrs))],
                                                 fc=dmrsFc[match(rownames(allGeneNearestDmr1mb_dmrdf), names(dmrsFc))],
                                                 dmrrank=na.omit(match(rownames(allGeneNearestDmr1mb_dmrdf),rownames(dmrs)))));

allGeneNearestDmr1mb_dmrdf_mult = subset(allGeneNearestDmr1mb_dmrdf, numgenes > 1);

par(mfrow=c(2,2)); breaks=25; col='grey'; border='grey'
hist(allGeneNearestDmr1mb_dmrdf$numgenes, breaks=breaks, col=col, border=border, 
     main=paste0(.quantile.tails(allGeneNearestDmr1mb_dmrdf$numgenes),collapse=', '));
abline(v=.quantile.tails(allGeneNearestDmr1mb_dmrdf$numgenes), col='red')
hist(allGeneNearestDmr1mb_dmrdf$meandist, breaks=breaks, col=col, border=border, 
     main=paste0(round(.quantile.tails(allGeneNearestDmr1mb_dmrdf$meandist)),collapse=', '));
abline(v=.quantile.tails(allGeneNearestDmr1mb_dmrdf$meandist), col='red');

hist(allGeneNearestDmr1mb_dmrdf_mult$numgenes, breaks=breaks, col=col, border=border, 
     main=paste0(.quantile.tails(allGeneNearestDmr1mb_dmrdf_mult$numgenes),collapse=', '));
abline(v=.quantile.tails(allGeneNearestDmr1mb_dmrdf_mult$numgenes), col='red')
hist(allGeneNearestDmr1mb_dmrdf_mult$meandist, breaks=breaks, col=col, border=border, 
     main=paste0(round(.quantile.tails(allGeneNearestDmr1mb_dmrdf_mult$meandist)),collapse=', '));
abline(v=.quantile.tails(allGeneNearestDmr1mb_dmrdf_mult$meandist), col='red')




par(oma=c(0,0,2,0))
.verboseScatterplotAllColumnPairs(allGeneNearestDmr1mb_dmrdf[,-8],
                                  bg=labels2colors(allGeneNearestDmr1mb_dmrdf$direction),
                                  mfrow=c(4,9), pch=21, col='black');
title('all dmrs near >0 genes', outer=T);

par(oma=c(0,0,2,0))
.verboseScatterplotAllColumnPairs(abs(allGeneNearestDmr1mb_dmrdf[,-8]),
                                  bg=labels2colors(allGeneNearestDmr1mb_dmrdf$direction),
                                  mfrow=c(4,9), pch=21, col='black');
title('all dmrs near >0 genes, absolute values', outer=T);

par(oma=c(0,0,2,0))
.verboseScatterplotAllColumnPairs(allGeneNearestDmr1mb_dmrdf_mult[,-8],
                                  bg=labels2colors(allGeneNearestDmr1mb_dmrdf_mult$direction),
                                  mfrow=c(4,9), pch=21, col='black');
title('all dmrs near >1 genes', outer=T);

par(oma=c(0,0,2,0))
.verboseScatterplotAllColumnPairs(abs(allGeneNearestDmr1mb_dmrdf_mult[,-8]),
                                  bg=labels2colors(allGeneNearestDmr1mb_dmrdf_mult$direction),
                                  mfrow=c(4,9), pch=21, col='black');
title('all dmrs near >1 genes, absolute values', outer=T);

ngenerows = allGeneNearestDmr1mb_dmrdf_mult$numgenes > 15;

par(oma=c(0,0,2,0))
.verboseScatterplotAllColumnPairs(allGeneNearestDmr1mb_dmrdf_mult[ngenerows,-8],
                                  bg=labels2colors(allGeneNearestDmr1mb_dmrdf_mult$direction[ngenerows]),
                                  mfrow=c(4,9), pch=21, col='black');
title('all dmrs near >15 genes', outer=T);

par(oma=c(0,0,2,0))
.verboseScatterplotAllColumnPairs(allGeneNearestDmr1mb_dmrdf_mult[!ngenerows,-8],
                                  bg=labels2colors(allGeneNearestDmr1mb_dmrdf_mult$direction[!ngenerows]),
                                  mfrow=c(4,9), pch=21, col='black');
title('all dmrs near <=15 genes', outer=T);



par(oma=c(0,0,2,0))
.verboseScatterplotAllColumnPairs(abs(allGeneNearestDmr1mb_dmrdf_mult[ngenerows,-8]),
                                  bg=labels2colors(allGeneNearestDmr1mb_dmrdf_mult$direction[ngenerows]),
                                  mfrow=c(4,9), pch=21, col='black');
title('all dmrs near >15 genes, absolute values', outer=T);




# par(oma=c(0,0,2,0))
# .verboseScatterplotAllColumnPairs(subset(allGeneNearestDmr1mb_dmrdf, direction=='hypo')[,-8],
#                                   bg=labels2colors(sapply(strsplit(rownames(subset(allGeneNearestDmr1mb_dmrdf, 
#                                                                                    direction=='hypo')),':'), 
#                                                           function(x) x[1])),
#                                   mfrow=c(4,9), pch=21, col='black');
# title('all dmrs near >0 genes, D-hyper', outer=T);

dmrs_ranks = .orderDMRsByCompositeRank(dmrs, rankCols=c('n','width','areaStat','maxStat','meanDiff'))$ranks;
dmrs_ranks = dmrs_ranks[rownames(dmrs_ranks) %in% rownames(allGeneNearestDmr1mb_dmrdf), ]
dmrs_ranks = as.data.frame(cbind(dmrs_ranks, 
                                 allGeneNearestDmr1mb_dmrdf[na.omit(match(rownames(dmrs_ranks), 
                                                                          rownames(allGeneNearestDmr1mb_dmrdf))), 
                                                            c(1,2,8,9)]));

par(oma=c(0,0,2,0))
.verboseScatterplotAllColumnPairs(dmrs_ranks[,-9],
                                  bg=labels2colors(dmrs_ranks$direction),
                                  mfrow=c(4,9), pch=21, col='black');
title('all dmrs near >0 genes, ranks', outer=T);


allGeneNearestDmr1mb_dmrdf = as.data.frame(cbind(allGeneNearestDmr1mb_dmrdf,
                                                 'abs.areaStat'=abs(allGeneNearestDmr1mb_dmrdf$areaStat),
                                                 'abs.meanDiff'=abs(allGeneNearestDmr1mb_dmrdf$meanDiff),
                                                 'abs.fc'=abs(allGeneNearestDmr1mb_dmrdf$fc)));
allGeneNearestDmr1mb_dmrdf_mult = as.data.frame(cbind(allGeneNearestDmr1mb_dmrdf_mult,
                                                 'abs.areaStat'=abs(allGeneNearestDmr1mb_dmrdf_mult$areaStat),
                                                 'abs.meanDiff'=abs(allGeneNearestDmr1mb_dmrdf_mult$meanDiff),
                                                 'abs.fc'=abs(allGeneNearestDmr1mb_dmrdf_mult$fc)));
par(oma=c(0,0,2,0));
.verboseBoxplotColumns(allGeneNearestDmr1mb_dmrdf[-8], 
                       allGeneNearestDmr1mb_dmrdf$direction=='hypo',
                       names=c('ND-hyper','D-hyper'), frame.plot=F, col=c('lightblue','greenyellow'))
title('all dmrs near >0 genes', outer=T);

par(oma=c(0,0,2,0));
.verboseBoxplotColumns(allGeneNearestDmr1mb_dmrdf_mult[-8], 
                       allGeneNearestDmr1mb_dmrdf_mult$direction=='hypo',
                       names=c('ND-hyper','D-hyper'), frame.plot=F, col=c('lightblue','greenyellow'))
title('all dmrs near >1 genes', outer=T);

ngenerows = allGeneNearestDmr1mb_dmrdf_mult$numgenes > 15;
par(oma=c(0,0,2,0));
.verboseBoxplotColumns(allGeneNearestDmr1mb_dmrdf_mult[ngenerows,-8], 
                       allGeneNearestDmr1mb_dmrdf_mult$direction[ngenerows]=='hypo',
                       names=c('ND-hyper','D-hyper'), frame.plot=F, col=c('lightblue','greenyellow'))
title('all dmrs near >15 genes', outer=T);

par(oma=c(0,0,2,0));
.verboseBoxplotColumns(allGeneNearestDmr1mb_dmrdf_mult[!ngenerows,-8], 
                       allGeneNearestDmr1mb_dmrdf_mult$direction[!ngenerows]=='hypo',
                       names=c('ND-hyper','D-hyper'), frame.plot=F, col=c('lightblue','greenyellow'))
title('all dmrs near <=15 genes', outer=T);



par(mfrow=c(1,3));
verboseBoxplot(abs(allGeneNearestDmr1mb_dmrdf_mult$areaStat), 
               allGeneNearestDmr1mb_dmrdf_mult$direction=='hypo',ylab='',xlab='')
verboseBoxplot(abs(allGeneNearestDmr1mb_dmrdf_mult$areaStat)[ngenerows], 
               allGeneNearestDmr1mb_dmrdf_mult$direction[ngenerows]=='hypo',ylab='',xlab='');
verboseBoxplot(abs(allGeneNearestDmr1mb_dmrdf_mult$areaStat)[!ngenerows], 
               allGeneNearestDmr1mb_dmrdf_mult$direction[!ngenerows]=='hypo',ylab='',xlab='');



par(mfrow=c(2,8))
verboseBoxplot(abs(allGeneNearestDmr1mb_dmrdf_mult$areaStat)[ngenerows], 
               allGeneNearestDmr1mb_dmrdf_mult$direction[ngenerows]=='hypo',ylab='',xlab='');
verboseBoxplot(abs(allGeneNearestDmr1mb_dmrdf_mult$areaStat)[!ngenerows], 
               allGeneNearestDmr1mb_dmrdf_mult$direction[!ngenerows]=='hypo',ylab='',xlab='');
verboseBoxplot(abs(allGeneNearestDmr1mb_dmrdf_mult$meanDiff)[ngenerows], 
               allGeneNearestDmr1mb_dmrdf_mult$direction[ngenerows]=='hypo',ylab='',xlab='');
verboseBoxplot(abs(allGeneNearestDmr1mb_dmrdf_mult$meanDiff)[!ngenerows], 
               allGeneNearestDmr1mb_dmrdf_mult$direction[!ngenerows]=='hypo',ylab='',xlab='');
verboseBoxplot(abs(allGeneNearestDmr1mb_dmrdf_mult$fc)[ngenerows], 
               allGeneNearestDmr1mb_dmrdf_mult$direction[ngenerows]=='hypo',ylab='',xlab='');
verboseBoxplot(abs(allGeneNearestDmr1mb_dmrdf_mult$fc)[!ngenerows], 
               allGeneNearestDmr1mb_dmrdf_mult$direction[!ngenerows]=='hypo',ylab='',xlab='');
verboseBoxplot(abs(allGeneNearestDmr1mb_dmrdf_mult$numgenes)[ngenerows], 
               allGeneNearestDmr1mb_dmrdf_mult$direction[ngenerows]=='hypo',ylab='',xlab='');
verboseBoxplot(abs(allGeneNearestDmr1mb_dmrdf_mult$numgenes)[!ngenerows], 
               allGeneNearestDmr1mb_dmrdf_mult$direction[!ngenerows]=='hypo',ylab='',xlab='');
verboseBoxplot(abs(allGeneNearestDmr1mb_dmrdf_mult$meandist)[ngenerows], 
               allGeneNearestDmr1mb_dmrdf_mult$direction[ngenerows]=='hypo',ylab='',xlab='');
verboseBoxplot(abs(allGeneNearestDmr1mb_dmrdf_mult$meandist)[!ngenerows], 
               allGeneNearestDmr1mb_dmrdf_mult$direction[!ngenerows]=='hypo',ylab='',xlab='');

verboseBoxplot(abs(allGeneNearestDmr1mb_dmrdf_mult$n)[ngenerows], 
               allGeneNearestDmr1mb_dmrdf_mult$direction[ngenerows]=='hypo',ylab='',xlab='');
verboseBoxplot(abs(allGeneNearestDmr1mb_dmrdf_mult$n)[!ngenerows], 
               allGeneNearestDmr1mb_dmrdf_mult$direction[!ngenerows]=='hypo',ylab='',xlab='');

verboseBoxplot(abs(allGeneNearestDmr1mb_dmrdf_mult$width)[ngenerows], 
               allGeneNearestDmr1mb_dmrdf_mult$direction[ngenerows]=='hypo',ylab='',xlab='');
verboseBoxplot(abs(allGeneNearestDmr1mb_dmrdf_mult$width)[!ngenerows], 
               allGeneNearestDmr1mb_dmrdf_mult$direction[!ngenerows]=='hypo',ylab='',xlab='');

verboseBoxplot(abs(allGeneNearestDmr1mb_dmrdf_mult$invdensity)[ngenerows], 
               allGeneNearestDmr1mb_dmrdf_mult$direction[ngenerows]=='hypo',ylab='',xlab='');
verboseBoxplot(abs(allGeneNearestDmr1mb_dmrdf_mult$invdensity)[!ngenerows], 
               allGeneNearestDmr1mb_dmrdf_mult$direction[!ngenerows]=='hypo',ylab='',xlab='');



noDeGeneNearestDmrDist = as.data.frame(distanceToNearest(geneGR[!(mcols(geneGR)[,1] %in% rownames(res.1))], dmrsGR))



res.1 = subset(res[order(res$padj), ], padj<.1);
res.1GR = GRanges(seqnames=res.1$seqnames, 
                  ranges=IRanges(start=res.1$start, end=res.1$end), 
                  strand=res.1$strand,
                  mcols=res.1[,c(9,1:3,10:13)]);
res.1nearestDmrDist = as.data.frame(distanceToNearest(res.1GR, dmrsGR));


# check for DMRs overlapping differentially expressed genes
res.1OV = .searchForGenesAcrossNestedDmrHitList(dmrsGRhitsOV, rownames(res.1), 6);

# check for DMRs with 50kb, 100kb of differentially expressed genes

res.1NNup50kb = .searchForGenesAcrossNestedDmrHitList(dmrsGRhitsNNup50kb, rownames(res.1), 1);
res.1NNdown50kb = .searchForGenesAcrossNestedDmrHitList(dmrsGRhitsNNdown50kb, rownames(res.1), 1);

res.1NNup100kb = .searchForGenesAcrossNestedDmrHitList(dmrsGRhitsNNup100kb, rownames(res.1), 1);
res.1NNdown100kb = .searchForGenesAcrossNestedDmrHitList(dmrsGRhitsNNdown100kb, rownames(res.1), 1);

#
DEgenes50kbFromDMR = unique(c(names(res.1OV), names(res.1NNdown50kb), names(res.1NNup50kb)));
allGeneNearestDmr1mbDEgenes50kb = subset(allGeneNearestDmr1mb, queryHits %in% DEgenes50kbFromDMR);
allGeneNearestDmr1mbDEgenes50kbNearby = subset(allGeneNearestDmr1mb, subjectHits %in% allGeneNearestDmr1mbDEgenes50kb$subjectHits & distance <= 50000);

# compute fold-change for dmrs
dmrsFc = log2(dmrs$group2.mean / dmrs$group1.mean);
names(dmrsFc) = rownames(dmrs);

# get results for lncRNAs
reslnc = subset(res, biotype=='lncRNA')
# > nrow(reslnc)
# [1] 2058
# > sum(is.na(reslnc$padj))
# [1] 1490

reslncFc = reslnc[!is.na(reslnc$padj), ];
reslncExpr = subset(reslnc, baseMean > 0);

rescoding = subset(res, biotype=='protein_coding');

res_noNApadj = res[!is.na(res$padj), ];
resexpr = subset(res, baseMean>0);


resexprdmrinds = rownames(resexpr) %in% names(dmrOVbyGene);
# compare expression levels and fc in dmr vs non-dmr genes
par(mfrow=c(2,3))
verboseBoxplot(resexpr$baseMean, resexprdmrinds,  
               ylab='mean expr', xlab='', names=c('no dmr','dmr'), main='resexpr',
               frame.plot=F, ylim=c(0,5000))
verboseBoxplot(resexpr$log2FoldChange, resexprdmrinds,  
               ylab='log2fc', xlab='', names=c('no dmr','dmr'), main='resexpr',
               frame.plot=F, ylim=c(-.2,.2));
verboseBoxplot(abs(resexpr$log2FoldChange), resexprdmrinds,  
               ylab='abs(log2fc)', xlab='', names=c('no dmr','dmr'), main='resexpr',
               frame.plot=F, ylim=c(0,.3))


resexprDmrScaffoldsdmrinds = rownames(resexprDmrScaffolds) %in% names(dmrOVbyGene);
# compare expression levels and fc in dmr vs non-dmr genes
#par(mfrow=c(1,3))
verboseBoxplot(resexprDmrScaffolds$baseMean, resexprDmrScaffoldsdmrinds,  
               ylab='mean expr', xlab='', names=c('no dmr','dmr'), main='resexprDmrScaffolds',
               frame.plot=F, ylim=c(0,5000))
verboseBoxplot(resexprDmrScaffolds$log2FoldChange, resexprDmrScaffoldsdmrinds,  
               ylab='log2fc', xlab='', names=c('no dmr','dmr'), main='resexprDmrScaffolds',
               frame.plot=F, ylim=c(-.2,.2));
verboseBoxplot(abs(resexprDmrScaffolds$log2FoldChange), resexprDmrScaffoldsdmrinds,  
               ylab='abs(log2fc)', xlab='', names=c('no dmr','dmr'), main='resexprDmrScaffolds',
               frame.plot=F, ylim=c(0,.3))


# compare fc corrected for mean expression level in dmr vs non-dmr genes
par(mfrow=c(2,2));
verboseBoxplot(lm(resexpr$log2FoldChange ~ resexpr$baseMean)$residuals, resexprdmrinds,  
               ylab='log2fc corrected for baseMean', xlab='', names=c('no dmr','dmr'), main='resexpr',
               frame.plot=F, ylim=c(-.2,.2));
verboseBoxplot(lm(abs(resexpr$log2FoldChange) ~ resexpr$baseMean)$residuals, resexprdmrinds,  
               ylab='abs(log2fc) corrected for baseMean', xlab='', names=c('no dmr','dmr'), main='resexpr',
               frame.plot=F, ylim=c(-.1,.3))

verboseBoxplot(lm(resexprDmrScaffolds$log2FoldChange ~ resexprDmrScaffolds$baseMean)$residuals, resexprDmrScaffoldsdmrinds,  
               ylab='log2fc corrected for baseMean', xlab='', names=c('no dmr','dmr'), main='resexprDmrScaffolds',
               frame.plot=F, ylim=c(-.2,.2));
verboseBoxplot(lm(abs(resexprDmrScaffolds$log2FoldChange) ~ resexprDmrScaffolds$baseMean)$residuals, resexprDmrScaffoldsdmrinds,  
               ylab='abs(log2fc) corrected for baseMean', xlab='', names=c('no dmr','dmr'), main='resexprDmrScaffolds',
               frame.plot=F, ylim=c(-.1,.3))

# compare expression levels and fc in D vs ND dmr genes
resexprDmr = resexpr[resexprdmrinds,];
resexprdmrhypoinds = rownames(resexprDmr) %in% geneOVbyDmrSymsUni.hypo;
par(mfrow=c(1,3));
verboseBoxplot(resexprDmr$baseMean, resexprdmrhypoinds,  
               ylab='mean expr', xlab='', names=c('ND-hyper','D-hyper'), main='resexprDmr',
               frame.plot=F, ylim=c(0,10000))
verboseBoxplot(resexprDmr$log2FoldChange, resexprdmrhypoinds,  
               ylab='log2fc', xlab='', names=c('ND-hyper','D-hyper'), main='resexprDmr',
               frame.plot=F, ylim=c(-.2,.2))
verboseBoxplot(abs(resexprDmr$log2FoldChange), resexprdmrhypoinds,  
               ylab='abs(log2fc)', xlab='', names=c('ND-hyper','D-hyper'), main='resexprDmr',
               frame.plot=F, ylim=c(0,.2));

#
par(mfrow=c(3,2))
verboseScatterplot(resexpr$baseMean, resexpr$log2FoldChange, abline=T, abline.col='red',
                   frame.plot=F, xlab='mean expr', ylab='log2fc',
                   corOptions="method='s'", main='spearman, resexpr',
                   bg=numbers2colors(as.numeric(resexprdmrinds)), col='black', pch=21);
verboseScatterplot(resexpr$baseMean, abs(resexpr$log2FoldChange), abline=T, abline.col='red',
                   frame.plot=F, xlab='mean expr', ylab='abs(log2fc)',
                   corOptions="method='s'", main='spearman, resexpr',
                   bg=numbers2colors(as.numeric(resexprdmrinds)), col='black', pch=21);

verboseScatterplot(resexprDmr$baseMean, resexprDmr$log2FoldChange, abline=T, abline.col='red',
                   frame.plot=F, xlab='mean expr', ylab='log2fc',
                   corOptions="method='s'", main='dmr genes only, spearman, resexpr',
                   bg=numbers2colors(as.numeric(resexprdmrhypoinds)), col='black', pch=21);
verboseScatterplot(resexprDmr$baseMean, abs(resexprDmr$log2FoldChange), abline=T, abline.col='red',
                   frame.plot=F, xlab='mean expr', ylab='abs(log2fc)',
                   corOptions="method='s'", main='dmr genes only, spearman, resexpr',
                   bg=numbers2colors(as.numeric(resexprdmrhypoinds)), col='black', pch=21);

verboseScatterplot(resexprDmrScaffolds$baseMean, resexprDmrScaffolds$log2FoldChange, abline=T, abline.col='red',
                   frame.plot=F, xlab='mean expr', ylab='log2fc',
                   corOptions="method='s'", main='spearman, resexprDmrScaffolds',
                   bg=numbers2colors(as.numeric(resexprDmrScaffoldsdmrinds)), col='black', pch=21);
verboseScatterplot(resexprDmrScaffolds$baseMean, abs(resexprDmrScaffolds$log2FoldChange), abline=T, abline.col='red',
                   frame.plot=F, xlab='mean expr', ylab='abs(log2fc)',
                   corOptions="method='s'", main='spearman, resexprDmrScaffolds',
                   bg=numbers2colors(as.numeric(resexprDmrScaffoldsdmrinds)), col='black', pch=21);




# plot expression and dmr fc against each other
dmrFCbyGeneMeans = sapply(dmrFCbyGene, mean);
resexprDmr = as.data.frame(cbind(resexprDmr, 
                                 dmrFCmean=dmrFCbyGeneMeans[match(rownames(resexprDmr), 
                                                                  names(dmrFCbyGeneMeans))]))
resexprDmrHigh = subset(resexprDmr, baseMean >= median(resexprDmr$baseMean));
resexprDmrLow = subset(resexprDmr, baseMean < median(resexprDmr$baseMean));

resexprDmr_brCGhigh = subset(resexprDmr, brCG >= median(resexprDmr$brCG));
resexprDmr_brCGlow = subset(resexprDmr, brCG < median(resexprDmr$brCG));

resexprDmr1brCG = subset(resexprDmr, brCG <= quantile(resexprDmr$brCG)[2]);
resexprDmr2brCG = subset(resexprDmr, brCG > quantile(resexprDmr$brCG)[2] & brCG <= quantile(resexprDmr$brCG)[3]);
resexprDmr3brCG = subset(resexprDmr, brCG > quantile(resexprDmr$brCG)[3]);

par(mfrow=c(2,3));
for (tmpres in c('resexprDmr','resexprDmrHigh','resexprDmrLow')) {
  thisbg = numbers2colors(as.numeric(rownames(get(tmpres)) %in% geneOVbyDmrSyms.hypo & rownames(get(tmpres)) %in% geneOVbyDmrSyms.hyper),
                          colors=c('lightblue','yellow'))
  verboseScatterplot(get(tmpres)$log2FoldChange, 
                     dmrFCbyGeneMeans[match(rownames(get(tmpres)), names(dmrFCbyGeneMeans))],
                     xlab='exprFc', ylab='dmrFc', main=paste0(tmpres,'\n'),
                     frame.plot=F, abline=T, abline.col='red',
                     bg=thisbg, pch=21, col='black', cex=1.5);
  abline(h=0,v=0,col='grey');
}; rm(tmpres,thisbg);

#par(mfrow=c(1,3));
for (tmpres in c('resexprDmr1brCG','resexprDmr2brCG','resexprDmr3brCG')) {
  thisbg = numbers2colors(as.numeric(rownames(get(tmpres)) %in% geneOVbyDmrSyms.hypo & rownames(get(tmpres)) %in% geneOVbyDmrSyms.hyper),
                          colors=c('lightblue','yellow'))
  verboseScatterplot(get(tmpres)$log2FoldChange, 
                     dmrFCbyGeneMeans[match(rownames(get(tmpres)), names(dmrFCbyGeneMeans))],
                     xlab='exprFc', ylab='dmrFc', main=paste0(tmpres,'\n'),
                     frame.plot=F, abline=T, abline.col='red',
                     bg=thisbg, pch=21, col='black', cex=1.5);
  abline(h=0,v=0,col='grey');
}; rm(tmpres,thisbg)


#
resexprDmr_intron = resexprDmr[rownames(resexprDmr) %in% intron$byDmrSyms, ]

fcsametest = apply(resexprDmr[,names(resexprDmr) %in% c('log2FoldChange','dmrFCmean')], 1, 
                   function(f) all(f<0) | all(f>0));
resexprDmrAnti = resexprDmr[!fcsametest, ];
resexprDmrSame = resexprDmr[fcsametest, ];


# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# examine gene overlaps
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# make list of gene hits only
#  names are dmrs, elements are data frames where each row is a gene
geneOVbyDmr = dmrsGR_geneOV[!is.na(dmrsGR_geneOV)];
geneOVbyDmr.hypo = dmrsGR_geneOV.hypo[!is.na(dmrsGR_geneOV.hypo)];
geneOVbyDmr.hyper = dmrsGR_geneOV.hyper[!is.na(dmrsGR_geneOV.hyper)];
# > length(geneOVbyDmr)
# [1] 522
# > length(geneOVbyDmr.hypo)
# [1] 365
# > length(geneOVbyDmr.hyper)
# [1] 157
# > table(sapply(geneOVbyDmr, nrow))
#   1   2 
# 507  15 
# > table(sapply(geneOVbyDmr.hypo, nrow))
#   1   2 
# 354  11 
# > table(sapply(geneOVbyDmr.hyper, nrow))
#   1   2 
# 153   4 

# test likelihood of 709 dmrs hitting 522 genes
# very slow
tmpnull2 = list();
geneScaffolds = unique(as.character(seqnames(geneGR)));
dmrGeneScaffolds = geneScaffolds[geneScaffolds %in% dmrs$chr];
dmrwidths = round(dmrs$width);
for (i in 1:10000) {
  if(i %% 100 == 0) {cat(i,'...')}
  tmp = c();
  for (s in dmrGeneScaffolds) {
    # cat(s,'...');
    dmrsNullScaffold = .generateIntervalsFromScaffold(s, 
                                                      slGeneStats$length[rownames(slGeneStats)==s], 
                                                      dmrwidths[dmrs$chr==s]);
    tmpov = findOverlaps(dmrsNullScaffold, geneGR[seqnames(geneGR )==s]);
    tmp = c(tmp, length(unique(as.data.frame(tmpov)$queryHits)))
  }
  tmpnull2[[i]] = tmp;
}; rm(i,tmp,dmrsNullScaffold);
summary(sapply(tmpnull2, sum));
sum(sapply(tmpnull2, sum) >= 522);

permtestnumDmrGenes = tmpnull2

###### testing
dmrwidthlist = lapply(split(dmrs, dmrs$chr), 
                      function(f) list(chr=strsplit(rownames(f)[1],':')[[1]][1], 
                                       widths=round(f$width)));


nullGRlist = list();
nullOVlist = list();
for (i in 1:100) {
  cat(i,'...')
  dmrnulls = lapply(dmrwidthlist, 
                    function(f) .generateIntervalsFromScaffold(f$chr, 
                                                               slGeneStats$length[rownames(slGeneStats)==f$chr], 
                                                               f$widths));
  dmrnulls = Reduce(c, GRangesList(dmrnulls));
  dmrnullsOV = suppressWarnings(findOverlaps(dmrnulls, geneGR));
  nullGRlist[[i]] = geneGR[as.data.frame(dmrnullsOV)$subjectHits];
  nullOVlist[[i]] = dmrnullsOV;
}; rm(i,dmrnulls,dmrnullsOV)

nullDMRstats = as.data.frame(cbind(genes = sapply(nullGRlist, length),
                                   lncRNAs = sapply(nullGRlist, 
                                                    function(f) sum(as.data.frame(f)$mcols.biotype=='lncRNA')),
                                   multDmrGenes = sapply(nullGRlist, 
                                                         function(f) sum(table(as.data.frame(f)$mcols.geneSym) > 1)),
                                   multGeneDmrs = sapply(nullOVlist, 
                                                         function(f) sum(duplicated(as.data.frame(f)$queryHits))),
                                   dmrGenesOVlnc = sapply(nullGRlist, 
                                                          function(f) sum(as.data.frame(f)$mcols.geneSym %in% gene_lnc_OV$queryHits)),
                                   t(sapply(nullGRlist, 
                                            function(f) apply(resexpr[,c(1,2,7,12:14)], 2, 
                                                              function(ff) kruskal.test(ff, 
                                                                                        rownames(resexpr) %in% as.data.frame(f)$mcols.geneSym)$statistic)))))
names(nullDMRstats)[6:ncol(nullDMRstats)] = paste0(names(nullDMRstats)[6:ncol(nullDMRstats)],
                                                   '.kwstat');
nullDMRstats = as.data.frame(cbind(nullDMRstats,
                                   abs.log2FoldChange.kwstat = sapply(nullGRlist, 
                                                                      function(f) kruskal.test(abs(resexpr$log2FoldChange), 
                                                                                               rownames(resexpr) %in% as.data.frame(f)$mcols.geneSym)$statistic),
                                   GS_brCGdiff.kwstat = sapply(nullGRlist, 
                                                               function(f) kruskal.test(resexpr$GC-resexpr$brCG, 
                                                                                        rownames(resexpr) %in% as.data.frame(f)$mcols.geneSym)$statistic)))

save(nullGRlist, nullOVlist, nullDMRstats, file='null_GR_OV_DMRstats_Oct2016.RData')


dmrrealstats = c(apply(resexpr[,c(1,2,7,12:14)], 2, 
                       function(f) kruskal.test(f, 
                                                rownames(resexpr) %in% names(dmrOVbyGene))$statistic),
                 kruskal.test(abs(resexpr$log2FoldChange), rownames(resexpr) %in% names(dmrOVbyGene))$statistic,
                 kruskal.test(resexpr$GC-resexpr$brCG, rownames(resexpr) %in% names(dmrOVbyGene))$statistic);
names(dmrrealstats)[7:8] = c('abs.log2FoldChange','GS_brCGdiff');
names(dmrrealstats) = paste0(names(dmrrealstats), '.kwstat');
dmrrealstats = dmrrealstats[match(names(nullDMRstats)[6:ncol(nullDMRstats)], names(dmrrealstats))];

for (i in 1:length(dmrrealstats)) {
  print(names(dmrrealstats)[i])
  print(sum(nullDMRstats[,match(names(dmrrealstats)[i], names(nullDMRstats))]  >= dmrrealstats[i]) / nrow(nullDMRstats));
}; rm(i);

mdiffs = c();
thiscol = which(names(resexpr) == 'log2FoldChange')
for (n in 1:length(nullGRlist)) {
  thesegenes = rownames(resexpr) %in% as.data.frame(nullGRlist[[n]])$mcols.geneSym;
  mdiffs = c(mdiffs, median(resexpr[thesegenes, thiscol]) - median(resexpr[!thesegenes, thiscol]));
}
rm(n,thesegenes)
###### end testing

###### testing

dmrscaffrows = rownames(slGeneStats) %in% dmrscaffs;
dmrscaffs = rownames(slGeneStats)[dmrscaffrows];

dmrwidths = round(dmrs$width);
dmrScaffoldLengths = slGeneStats[dmrscaffrows, ]$length;
names(dmrScaffoldLengths) = rownames(slGeneStats)[rownames(slGeneStats) %in% dmrGeneScaffolds];

dmrscaffcounts = table(dmrs$chr);
dmrscaffcounts = dmrscaffcounts[match(dmrscaffs, names(dmrscaffcounts))];



shuffGRlist = list();
shuffOVlist = list();
shuffOVlistGenes = list();
for (i in 1:10) {
  cat(i,'...');
  dmrwidthsShuffled = sample(dmrwidths, size=length(dmrwidths), replace=F);
  
#   dmrscaffsShuffled = sample(dmrs$chr, size=length(dmrwidths), replace=F, 
#                              prob=(dmrscaffcounts[match(dmrs$chr, names(dmrscaffcounts))]));
  
  dmrscaffsShuffled = sample(names(dmrscaffcounts), size=length(dmrwidths), replace=T,
                             prob=dmrscaffcounts);
  
  dmrshuff = as.data.frame(cbind(dmrscaffsShuffled, dmrwidthsShuffled));
  dmrshuff = split(dmrshuff, dmrshuff[,1]);
  dmrshufflist = lapply(dmrshuff, function(f) list(chr=f[1,1], widths=as.numeric(f[,2])));
  dmrshuffGR = lapply(dmrshufflist, 
                      function(f) .generateIntervalsFromScaffold(f$chr, 
                                                                 slGeneStats$length[rownames(slGeneStats)==f$chr], 
                                                                 f$widths));
  shuffGRlist[[i]] = Reduce(c, GRangesList(dmrshuffGR));
  shuffOVlist[[i]] = suppressWarnings(findOverlaps(shuffGRlist[[i]], geneGR));
  shuffOVlistGenes[[i]] = geneGR[as.data.frame(shuffOVlist[[i]])$subjectHits];
  
}; rm(i, dmrshuffGR);

shuffDMRstats = as.data.frame(cbind(genes = sapply(shuffOVlistGenes, length),
                                   lncRNAs = sapply(shuffOVlistGenes, 
                                                    function(f) sum(as.data.frame(f)$mcols.biotype=='lncRNA')),
                                   multDmrGenes = sapply(shuffOVlistGenes, 
                                                         function(f) sum(table(as.data.frame(f)$mcols.geneSym) > 1)),
                                   multGeneDmrs = sapply(shuffOVlist, 
                                                         function(f) sum(duplicated(as.data.frame(f)$queryHits))),
                                   dmrGenesOVlnc = sapply(shuffOVlistGenes, 
                                                          function(f) sum(as.data.frame(f)$mcols.geneSym %in% gene_lnc_OV$queryHits)),
                                   t(sapply(shuffOVlistGenes, 
                                            function(f) apply(resexprDmrScaffolds[,c(1,2,7,12:14)], 2, 
                                                              function(ff) kruskal.test(ff, 
                                                                                        rownames(resexprDmrScaffolds) %in% as.data.frame(f)$mcols.geneSym)$statistic)))))
names(shuffDMRstats)[6:ncol(shuffDMRstats)] = paste0(names(shuffDMRstats)[6:ncol(shuffDMRstats)],
                                                   '.kwstat');
shuffDMRstats = as.data.frame(cbind(shuffDMRstats,
                                   abs.log2FoldChange.kwstat = sapply(shuffOVlistGenes, 
                                                                      function(f) kruskal.test(abs(resexprDmrScaffolds$log2FoldChange), 
                                                                                               rownames(resexprDmrScaffolds) %in% as.data.frame(f)$mcols.geneSym)$statistic),
                                   GS_brCGdiff.kwstat = sapply(shuffOVlistGenes, 
                                                               function(f) kruskal.test(resexprDmrScaffolds$GC-resexprDmrScaffolds$brCG, 
                                                                                        rownames(resexprDmrScaffolds) %in% as.data.frame(f)$mcols.geneSym)$statistic)))
for (i in 1:length(dmrrealstats)) {
  print(names(dmrrealstats)[i])
  print(sum(shuffDMRstats[,match(names(dmrrealstats)[i], names(shuffDMRstats))]  >= dmrrealstats[i]) / nrow(shuffDMRstats));
}; rm(i);


########



nullintsGRlist = list();
nullintsOVlist = list();
nullintsOVlistGenes = list();
for (i in 1:1000) {
  cat(i,'...');
  nullintsSample = sample(nullints, size=length(dmrwidths), replace=F)
  nullintssplit = strsplit(nullintsSample,':');
  nullintssplit2 = strsplit(sapply(nullintssplit, function(f) f[2]),'-')
  nullintsGRlist[[i]] = GRanges(seqnames=sapply(nullintssplit, function(f) f[1]), 
                                ranges=IRanges(start=as.numeric(sapply(nullintssplit2, function(f) f[1])), 
                                               end=as.numeric(sapply(nullintssplit2, function(f) f[2]))), 
                             strand='*');
  nullintsOVlist[[i]] = suppressWarnings(findOverlaps(nullintsGRlist[[i]], geneGR));
  nullintsOVlistGenes[[i]] = geneGR[as.data.frame(nullintsOVlist[[i]])$subjectHits];
}


nullintsDMRstats = as.data.frame(cbind(genes = sapply(nullintsOVlistGenes, length),
                                    lncRNAs = sapply(nullintsOVlistGenes, 
                                                     function(f) sum(as.data.frame(f)$mcols.biotype=='lncRNA')),
                                    multDmrGenes = sapply(nullintsOVlistGenes, 
                                                          function(f) sum(table(as.data.frame(f)$mcols.geneSym) > 1)),
                                    multGeneDmrs = sapply(nullintsOVlist, 
                                                          function(f) sum(duplicated(as.data.frame(f)$queryHits))),
                                    dmrGenesOVlnc = sapply(nullintsOVlistGenes, 
                                                           function(f) sum(as.data.frame(f)$mcols.geneSym %in% gene_lnc_OV$queryHits)),
                                    t(sapply(nullintsOVlistGenes, 
                                             function(f) apply(resexprDmrScaffolds[,c(1,2,7,12:14)], 2, 
                                                               function(ff) kruskal.test(ff, 
                                                                                         rownames(resexprDmrScaffolds) %in% as.data.frame(f)$mcols.geneSym)$statistic)))))
names(nullintsDMRstats)[6:ncol(nullintsDMRstats)] = paste0(names(nullintsDMRstats)[6:ncol(nullintsDMRstats)],
                                                     '.kwstat');

########


# make list where names are dmrs, elements are vectors of gene symbols
geneOVbyDmrSyms = lapply(geneOVbyDmr, function(f) f$mcols.geneSym);
geneOVbyDmrSyms.hypo = lapply(geneOVbyDmr.hypo, function(f) f$mcols.geneSym);
geneOVbyDmrSyms.hyper = lapply(geneOVbyDmr.hyper, function(f) f$mcols.geneSym);

# make unnamed vector of genes overlapped by dmrs
geneOVbyDmrSymsUni = unique(unlist(geneOVbyDmrSyms));
geneOVbyDmrSymsUni.hypo = unique(unlist(geneOVbyDmrSyms.hypo));
geneOVbyDmrSymsUni.hyper = unique(unlist(geneOVbyDmrSyms.hyper));
# > length(geneOVbyDmrSymsUni)
# [1] 484
# > length(geneOVbyDmrSymsUni.hypo)
# [1] 350
# > length(geneOVbyDmrSymsUni.hyper)
# [1] 151

# make list where names are genes and elements are data frames where each rowname is a dmr
#  careful of rownames appended by '.1', '.2', etc. 
#   those are dmrs that hit >1 gene
dmrOVbyGene = do.call('rbind', geneOVbyDmr);
dmrOVbyGene = split(dmrOVbyGene, dmrOVbyGene$mcols.geneSym);
# replace data frames with dmr info instead of gene info
for (g in 1:length(dmrOVbyGene)) {
  xdmrnames = strsplit(rownames(dmrOVbyGene[[g]]), '.', fixed=T);
  xdmrnames = sapply(xdmrnames, function(f) f[1]);
  dmrOVbyGene[[g]] = dmrs[rownames(dmrs) %in% xdmrnames, ];
}; rm(g,xdmrnames);
# > table(sapply(dmrOVbyGene, nrow))
#   1   2   3   4 
# 442  34   5   3 
# > table(sapply(dmrOVbyGene[names(dmrOVbyGene) %in% geneOVbyDmrSymsUni.hypo], nrow))
#   1   2   3   4 
# 317  26   4   3 
# > table(sapply(dmrOVbyGene[names(dmrOVbyGene) %in% geneOVbyDmrSymsUni.hyper], nrow))
#   1   2   3   4 
# 125  21   3   2

#
dmr_gffGenes = subset(an$gffGenesDF, geneSym %in% names(dmrOVbyGene));
WGCNA::verboseBoxplot(dmr_gffGenes$isoforms, 
                      dmr_gffGenes$geneSym %in% unlist(geneOVbyDmrSyms[sapply(geneOVbyDmrSyms, length) > 1])
                      );

# MAYBE WRONG
dmrFCbyGene = dmrOVbyGene;
for (i in 1:length(dmrFCbyGene)) {
  dmrFCbyGene[[i]] = dmrsFc[match(rownames(dmrFCbyGene[[i]]), names(dmrsFc))];
}; rm(i);

dmr_areaStatbyGene = sapply(dmrOVbyGene, function(f) mean(f$areaStat));


# any relationship between DMR width and overall gene expression level?
cor.test(sapply(dmrOVbyGene, function(f) mean(f$width)),
         resexpr$baseMean[match(names(dmrOVbyGene), rownames(resexpr))])

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# examine genes within 50kb, 100kb
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
dmrsGRhitsNNup50kb = .subsetNestedNNHitListOnStreamAndDistance(dmrsGRhitsNN,'up',50000);
dmrsGRhitsNNdown50kb = .subsetNestedNNHitListOnStreamAndDistance(dmrsGRhitsNN,'down',50000);

dmrsGRhitsNNup100kb = .subsetNestedNNHitListOnStreamAndDistance(dmrsGRhitsNN,'up',100000);
dmrsGRhitsNNdown100kb = .subsetNestedNNHitListOnStreamAndDistance(dmrsGRhitsNN,'down',100000);

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# study 2-gene dmrs
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------

# test likelihood of 15 DMRs hitting two genes
gene_gene_OVlist2strand = gene_gene_OVlist[sapply(gene_gene_OVlist, 
                                                  function(f) length(unique(strand(f)))) 
                                           == 2];
# make GRanges object of all intervals with overlapping genes
gene_gene_OVranges = do.call('rbind', lapply(gene_gene_OVlist2strand, 
                                    function(f) as.data.frame(disjoin(ranges(f))[2,]) ) );
gene_gene_OVscaffolds = sapply(gene_gene_OVlist2strand, function(f) unique(as.character(seqnames(f))));
gene_gene_OV_GR = GRanges(seqnames=gene_gene_OVscaffolds, 
                          ranges=IRanges(start=gene_gene_OVranges$start, end=gene_gene_OVranges$end));
# get all dmr widths
dmrwidths = round(dmrs$width);

# loop through scaffolds with dmrs, for each scaffold:
#  generate random intervals the same width as existing dmrs from the scaffold
#  overlap these intervals against intervals with overlapping genes on the scaffold
tmpnull = list();
for (i in 1:10000) {
  if(i %% 100 == 0) {cat(i,'...')}
  tmp = c();
  for (s in unique(gene_gene_OVscaffolds)[unique(gene_gene_OVscaffolds) %in% dmrs$chr]) {
   # cat(s,'...');
    tmpov = findOverlaps(.generateIntervalsFromScaffold(s, 
                                                        slGeneStats$length[rownames(slGeneStats)==s], 
                                                        dmrwidths[dmrs$chr==s]),
                         gene_gene_OV_GR[seqnames(gene_gene_OV_GR )==s])
    tmp = c(tmp, length(unique(as.data.frame(tmpov)$queryHits)))
  }
  tmpnull[[i]] = tmp;
}
rm(s,tmpov);
x=c(); bin=1;
for (i in seq(1,length(tmpnull),bin)) {
  x = c(x, sum(sapply(tmpnull[1:i], sum) >= 15) / i);
}; rm(i);
plot(seq(1,length(tmpnull),bin),x,
     type='l', frame.plot=F, 
     xlab='run',ylab='cumulative p-value'); 
abline(h=.05,col='red')

# make ov lists for 2-gene dmrs
geneOVbyDmr2gene = geneOVbyDmr[sapply(geneOVbyDmr,nrow) > 1];
geneOVbyDmr2gene_lnc = geneOVbyDmr2gene[sapply(geneOVbyDmr2gene, function(f) any(grepl('lncRNA',f)))];
geneOVbyDmr2gene_coding = geneOVbyDmr2gene[sapply(geneOVbyDmr2gene, function(f) !any(grepl('lncRNA',f)))];

geneOVbyDmr2geneSyms = unique(unlist(lapply(geneOVbyDmr2gene, function(f) f$mcols.geneSym)));

geneOVbyDmr2geneSyms_paired_with_lnc = sapply(geneOVbyDmr2gene_lnc, 
                                              function(f) f$mcols.geneSym[f$mcols.biotype=='protein_coding']);
geneOVbyDmr2geneSyms_paired_with_lnc = c(geneOVbyDmr2geneSyms_paired_with_lnc, 'b3gnt9');
geneOVbyDmr2gene_lncSyms = unique(unlist(lapply(geneOVbyDmr2gene_lnc, 
                                                function(f) f$mcols.geneSym[f$mcols.biotype=='lncRNA'])));
geneOVbyDmr2gene_lncSyms = c(geneOVbyDmr2gene_lncSyms, 'LOC102307825');


# get expression fc for 2-gene dmr genes
deseq_dmrOVbyGeneFc2gene = deseq_dmrOVbyGeneFc[names(deseq_dmrOVbyGeneFc) %in% geneOVbyDmr2geneSyms];
# get dmr fc for 2-gene dmrs
dmrsFc2gene = dmrsFc[names(dmrsFc) %in% names(geneOVbyDmr2gene)];
# make vector of dmr fc for 2-gene dmr genes
geneOVbyDmrSyms2gene = geneOVbyDmrSyms[match(names(dmrsFc2gene), names(geneOVbyDmrSyms))];
dmrsFc2geneByGene = c();
for (i in 1:length(geneOVbyDmrSyms2gene)) {
  tmpfc = dmrsFc2gene[match(names(geneOVbyDmrSyms2gene)[i], names(dmrsFc2gene))];
  dmrsFc2geneByGene = c(dmrsFc2geneByGene, rep(tmpfc,2));
}
rm(i,tmpfc);
names(dmrsFc2geneByGene) = unlist(geneOVbyDmrSyms2gene)
dmrsFc2geneByGene = dmrsFc2geneByGene[match(names(deseq_dmrOVbyGeneFc2gene), names(dmrsFc2geneByGene))];

# compare expression fc of genes paired with lncs to genes paired with protein coding genes
deseq_dmrOVbyGeneFc2gene_coding = deseq_dmrOVbyGeneFc2gene[!(names(deseq_dmrOVbyGeneFc2gene) %in% geneOVbyDmr2gene_lncSyms)];

fc2gene_coding = as.data.frame(cbind(expr=deseq_dmrOVbyGeneFc2gene_coding, 
                                     dmr=dmrsFc2geneByGene[match(names(deseq_dmrOVbyGeneFc2gene_coding), names(dmrsFc2geneByGene))]));
par(mfrow=c(1,3));
verboseBoxplot(deseq_dmrOVbyGeneFc2gene_coding, 
               names(deseq_dmrOVbyGeneFc2gene_coding) %in% geneOVbyDmr2geneSyms_paired_with_lnc,
               ylab='log2(Expression fold-change)', xlab='',
               names=c('paired with protein-coding','paired with lncRNA'),
               notch=F, frame.plot=F,
               col='grey', border='darkgrey');
abline(h=0,col='red')
x = fc2gene_coding[!(rownames(fc2gene_coding) %in% geneOVbyDmr2geneSyms_paired_with_lnc), ]
verboseScatterplot(x$dmr,x$expr,abline=T,abline.col='red',frame.plot=F,ylab='exprfc',xlab='dmrfc',main='spearman',type='n',corOptions="method='s'");
text(x$dmr,x$expr,labels=rownames(x))
abline(h=0,v=0,col='grey')
x = fc2gene_coding[rownames(fc2gene_coding) %in% geneOVbyDmr2geneSyms_paired_with_lnc, ]
verboseScatterplot(x$dmr,x$expr,abline=T,abline.col='red',frame.plot=F,ylab='exprfc',xlab='dmrfc',main='spearman',type='n',corOptions="method='s'");
text(x$dmr,x$expr,labels=rownames(x))
abline(h=0,v=0,col='grey')
rm(x);


# > deseq_dmrOVbyGeneFc2gene
#         aff3        alg10       b3gnt9        cenph        dimt1        ephb6        flrt1       gabra5       gabrb3 LOC102294049 LOC102296986 LOC102297721 LOC102298438 LOC102302092 LOC102307825 LOC102310329 
# -0.032523125 -0.017330730 -0.097675966 -0.010101864 -0.031251942  0.008330519 -0.048222694  0.065617119  0.061826690 -0.120582326 -0.065005583  0.128975401 -0.131381951  0.014146160 -0.076033104  0.026970962 
# LOC102310435 LOC102310532 LOC102311396 LOC102314293 LOC106632365 LOC106632423 LOC106633225 LOC106633436        lrrn1      macrod1         mob2       ralgds        samd7        tprkb 
# 0.037236080 -0.039580264 -0.003824029  0.048148538  0.115002323 -0.058550695 -0.005715247           NA -0.145378198  0.168807856 -0.005746078 -0.057685741 -0.084089767  0.113232746


par(mfrow=c(1,2));
# compute correlation of expr fc in lncs and the protein-coding gene they're paired with
x = as.data.frame(matrix(nrow=length(geneOVbyDmr2gene_lnc),ncol=2,
                         dimnames=list(names(geneOVbyDmr2gene_lnc),c('lnc','coding'))));
for (i in 1:nrow(x)) {
  this = geneOVbyDmr2gene_lnc[[i]];
  x$lnc[i] = deseq_dmrOVbyGeneFc2gene[match(this$mcols.geneSym[this$mcols.biotype=='lncRNA'], 
                                            names(deseq_dmrOVbyGeneFc2gene))];
  x$coding[i] = deseq_dmrOVbyGeneFc2gene[match(this$mcols.geneSym[this$mcols.biotype=='protein_coding'], 
                                            names(deseq_dmrOVbyGeneFc2gene))];
}
rm(i)
x[8, ] = c(-0.076033104, -0.097675966);
rownames(x)[8] = 'scaffold_437:161994-162246';
verboseScatterplot(x$lnc,x$coding,abline=T,abline.col='red',
                   frame.plot=F,
                   ylab='coding expr fc',xlab='lnc expr fc',pch=19);
text(x$lnc,x$coding,labels=rownames(x))
abline(h=0,v=0,col='grey');

# compute correlation of expr fc in protein-coding gene pairs
x = as.data.frame(matrix(nrow=length(geneOVbyDmr2gene_coding),ncol=2,
                         dimnames=list(names(geneOVbyDmr2gene_coding),c('coding1','coding2'))));
for (i in 1:nrow(x)) {
  this = geneOVbyDmr2gene_coding[[i]];
  x[i,] = deseq_dmrOVbyGeneFc2gene[match(this$mcols.geneSym, names(deseq_dmrOVbyGeneFc2gene))];
}; rm(i)

x = x[-which(rownames(x)=='scaffold_437:161994-162246'),]
verboseScatterplot(x$coding1,x$coding2,abline=T,abline.col='red',
                   frame.plot=F,
                   ylab='coding2',xlab='coding1',pch=19);
text(x$coding1,x$coding2,labels=rownames(x))
abline(h=0,v=0,col='grey');

# compare the above results to correlations between gene-pair expr fc across genome
lncinds = sapply(gene_gene_OVlist, function(f) paste0(sort(f$mcols.biotype),collapse=':')) == 'lncRNA:protein_coding'
gene_gene_OVlist_lnc = gene_gene_OVlist[lncinds];
codinginds = sapply(gene_gene_OVlist, function(f) paste0(sort(f$mcols.biotype),collapse=':')) == 'protein_coding:protein_coding'
gene_gene_OVlist_coding = gene_gene_OVlist[codinginds];


x = as.data.frame(matrix(nrow=length(gene_gene_OVlist_lnc),ncol=2,
                         dimnames=list(1:length(gene_gene_OVlist_lnc),c('lnc','coding'))));
for (i in 1:length(gene_gene_OVlist_lnc)) {
  this = as.data.frame(gene_gene_OVlist_lnc[[i]]);#print(i)
  x$lnc[i] = deseq_dmrOVbyGeneFc[match(this$mcols.geneSym[this$mcols.biotype=='lncRNA'], 
                                            names(deseq_dmrOVbyGeneFc))];
  x$coding[i] = deseq_dmrOVbyGeneFc[match(this$mcols.geneSym[this$mcols.biotype=='protein_coding'], 
                                               names(deseq_dmrOVbyGeneFc))];
}
rm(i)








# examine expression levels for genes overlapped by dmrs
deseq_dmrOVbyGene = res[match(names(dmrOVbyGene), rownames(res)), ];
deseq_dmrOVbyGeneFc = deseq_dmrOVbyGene$log2FoldChange;
names(deseq_dmrOVbyGeneFc) = rownames(deseq_dmrOVbyGene);

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# examine gene basal regulatory region overlaps
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# make list of gene br hits only
#  names are dmrs, elements are data frames where each row is a gene
brOVbyDmr = dmrsGR_brOV[!is.na(dmrsGR_brOV)];
brOVbyDmr.hypo = dmrsGR_brOV.hypo[!is.na(dmrsGR_brOV.hypo)];
brOVbyDmr.hyper = dmrsGR_brOV.hyper[!is.na(dmrsGR_brOV.hyper)];
# > length(brOVbyDmr)
# [1] 147
# > length(brOVbyDmr.hypo)
# [1] 99
# > length(brOVbyDmr.hyper)
# [1] 48

# make list where names are dmrs, elements are vectors of gene symbols
brOVbyDmrSyms = lapply(brOVbyDmr, function(f) f$mcols.geneSym);
brOVbyDmrSyms.hypo = lapply(brOVbyDmr.hypo, function(f) f$mcols.geneSym);
brOVbyDmrSyms.hyper = lapply(brOVbyDmr.hyper, function(f) f$mcols.geneSym);

# ==========================================================================================================
# ==========================================================================================================
# scan dmrs for TFBSs
# ==========================================================================================================
# ==========================================================================================================
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# 
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
#

load('~/Documents/JASPAR_hormone_receptor_IDs_wholeGenomeHits_dataFrame.RData')
tfgenome = tfhitsMultiHormoneByScaffoldDf; rm(tfhitsMultiHormoneByScaffoldDf)
tfgenomeGR = GRanges(seqnames=tfgenome$seqnames, 
                     ranges=IRanges(start=tfgenome$start, end=tfgenome$end), 
                     strand=tfgenome$strand, 
                     mcols=tfgenome[,c(6,7,9:12)]);
dmrsGR_tfgenomeOV = as.data.frame(findOverlaps(dmrsGR, tfgenomeGR));
dmrsGR_tfgenomeOV_tfDf = as.data.frame(tfgenomeGR[dmrsGR_tfgenomeOV$subjectHits]);
dmrsGR_tfgenomeOV_tfDf_intervals = paste0(dmrsGR_tfgenomeOV_tfDf$seqnames, ':',
                                          dmrsGR_tfgenomeOV_tfDf$start, '-',
                                          dmrsGR_tfgenomeOV_tfDf$end);
names(dmrsGR_tfgenomeOV_tfDf_intervals) = dmrsGR_tfgenomeOV_tfDf$mcols.TF;
dmrsGR_tfgenomeOV_tfbed = .write_bedfile_for_intervals(dmrsGR_tfgenomeOV_tfDf_intervals,
                                                       'JASPAR_hormone_receptor_IDs_wholeGenomeHits_dataFrame_dmrHits.bed',
                                                       TRUE);
dmrsGR_tfgenomeOV_dmrs = names(dmrsGR)[unique(dmrsGR_tfgenomeOV$queryHits)];

numOVyn_dmr_tf_OV = numOVyn[rownames(numOVyn) %in% dmrsGR_tfgenomeOV_dmrs & numOVyn$br>0, -c(1:5,9)]


tfgenomeGRzOV_GR = tfgenomeGR[unique(as.data.frame(findOverlaps(tfgenomeGR,zGR))$queryHits)]


genomeFasta = FaFile('/Volumes/fishstudies-1/_Burtoni_genome_files/H_burtoni_v1.assembly.fa');
dmrSeqs = getSeq(genomeFasta, dmrsGR);
names(dmrSeqs) = dmrsints;

# 
dmrSeqsSubset = dmrSeqs[names(dmrSeqs) %in% goodints];
ID = 'MA0007.2';
Sys.time()
tfhits = .searchSeqsForTFBS(ID, dmrSeqsSubset, min.score='90%');
Sys.time()

# took ~5hrs with 'ok' dmrs (n=752)
IDvec = names(JASPAR2014SitesSeqs);
Sys.time()
tfhitsDmrsMulti = .searchSeqsForMultipleTFBS(IDvec, dmrSeqs, min.score='90%');
Sys.time()
save(tfhitsDmrsMulti, file='tfhitsDmrsMulti_dmrs_allJASPAR.RData')
tfhitsDmrsMultiByDmr = .groupMultiTfbsHitsByDmrs(tfhitsDmrsMulti, verbose=100);
tfhitsDmrsMultiByDmr_filt = lapply(tfhitsDmrsMultiByDmr, function(f) f[nchar(f$siteSeqs) > 9, ]);
tfhitsDmrsMultiByDmr_filt = tfhitsDmrsMultiByDmr_filt[sapply(tfhitsDmrsMultiByDmr_filt, nrow) > 0];

#
tfhitsDmrsMultiByDmr_filt_summaries = lapply(tfhitsDmrsMultiByDmr_filt, .summarizeDmrTfHits);
# tfhitsDmrsMultiByDmr_filt_summaries = list();
# for (i in 1:length(tfhitsDmrsMultiByDmr_filt)) {
#   cat(i,'...');
#   tfhitsDmrsMultiByDmr_filt_summaries[[names(tfhitsDmrsMultiByDmr_filt)[i]]] = .summarizeDmrTfHits(tfhitsDmrsMultiByDmr_filt[[i]])
# }; rm(i);

tfhitsDmrsMultiByDmr_filt_abspos = lapply(tfhitsDmrsMultiByDmr_filt[sapply(tfhitsDmrsMultiByDmr_filt, nrow) > 0], 
                                          .computeTfHitAbsolutePositions);
tfhitsDmrsMultiByDmr_filt_abspos_intervals = unlist(sapply(tfhitsDmrsMultiByDmr_filt_abspos, 
                                                           function(x) paste0(strsplit(x$seqnames[1], 
                                                                                       '[:-]')[[1]][1], 
                                                                              ':', x$start, '-', x$end)));
names(tfhitsDmrsMultiByDmr_filt_abspos_intervals) = unlist(lapply(tfhitsDmrsMultiByDmr_filt_abspos, 
                                                                  function(x) paste0(x$TF,x$strand)))

tfhitsDmrsMultiByDmr_filt_abspos_bed = .write_bedfile_for_intervals(tfhitsDmrsMultiByDmr_filt_abspos_intervals,
                                                                    'JASPAR_all_ids_dmrHits.bed', 
                                                                    TRUE);

#
tfidlookup = .buildTfIdNameLookup(tfhitsDmrsMulti);

#
tf = .analyzeTfHitsByDmr(tfhitsDmrsMulti, cpgcommon);

.quantile.5 = function (vec) { return(quantile(vec, seq(.05,1,.05))) }
.order2 = function(df, columns, decreasing=T) {
  if (is.character(columns)) {column=match(columns, colnames(df))}
  return(df[order(df[,columns], decreasing=decreasing), ])
}

tmp = subset(tf$byDmr[[thisdmr]], strand=='-');
tmp = subset(tmp, nchar(siteSeqs) > 9)
tmpnamesplit = strsplit(thisdmr,'[:-]')[[1]];
tmptfints = paste0(rep(tmpnamesplit[1],nrow(tmp)), ':',
                   as.numeric(tmpnamesplit[2]) + tmp$start, '-',
                   as.numeric(tmpnamesplit[2]) + tmp$end);
names(tmptfints) = tmp$TF;
tmptfintsGR = GRanges(seqnames=sapply(strsplit(tmptfints, '[:-]'), function(f) f[1]), 
                      ranges=IRanges(start=as.numeric(sapply(strsplit(tmptfints, '[:-]'), function(f) f[2])), 
                                     end=as.numeric(sapply(strsplit(tmptfints, '[:-]'), function(f) f[3]))), 
                      strand=tmp$strand, mcols=tmp[,9:13])

ar_hits = do.call('rbind',lapply(tf$hitsWithCpGs, function(f) f[f$TF=='AR',]))
ar_hits_numOV = numOV[match(ar_hits$seqnames, rownames(numOV)), ];

.verboseScatterplotAllColumnPairs(tf$stats[,4:7],mfrow=c(2,3),col=labels2colors(tf$stats$class),pch=20,cex=3);
.verboseBoxplotColumns(tf$stats[4:7],tf$stats$class,c(4,1),col=names(table(labels2colors(tf$stats$class))));
.verboseScatterplotVecAgainstColumns(abs(dmrs[,7:15]), tf$maxCpgByHit, col=labels2colors(dmrs$direction),cex=2,pch=20);
.verboseBoxplotColumns((dmrs[,7:15]), as.factor(tf$maxCpgByHit), col='grey',frame.plot=F);

#
tfgeneOV = tf$hitsWithCpGs[names(geneOVbyDmr)];
tfgeneOV = lapply(tf$hitsWithCpGs[names(geneOVbyDmr)], function(f) f[nchar(f$siteSeqs) > 10, ]);
tfgeneOV = tfgeneOV[!sapply(tfgeneOV, is.null) & unlist(sapply(tfgeneOV, nrow)) > 0];
for (i in 1:length(tfgeneOV)) {#print(i)
  tfgeneOV[[i]] = subset(tfgeneOV[[i]], strand %in% as.character(geneOVbyDmr[[names(tfgeneOV)[i]]]$strand))
}; rm(i)

tfbrOV = lapply(tf$hitsWithCpGs[names(brOVbyDmr)], function(f) f[nchar(f$siteSeqs) > 10, ]);
tfbrOV = tfbrOV[!sapply(tfbrOV, is.null) & unlist(sapply(tfbrOV, nrow)) > 0];
for (i in 1:length(tfbrOV)) {#print(i)
  tfbrOV[[i]] = subset(tfbrOV[[i]], strand %in% as.character(brOVbyDmr[[names(tfbrOV)[i]]]$strand))
}; rm(i)

tfbrOV_AR = do.call('rbind',lapply(tfbrOV, function(f) f[f$TF=='AR',]));
tfgeneOV_AR = do.call('rbind',lapply(tfgeneOV, function(f) f[f$TF=='AR',]))

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# compute GO and KEGG enrichments
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------

bg = unique(unlist(strsplit(annoCombo$new$hsaHomologEntrez, ',')));
org.Hs.egGO2ALLEGSlist = as.list(org.Hs.egGO2ALLEGS);

geneOVbyDmrSymsUniEntrez = .mapAbGenesHsEntrez(geneOVbyDmrSymsUni, abHsMap, 'ab')
geneOVbyDmrSymsUniEntrezLimit = .limitVectorsInList(geneOVbyDmrSymsUniEntrez, num=2, smartUnlist=T)$newvec;


geneOVbyDmrSymsUniHighEntrez = .mapAbGenesHsEntrez(rownames(resexprDmrHigh), abHsMap, 'ab')
geneOVbyDmrSymsUniHighEntrezLimit = .limitVectorsInList(geneOVbyDmrSymsUniHighEntrez, num=2, smartUnlist=T)$newvec;

sameFcEntrez = .mapAbGenesHsEntrez(rownames(resexprDmrSame), abHsMap, 'ab');
sameFcEntrezLimit = .limitVectorsInList(sameFcEntrez, num=2, smartUnlist=T)$newvec;

antiFcEntrez = .mapAbGenesHsEntrez(rownames(resexprDmrAnti), abHsMap, 'ab');
antiFcEntrezLimit = .limitVectorsInList(antiFcEntrez, num=2, smartUnlist=T)$newvec;

# IDS = list(gene_br.ids=unique(c(dmrsGR_geneOV.ids$ids$newvec, dmrsGR_brOV.ids$ids$newvec)),
#            gene_br.hypo.ids=unique(c(dmrsGR_geneOV.hypo.ids$ids$newvec, dmrsGR_brOV.hypo.ids$ids$newvec)),
#            gene_br.hyper.ids=unique(c(dmrsGR_geneOV.hyper.ids$ids$newvec, dmrsGR_brOV.hyper.ids$ids$newvec)));
IDS = list(all=geneOVbyDmrSymsUniEntrezLimit,
           same=sameFcEntrezLimit,
           anti=antiFcEntrezLimit
           )
for (fdrmeth in c('BY','BH')) {
  cat(paste0('------------------------- ', fdrmeth,  ' -------------------------'));
  for (fdrth in c(.01, .05, .1)) {
    cat(paste0('------------------------- ', fdrth,  ' -------------------------'));
    for (i in 1:length(IDS)) {
      cat(paste0('------------------------- ', names(IDS)[i],  ' -------------------------'));
      thisname = paste0('go.',names(IDS)[i],'_',fdrmeth,substr(fdrth, 2, nchar(fdrth)));
      thisres = .GOFunctionAllOntologies(IDS[[i]], bg, fdrmethod=fdrmeth, fdrth=fdrth);
      if (nrow(thisres) > 0) {
        assign(thisname, .addGenesToGOFunctionResults(thisres, modulegenes=IDS[[i]]));
      } else {
        assign(thisname, thisres); 
      }
    }
  }
}; rm(IDS,fdrmeth,fdrth,i,thisname,thisres);

for (i in grep('go.anti',ls(),fixed=T,value=T) ) { print(i);print(get(i)[-8]) }; rm(i);


.getGenesByTerm = function (GOFunctionOutput, geneCol, geneDelim, termCol, originalAbGenes=NULL) {
  x = strsplit(GOFunctionOutput[,geneCol], geneDelim, fixed=T);
  names(x) = GOFunctionOutput[,termCol];
  if (is.vector(originalAbGenes)) {
    x = lapply(x, function(f) originalAbGenes[originalAbGenes %in% f]);
  }
  return(x)
}

.getTermsByGene = function (getGenesByTermOutput) {
  x = names(unlist(getGenesByTermOutput));
  x = strsplit(x, '.', fixed=T);
  out = list();
  for (i in x) {
    if (i[2] %in% names(out)) {
      out[[i[2]]] = c(out[[i[2]]], i[1]);
    } else {
      out[[i[2]]] = i[1];
    }
  }
  return(out)
}

.parseGOhitGenes = function (GOFunctionOutput, geneCol, geneDelim, termCol, originalAbGenes=NULL) {
  genesByTerm = .getGenesByTerm(GOFunctionOutput, geneCol, geneDelim, termCol, originalAbGenes=originalAbGenes);
  termsByGene = .getTermsByGene(genesByTerm);
  return(list(byTerm=genesByTerm, byGene=termsByGene));
}

.buildGeneTermTable = function (parseGOhitGenesOutput) {
  bygene = parseGOhitGenesOutput$byGene;
  byterm = parseGOhitGenesOutput$byTerm;
  df = as.data.frame(matrix(nrow=length(bygene), ncol=length(byterm), 
                            dimnames=list(names(bygene),names(byterm))));
  for (g in 1:length(bygene)) {
    df[g, ] = as.numeric(names(df) %in% bygene[[g]]);
  }
  return(df)
}

x  = .parseGOhitGenes(go._BY.05, 8, ',', 3,geneOVbyDmrSymsUniEntrezLimit);
xt = .buildGeneTermTable(x);


fgcol = labels2colors(as.character(grepl('cds',fg)));
names(fgcol)=names(fg)
par(mfrow=c(3,8));
for (i in 1:22) {
  termgenes = names(x$byTerm[[i]]);
  termgenedmrs = unlist(geneOVbyDmrSyms[sapply(geneOVbyDmrSyms, function(f) any(f %in% termgenes))]);
  termgenedmrs = termgenedmrs[match(termgenes, termgenedmrs)];
  termgenedmrsfg = fgcol[match(names(termgenedmrs), names(fgcol))]
  verboseScatterplot(resexprDmr[match(termgenes, rownames(resexprDmr)), ]$log2FoldChange,
                     resexprDmr[match(termgenes, rownames(resexprDmr)), ]$dmrFCmean,
                     abline=T, abline.col='red', frame.plot=F, 
                     pch=21, col='black',bg=termgenedmrsfg, cex=2,
                     xlab='exprFc',ylab='dmrFc', main=paste0(names(x$byTerm)[i],'\n'));
  abline(h=0,v=0,col='grey');
};rm(i,termgenes,termgenedmrs,termgenedmrsfg)

lapply(x$byTerm, function(f) cor.test(resexprDmr[rownames(resexprDmr) %in% names(f), ]$log2FoldChange, resexprDmr[rownames(resexprDmr) %in% names(f), ]$dmrFCmean))




#########

GOterm_numTFs = as.data.frame(rbind(read.table('~/Documents/BrowseGOTcoF-DB.txt',header=F,sep='\t',quote=''),
                                    read.table('~/Documents/BrowseGOTcoF-DB_MF.txt',header=F,sep='\t',quote='')));

allTFs = read.table('~/Documents/BrowseTFTcoF-DB_allTFs.txt',header=T, sep='\t');
subset(GOterm_numTFs, V2 %in% go._BY.05$goid);

x  = .parseGOhitGenes(go._BY.05, 8, ',', 3,geneOVbyDmrSymsUniEntrezLimit);
x$byTerm$`receptor activity`