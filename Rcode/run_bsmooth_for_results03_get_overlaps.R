# clear workspace, set options and working directory, load functions
rm(list=ls()); options(stringsAsFactors=F); 
setwd('~/Documents/_BS-seq_analysis/');
library('GenomicRanges');
library('Rsamtools');
library('bsseq')
source('/Volumes/fishstudies-1/_code/run_bsmooth_functions.R');

###
# 
load('run_bsmooth_for_results02_load_annotations_and_prepWORKSPACE.RData');


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
rm(gr,featureGRnames);

# make version of overlap results that only contain gene names
dmrsOVsimple = lapply(dmrsOV[-which(names(dmrsOV) %in% c('te','z'))], 
                      function(f) split(f$mcols.geneSym, 
                                        f$queryHits));
dmrsOVsimple = lapply(dmrsOVsimple, function(f) lapply(f, unique));

# make a vector of all genes and brs overlapped by a dmr
alldmrgenes = unique(unlist(dmrsOVsimple));

# make version of overlap results where dmrs are organized by the genes they overlapped
dmrsOVbyGene = dmrsOV[names(dmrsOV) %in% c('br3prime','br','cds','exon','gene','intron')];
for (i in 1:length(dmrsOVbyGene)) {
  dmrsOVbyGene[[i]] = lapply(split(dmrsOVbyGene[[i]], 
                                   dmrsOVbyGene[[i]]$mcols.geneSym), 
                             function(f) dmrs[rownames(dmrs) %in% f$queryHits, ]);
}
rm(i);

# make df of all dmrs that hit genes or br
dmrsOVgenebrdf = do.call('rbind',lapply(dmrsOV[c('gene','br','br3prime')], function(f) f[,1:6]));

# make version with only dmrs and genes that map uniquely
x = table(dmrsOVgenebrdf$queryHits);
dmrsOVgenebrdf1 = subset(dmrsOVgenebrdf, queryHits %in% names(x)[x == 1]);
x = table(dmrsOVgenebrdf1$mcols.geneSym);
dmrsOVgenebrdf1 = subset(dmrsOVgenebrdf1, mcols.geneSym %in% names(x)[x == 1]);
rm(x);

# get names of dmrs that hit >1 gene and genes that were hit by >1 dmr
multiGeneDmrs = unique(dmrsOV$gene$queryHits[which(duplicated(dmrsOV$gene$queryHits))]);
multiDmrGenes = unique(dmrsOV$gene$mcols.geneSym[which(duplicated(dmrsOV$gene$mcols.geneSym))]);
multigenes = unique(c(multiDmrGenes, 
                      subset(dmrsOVgenebrdf, 
                             queryHits %in% multiGeneDmrs)$mcols.geneSym));

# make table of overlaps where rows are dmrs and columns are features
numOV = as.data.frame(matrix(0, nrow=nrow(dmrs), ncol=length(dmrsOV), 
                             dimnames=list(rownames(dmrs), names(dmrsOV))));
numOVrows = lapply(sapply(dmrsOV, 
                          function(f) f$queryHits), 
                   function(ff) match(unique(ff), rownames(numOV)));
for (i in 1:length(numOVrows)) {
  numOV[numOVrows[[i]], i] = 1;
}; rm(i,numOVrows);

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
}; rm(i,x,check,brows);

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
}; rm(i,x,check,b3rows);

# make a feature group vector
fg = as.character(apply(numOV[,-c(1,2,5)], 1, 
                        function(f) paste0(names(f)[f>0],collapse=':')));
names(fg) = rownames(numOV);

# test for enrichments of hypo/hyper dmrs for different features
for (i in 1:ncol(numOV)) {
  cat(names(numOV)[i],'...')
  print(.checkGeneListEnrichment(rownames(numOV)[numOV[,i]==1], 
                                 rownames(numOV)[dmrs$direction=='hypo'], 
                                 rownames(numOV)))
}; rm(i);

# -----------------------------------------
with_te = list(br=.getOVgenes(brGR,teGR),
               br3prime=.getOVgenes(br3primeGR,teGR),
               gene=.getOVgenes(geneGR,teGR),
               intron=unique(.getOVgenes(intronGR,teGR)),
               cds=unique(.getOVgenes(cdsGR,teGR)),
               exon=unique(.getOVgenes(exonGR,teGR)));

# geneGR_dmr_scaffolds = .fixSeqLevels(subset(geneGR, seqnames %in% scaffolds_with_dmrs));
# with_te_dmr_scaffolds = lapply(with_te, function(x) x[x %in% names(geneGR_dmr_scaffolds)]);

# gene_gene_OV = as.data.frame(findOverlaps(geneGR,drop.self=T,ignore.strand=T,drop.redundant=T));
# gene_gene_OVlist = apply(gene_gene_OV, 1, function(f) geneGR[unlist(f)]);
# gene_gene_OVlistGR = Reduce(c, gene_gene_OVlist);
# names(gene_gene_OVlistGR) = mcols(gene_gene_OVlistGR)[,1];
# gene_gene_OVlistGR_dmr_scaffolds = .fixSeqLevels(subset(gene_gene_OVlistGR, seqnames %in% scaffolds_with_dmrs));
# 
# gene_gene_OVlist_biotypes = lapply(gene_gene_OVlist, function(f) mcols(f)[,3]);
# rm(gene_gene_OV); gc();

lncrows = mcols(geneGR)[,3]=='lncRNA';
gene_lnc_OV = gene_lnc_OV=.findOverlapsNames(geneGR[!lncrows], 1, 
                                             geneGR[lncrows], 1, ignore.strand=T);
rm(lncrows)

gene_lnc_OVcodingGR = .fixSeqLevels(geneGR[names(geneGR) %in% unique(gene_lnc_OV$queryHits)]);
# gene_lnc_OVcodingGR_dmr_scaffolds = .fixSeqLevels(subset(gene_lnc_OVcodingGR, 
#                                                          seqnames %in% scaffolds_with_dmrs));

# gene_lnc_OVdmrsSyms = names(dmrOVbyGene)[names(dmrOVbyGene) %in% gene_lnc_OV$queryHits]

gffdesclookup = data.frame(sym=.parseGffMetaCol(subset(an$gff, V3=='mRNA'), 
                                                pattern='gene='), 
                           description=.parseGffMetaCol(subset(an$gff, V3=='mRNA'), 
                                                        pattern='product='));

# ----------------------------------------

save(list=ls(), file='run_bsmooth_for_results03_get_overlapsWORKSPACE.RData');


# collect objects into lists for better organization

# annotation and scaffold level info
annotations_and_scaffold_stats = list(abHsMap=abHsMap,
                                      an=an, 
                                      annoCombo=annoCombo, 
                                      burtoniGois=burtoniGois,
                                      gene_lnc_OV=gene_lnc_OV,
                                      gffdesclookup=gffdesclookup,
                                      slGeneStats=slGeneStats, 
                                      slGeneStats_in70=slGeneStats_in70,
                                      teOVwith_features=with_te);
save(annotations_and_scaffold_stats, 
     file='annotations_and_scaffold_stats_after03_get_overlaps.RData');
rm(abHsMap, an, annoCombo, burtoniGois, gene_lnc_OV, gffdesclookup, slGeneStats, slGeneStats_in70, with_te);

grs = list(br3primeGR=br3primeGR,
           brGR=brGR,
           cdsGR=cdsGR,
           dmrsGR=dmrsGR,
           exonGR=exonGR,
           geneGR=geneGR,
           intronGR=intronGR,
           teGR=teGR,
           zGR=zGR,
           gene_lnc_OVcodingGR=gene_lnc_OVcodingGR);
save(grs, 
     file='GRanges_objects_after03_get_overlaps.RData');
rm(list=grep('GR$',ls(),value=T));

dmrs_and_overlaps = list(dmrs=dmrs,
                         dmrsOV=dmrsOV,
                         dmrsOVbyGene=dmrsOVbyGene,
                         dmrsOVgenebrdf=dmrsOVgenebrdf,
                         dmrsOVgenebrdf1=dmrsOVgenebrdf1,
                         dmrsOVsimple=dmrsOVsimple,
                         numOV=numOV,
                         fg=fg);
save(dmrs_and_overlaps, 
     file='dmrs_and_overlaps_after03_get_overlaps.RData');
rm(list=grep('dmrs|OV|fg',ls(),value=T));

gene_and_dmr_lists = list(alldmrgenes=alldmrgenes,
                          multiDmrGenes=multiDmrGenes,
                          multiGeneDmrs=multiGeneDmrs,
                          multigenes=multigenes);
save(gene_and_dmr_lists, 
     file='gene_and_dmr_lists_after03_get_overlaps.RData');
rm(alldmrgenes, multiDmrGenes, multiGeneDmrs, multigenes);

#####
## plots

numOVtoplot = numOV[,-which(names(numOV) %in% c('gene','br','br3prime'))];
numOVtoplot = numOVtoplot[apply(numOVtoplot,1,sum) > 0, ];
numOVtoplot = as.data.frame(cbind(numOVtoplot, direction=rep(0,nrow(numOVtoplot))));
#numOVtoplot$direction[dmrs$direction=='hypo'] = 1;