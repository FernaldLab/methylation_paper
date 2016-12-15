# clear workspace, set options and working directory, load functions
rm(list=ls()); options(stringsAsFactors=F); 
library(org.Hs.eg.db);
library(GOFunction);
library(clusterProfiler);

#setwd('~/Documents/_BS-seq_analysis/');
#source('/Volumes/fishstudies-1/_code/run_bsmooth_functions.R');
setwd('/Volumes/FISHSTUDIES/_methylation_paper_Nov16/');
source('Rcode/run_bsmooth_functions.R');

# ---------------------------------
# load lists from end of run_bsmooth_for_results03_get_overlaps.R
#for (f in grep('after03',list.files(),value=T)) { load(f) }; rm(f);
for (f in grep('after03',list.files('RData'),value=T)) { load(paste0('RData/',f)) }; rm(f);

for (l in ls()) { cat(l,':\n',names(get(l)),'\n') }; rm(l);
# annotations_and_scaffold_stats :
#   abHsMap an annoCombo burtoniGois gene_lnc_OV gffdesclookup slGeneStats slGeneStats_in70 teOVwith_features 
# dmrs_and_overlaps :
#   dmrs dmrsOV dmrsOVbyGene dmrsOVgenebrdf dmrsOVgenebrdf1 dmrsOVsimple numOV fg 
# gene_and_dmr_lists :
#   alldmrgenes multiDmrGenes multiGeneDmrs multigenes 
# grs :
#   br3primeGR brGR cdsGR dmrsGR exonGR geneGR intronGR teGR zGR gene_lnc_OVcodingGR 

# # extract all list subcomponents into workspace
# dlists = ls();
# for (l in dlists) {
#   thisl = get(l); for (i in 1:length(thisl)) { assign(names(thisl)[i], thisl[[i]]) }
# }; rm(list=dlists); rm(l,thisl,i,dlists); gc();
# 
# ls()
# # [1]  "abHsMap"   "alldmrgenes"    "an"             "annoCombo"       "br3primeGR"       "brGR"          "burtoniGois"       "cdsGR"  "dmrs"            
# # [10] "dmrsGR"    "dmrsOV"         "dmrsOVbyGene"   "dmrsOVgenebrdf"  "dmrsOVgenebrdf1"  "dmrsOVsimple"  "exonGR"            "fg"     "geneGR"          
# # [19] "intronGR"  "multiDmrGenes"  "multiGeneDmrs"  "multigenes"      "numOV"            "slGeneStats"   "slGeneStats_in70"  "teGR"   "zGR"

# extract list elements that wil be frequently used
dmrs = dmrs_and_overlaps$dmrs;
dmrsOV = dmrs_and_overlaps$dmrsOV;
dmrsOVsimple = dmrs_and_overlaps$dmrsOVsimple;
numOV = dmrs_and_overlaps$numOV;
numOV2 = numOV[, -which(names(numOV) %in% c('gene','te','z'))];

hypodmrs = rownames(subset(dmrs, direction=='hypo'));
hyperdmrs = rownames(subset(dmrs, direction=='hyper'));

abHsMap = annotations_and_scaffold_stats$abHsMap;
annoCombo = annotations_and_scaffold_stats$annoCombo;

# # make other useful objects
# dmrsfc = log2(dmrs$group2.mean / dmrs$group1.mean);
# names(dmrsfc) = rownames(dmrs);
# 
# dmrsOVbyGeneSimple = sapply(dmrs_and_overlaps$dmrsOVbyGene$gene, rownames);
# 
# dmrsfcByGene = lapply(dmrsOVbyGeneSimple, function(f) dmrsfc[match(f, names(dmrsfc))]);


# compute GO and KEGG enrichments
# ---------------------------------

# build background 
bgGO = unique(unlist(strsplit(annoCombo$new$hsaHomologEntrez, ',')));
org.Hs.egGO2ALLEGSlist = as.list(org.Hs.egGO2ALLEGS);

# build gene lists to test
# ---------------------------------
feature_entrez_lists = .buildEntrezListsFromAbGeneLists(dmrsOVsimple, abHsMap, 
                                                        1, method='first.sorted');
dmrsOVsimple.hypo = lapply(dmrsOVsimple, 
                           function(f) f[names(f) %in% hypodmrs]);
dmrsOVsimple.hyper = lapply(dmrsOVsimple, 
                            function(f) f[names(f) %in% hyperdmrs]);
feature_entrez_lists.hypo = .buildEntrezListsFromAbGeneLists(dmrsOVsimple.hypo, abHsMap, 
                                                             1, method='first.sorted');
feature_entrez_lists.hyper = .buildEntrezListsFromAbGeneLists(dmrsOVsimple.hyper, abHsMap, 
                                                              1, method='first.sorted');

# GOFunction
# ---------------------------------
featureGO = .screenGOFunctionAllOntologies(feature_entrez_lists, bgGO, fdrmethvec=c('BY'));
featureGO.numhits = as.data.frame(sapply(featureGO, function(f) sapply(f, nrow)));

featureGO.hypo = .screenGOFunctionAllOntologies(feature_entrez_lists.hypo, bgGO, fdrmethvec=c('BY'));
featureGO.hypo.numhits = as.data.frame(sapply(featureGO.hypo, function(f) sapply(f, nrow)));

featureGO.hyper = .screenGOFunctionAllOntologies(feature_entrez_lists.hyper, bgGO, fdrmethvec=c('BY'));
featureGO.hyper.numhits = as.data.frame(sapply(featureGO.hyper, function(f) sapply(f, nrow)));

save(list=grep('featureGO|entrez',ls(),value=T),
     file='featureGO_and_input_entrez_lists.RData');

# enrichGO
# ---------------------------------
feature_enrichGO = .screen_enrichGOAllOntologies(feature_entrez_lists, bgGO, fdrmethvec=c('BY'));
feature_enrichGO.numhits = as.data.frame(sapply(feature_enrichGO, function(f) sapply(f, nrow)));

feature_enrichGO.hypo = .screen_enrichGOAllOntologies(feature_entrez_lists.hypo, bgGO, fdrmethvec=c('BY'));
feature_enrichGO.hypo.numhits = as.data.frame(sapply(feature_enrichGO.hypo, function(f) sapply(f, nrow)));

feature_enrichGO.hyper = .screen_enrichGOAllOntologies(feature_entrez_lists.hyper, bgGO, fdrmethvec=c('BY'));
feature_enrichGO.hyper.numhits = as.data.frame(sapply(feature_enrichGO.hyper, function(f) sapply(f, nrow)));

# enrichKEGG
# ---------------------------------
feature_enrichKEGG = .screen_enrichKEGG(feature_entrez_lists, bgGO);
feature_enrichKEGG.numhits = as.data.frame(do.call('rbind',lapply(feature_enrichKEGG, function(f) sapply(f, nrow))));

feature_enrichKEGG.hypo = .screen_enrichKEGG(feature_entrez_lists.hypo, bgGO);
feature_enrichKEGG.hypo.numhits = as.data.frame(do.call('rbind',lapply(feature_enrichKEGG.hypo, function(f) sapply(f, nrow))));

feature_enrichKEGG.hyper = .screen_enrichKEGG(feature_entrez_lists.hyper, bgGO);
feature_enrichKEGG.hyper.numhits = as.data.frame(do.call('rbind',lapply(feature_enrichKEGG.hyper, function(f) sapply(f, nrow))));


# # all genes associated with dmr in any way
# # ---------------------------------
# alldmrgenes_entrez = .buildEntrezListsFromAbGeneLists(list(alldmrgenes=gene_and_dmr_lists$alldmrgenes),
#                                                       abHsMap,
#                                                       1, method='first.sorted');
# alldmrgenesGO = .screenGOFunctionAllOntologies(alldmrgenes_entrez, bgGO, fdrmethvec=c('BY'));
# alldmrgenesGO.numhits = as.data.frame(sapply(alldmrgenesGO, 
#                                              function(f) sapply(f, nrow)));
# 
# # make a vector of all genes and brs overlapped by a dmr
# alldmrgenes = unique(unlist(dmrsOVsimple));
# 
# dmrsOVsimple.z = lapply(dmrsOVsimple,
#                          function(f) f[names(f) %in% rownames(subset(dmrs_and_overlaps$numOV,
#                                                                      z > 0))]);
# subset(dmrs_and_overlaps$numOV, z > 0)

#
# ---------------------------------
x = featureGO$gene$BY.05
xx  = .parseGOhitGenes(x, 8, ',', 3,feature_entrez_lists$gene);
xt = .buildGeneTermTable(xx);
names(xt) = paste0(1:ncol(xt), '-', names(xt))
tmp=heatmap(1-cor(xt), symm=T)

xx.byTerm2 = sapply(xx$byTerm, names)


dmrpca = prcomp(cor(xt), scale.=T, center=T); summary(dmrpca);
dmrpcacor=cor(cor(xt), dmrpca$x); 
dmrpcacorpval=WGCNA::corPvalueFisher(dmrpcacor,ncol(xt));
numpc = which(summary(dmrpca)$importance[3,] >= .9)[1];
tmp=dmrpcacor[,1:numpc];
tmp=tmp[order((tmp[,1])),];
tmpmat=signif(tmp,2)
tmppval=dmrpcacorpval[match(rownames(tmp), rownames(dmrpcacorpval)), 1:numpc];
tmpmat[which(!(tmppval < .05/nrow(tmp)^2), arr.ind=T)] = '';
par(oma=c(0,5,0,0));WGCNA::labeledHeatmap(tmp, 
                                          xLabels=colnames(dmrpca$x)[1:numpc], 
                                          yLabels=rownames(tmp),
                                          cex.lab.y = .5, 
                                          colors = WGCNA::blueWhiteRed(100),
                                          textMatrix=tmpmat,
                                          cex.text=.8);

as.data.frame(enrichKEGG(gene=feature_entrez_lists$gene, universe=bgGO));


# add information to GO results with info about genes in each term:
#  dmrs, feature type, expression fc, methylation fc, expression level
# ---------------------------------
inputgenes = feature_entrez_lists$gene;
thisGO = featureGO$gene$BY.05;

z = .getDmrsAndAbGenesForTermEntrezGenes(thisGO, inputgenes, dmrsOV$gene); 
dmrsByTerm = lapply(z, function(f) unique(f$queryHits));
genesByTerm = lapply(z, function(f) unique(f$mcols.geneSym));

fcounts = t(sapply(genesByTerm, 
                   function(ff) 
                     sapply(dmrsOVbyGene, 
                            function(f) 
                              length(na.omit(match(ff, names(f)))))));

thisGO = as.data.frame(cbind(thisGO, fcounts));

# load expression data
load('RData/deseq_methylationDvsND_res_withDmrInfo_noDESeq_input_or_raw_expression.RData');


