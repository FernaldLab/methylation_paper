# clear workspace, set options and working directory, load functions
rm(list=ls()); options(stringsAsFactors=F); 
setwd('~/Documents/_BS-seq_analysis/');
source('/Volumes/fishstudies-1/_code/run_bsmooth_functions.R');
library(org.Hs.eg.db);
library(GOFunction);

# ---------------------------------
# load lists from end of run_bsmooth_for_results03_get_overlaps.R
for (f in grep('after03',list.files(),value=T)) { load(f) }; rm(f);

for (l in ls()) { cat(l,':\n',names(get(l)),'\n') }; rm(l);
# annotations_and_scaffold_stats :
#   abHsMap an annoCombo burtoniGois slGeneStats slGeneStats_in70 
# dmrs_and_overlaps :
#   dmrs dmrsOV dmrsOVbyGene dmrsOVgenebrdf dmrsOVgenebrdf1 dmrsOVsimple numOV fg 
# gene_and_dmr_lists :
#   alldmrgenes multiDmrGenes multiGeneDmrs multigenes 
# grs :
#   br3primeGR brGR cdsGR dmrsGR exonGR geneGR intronGR teGR zGR 

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


# compute GO and KEGG enrichments
# ---------------------------------

# build background 
bgGO = unique(unlist(strsplit(annoCombo$new$hsaHomologEntrez, ',')));
org.Hs.egGO2ALLEGSlist = as.list(org.Hs.egGO2ALLEGS);

# build gene lists to test

# each feature type separately
# ---------------------------------

feature_entrez_lists = .buildEntrezListsFromAbGeneLists(dmrsOVsimple, abHsMap, 
                                                        1, method='first.sorted');
featureGO = .screenGOFunctionAllOntologies(feature_entrez_lists, bgGO, fdrmethvec=c('BY'));
featureGO.numhits = as.data.frame(sapply(featureGO, 
                                         function(f) sapply(f, nrow)));

feature_enrichGO = .screen_enrichGOAllOntologies(feature_entrez_lists, bgGO, fdrmethvec=c('BY'))




# each feature type separately for hypo vs hyper dmrs
dmrsOVsimple.hypo = lapply(dmrsOVsimple, 
                           function(f) f[names(f) %in% hypodmrs]);
dmrsOVsimple.hyper = lapply(dmrsOVsimple, 
                            function(f) f[names(f) %in% hyperdmrs]);

feature_entrez_lists.hypo = .buildEntrezListsFromAbGeneLists(dmrsOVsimple.hypo, abHsMap, 
                                                             1, method='first.sorted');
featureGO.hypo = .screenGOFunctionAllOntologies(feature_entrez_lists.hypo, bgGO, fdrmethvec=c('BY'));
featureGO.hypo.numhits = as.data.frame(sapply(featureGO.hypo, 
                                              function(f) sapply(f, nrow)));

feature_entrez_lists.hyper = .buildEntrezListsFromAbGeneLists(dmrsOVsimple.hyper, abHsMap, 
                                                              1, method='first.sorted');
featureGO.hyper = .screenGOFunctionAllOntologies(feature_entrez_lists.hyper, bgGO, fdrmethvec=c('BY'));
featureGO.hyper.numhits = as.data.frame(sapply(featureGO.hyper, 
                                               function(f) sapply(f, nrow)));

save(list=grep('featureGO|entrez',ls(),value=T),
     file='featureGO_and_input_entrez_lists.RData');





# all genes associated with dmr in any way
# ---------------------------------
alldmrgenes_entrez = .buildEntrezListsFromAbGeneLists(list(alldmrgenes=gene_and_dmr_lists$alldmrgenes), 
                                                      abHsMap, 
                                                      1, method='first.sorted');
alldmrgenesGO = .screenGOFunctionAllOntologies(alldmrgenes_entrez, bgGO, fdrmethvec=c('BY'));
alldmrgenesGO.numhits = as.data.frame(sapply(alldmrgenesGO, 
                                             function(f) sapply(f, nrow)));

# # make a vector of all genes and brs overlapped by a dmr
# alldmrgenes = unique(unlist(dmrsOVsimple));

# dmrsOVsimple.z = lapply(dmrsOVsimple, 
#                         function(f) f[names(f) %in% rownames(subset(dmrs_and_overlaps$numOV, 
#                                                                     z > 0))]);
#   
#   subset(dmrs_and_overlaps$numOV, z > 0)

#
# ---------------------------------
x = featureGO$gene$BY.1
xx  = .parseGOhitGenes(x, 8, ',', 3,feature_entrez_lists$gene);
xt = .buildGeneTermTable(xx);
names(xt) = paste0(1:ncol(xt), '-', names(xt))
heatmap(1-cor(xt), symm=T)

xx.byTerm2 = sapply(xx$byTerm, names)


as.data.frame(enrichKEGG(gene=feature_entrez_lists$gene, universe=bgGO));


# add information to GO results with info about genes in each term:
#  dmrs, feature type, expression fc, methylation fc, expression level
# ---------------------------------
inputgenes = feature_entrez_lists$gene;
thisGO = featureGO$gene$BY.1;

z = .getDmrsAndAbGenesForTermEntrezGenes(thisGO, inputgenes, dmrsOV$gene); 
dmrsByTerm = lapply(z, function(f) unique(f$queryHits));
genesByTerm = lapply(z, function(f) unique(f$mcols.geneSym));

fcounts = data.frame(br3=sapply(z, function(f) length(unique(dmrsOV$br3prime$mcols.geneSym[na.omit(match(f$queryHits,dmrsOV$br3prime$queryHits))]))),
                     br=sapply(z, function(f) length(unique(dmrsOV$br$mcols.geneSym[na.omit(match(f$queryHits,dmrsOV$br$queryHits))]))),
                     cds=sapply(z, function(f) length(unique(dmrsOV$cds$mcols.geneSym[na.omit(match(f$queryHits,dmrsOV$cds$queryHits))]))),
                     exon=sapply(z, function(f) length(unique(dmrsOV$exon$mcols.geneSym[na.omit(match(f$queryHits,dmrsOV$exon$queryHits))]))),
                     intron=sapply(z, function(f) length(unique(dmrsOV$intron$mcols.geneSym[na.omit(match(f$queryHits,dmrsOV$intron$queryHits))]))));

dcounts = data.frame(dmrs=sapply(z, function(f) length(unique(f$queryHits))),
                     hypo=sapply(z, function(f) sum(unique(f$queryHits) %in% hypodmrs)),
                     hyper=sapply(z, function(f) sum(unique(f$queryHits) %in% hyperdmrs)));

hypofisher = .checkGeneListEnrichmentList(hypodmrs, dmrsByTerm, rownames(dmrs), order=F)$pvals
names(hypofisher) = paste0('hypo.fisher_',names(hypofisher));

ffisher = do.call('cbind', lapply(dmrsOVsimple, 
                                  function(f) .checkGeneListEnrichmentList(names(f), dmrsByTerm, rownames(dmrs), order=F)$pvals));
# use unique(unlist(dmrsByTerm)) for background?

numOVfisher = do.call('cbind',lapply(apply(numOV[-5], 2, 
                                           function(f) rownames(subset(numOV[-5], f>0))), 
                                     function(f) .checkGeneListEnrichmentList(f, dmrsByTerm, rownames(dmrs), order=F)$pvals))

dstats = do.call('rbind', lapply(z, 
                                 function(f) (.computeStatsForDMRs(dmrs[match(unique(f$queryHits), rownames(dmrs)), ], 
                                                                   c(7:15,17))$statmat[4,])));
dmrs.hypo = subset(dmrs, direction=='hypo');
dstats.hypo = do.call('rbind', lapply(z, 
                                      function(f) (.computeStatsForDMRs(dmrs.hypo[na.omit(match(unique(f$queryHits), rownames(dmrs.hypo))), ], 
                                                                        c(7:15,17))$statmat[4,])));
dstats.hypo = dstats.hypo[, !grepl('abs',names(dstats.hypo))];

dmrs.hyper = subset(dmrs, direction=='hyper');
# need to skip term 33 because only 1 hyper dmr (if using featureGO$gene$BY.1)
dstats.hyper = do.call('rbind', lapply(z[-33], 
                                      function(f) (.computeStatsForDMRs(dmrs.hyper[na.omit(match(unique(f$queryHits), rownames(dmrs.hyper))), ], 
                                                                        c(7:15,17))$statmat[4,])));
dstats.hyper = dstats.hyper[, !grepl('abs',names(dstats.hyper))];
dstats.hyper = as.data.frame(rbind(dstats.hyper, 
                                   dmrs[rownames(dmrs)==z[[33]]$queryHits[z[[33]]$queryHits %in% hyperdmrs],
                                        c(7:15,17)]));
rownames(dstats.hyper)[33] = rownames(dstats.hypo)[33];

par(mfrow=c(2,5));for(i in 1:ncol(dstats.hypo)) { 
  boxplot(abs(dstats.hypo[,i]), abs(dstats.hyper[,i]), notch=T, names=c('hypo','hyper'), 
          main=paste0(names(dstats.hypo)[i],'\n',
                      signif(kruskal.test(list(abs(dstats.hypo[,i]), abs(dstats.hyper[,i])))$p.value, 3))); 
  abline(h=median(dstats[,match(names(dstats.hypo)[i], names(dstats))]), col='red')  
}; rm(i)

#
dmrsfc = log2(dmrs$group2.mean / dmrs$group1.mean);
names(dmrsfc) = rownames(dmrs);

dfc = lapply(z, 
             function(f) log2(dmrs$group2.mean[match(unique(f$queryHits), 
                                                     rownames(dmrs))] / dmrs$group1.mean[match(unique(f$queryHits), 
                                                                                               rownames(dmrs))]));
dfc = data.frame(dmrfc=sapply(dfc, mean),
                 dmrfc.abs=sapply(dfc, function(f) mean(abs(f))));

dfc.hypo = lapply(z, 
                  function(f) log2(dmrs.hypo$group2.mean[match(unique(f$queryHits), rownames(dmrs.hypo))] / dmrs.hypo$group1.mean[match(unique(f$queryHits), 
                                                                                               rownames(dmrs.hypo))]));

dfc.hyper = lapply(z, 
                  function(f) log2(subset(dmrs,direction=='hyper')$group2.mean[match(unique(f$queryHits), rownames(subset(dmrs,direction=='hyper')))] / subset(dmrs,direction=='hyper')$group1.mean[match(unique(f$queryHits), 
                                                                                                                                                                                                       rownames(subset(dmrs,direction=='hyper')))]));

par(mfrow=c(3,9));
for(i in 1:length(dfc.hypo)) {
  boxplot(list(abs(dfc.hypo[[i]]), 
               abs(dfc.hyper[[i]])), 
          notch=T, 
          main=paste0(names(z)[i], '\n',
                      signif(kruskal.test(list(abs(dfc.hypo[[i]]), 
                                               abs(dfc.hyper[[i]])))$p.value, 3)), 
          names=c('hypo','hyper'))
}; rm(i)

dfc.hypo = data.frame(dmrfc.hypo=sapply(dfc.hypo, mean, na.rm=T),
                 dmrfc.hypo.abs=sapply(dfc.hypo, function(f) mean(abs(f),na.rm=T)));
# gene_part_cols = which(names(numOV) %in% c('br3prime','br','cds','exon','intron'));
# for (i in gene_part_cols) {
#   print(rownames(numOV)[which(  apply(numOV, 1, function(f) sum(f[gene_part_cols])==1 & f[i]==1  )   )])
# }; rm(i);

#thisGO = as.data.frame(cbind(thisGO,fcounts,dcounts,hypofisher,numOVfisher));
thisGO = as.data.frame(cbind(thisGO,fcounts,dcounts,dstats,dfc));

# load expression data

