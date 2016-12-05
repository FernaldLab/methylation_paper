# clear workspace, set options and working directory, load functions
rm(list=ls()); options(stringsAsFactors=F); 
setwd('~/Documents/_BS-seq_analysis/');
source('/Volumes/fishstudies-1/_code/run_bsmooth_functions.R');
library(org.Hs.eg.db);
library(GOFunction);

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
dmrsOVsimple = dmrs_and_overlaps$dmrsOVsimple;

hypodmrs = rownames(subset(dmrs, direction=='hypo'));
hyperdmrs = rownames(subset(dmrs, direction=='hyper'));

abHsMap = annotations_and_scaffold_stats$abHsMap;
annoCombo = annotations_and_scaffold_stats$annoCombo;


# compute GO and KEGG enrichments

# build background 
bgGO = unique(unlist(strsplit(annoCombo$new$hsaHomologEntrez, ',')));
org.Hs.egGO2ALLEGSlist = as.list(org.Hs.egGO2ALLEGS);

# build gene lists to test
# each feature type separately
feature_entrez_lists = .buildEntrezListsFromAbGeneLists(dmrsOVsimple, abHsMap, 
                                                        1, method='first.sorted');
featureGO = .screenGOFunctionAllOntologies(feature_entrez_lists, bgGO, fdrmethvec=c('BY'));
featureGO.numhits = as.data.frame(sapply(featureGO, 
                                         function(f) sapply(f, nrow)));

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


dmrsOVsimple.z = lapply(dmrsOVsimple, 
                        function(f) f[names(f) %in% rownames(subset(dmrs_and_overlaps$numOV, 
                                                                    z > 0))]);
  
  subset(dmrs_and_overlaps$numOV, z > 0)

x = featureGO$gene$BY.05
xx  = .parseGOhitGenes(x, 8, ',', 3,feature_entrez_lists$gene);
xt = .buildGeneTermTable(xx);
names(xt) = paste0(1:ncol(xt), '-', names(xt))
heatmap(1-cor(xt), symm=T)

xx.byTerm2 = sapply(xx$byTerm, names)


summary(clusterProfiler::enrichKEGG(gene=feature_entrez_lists$gene, universe=bgGO))