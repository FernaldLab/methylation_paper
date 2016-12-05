
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






# get distances to genes for shuffled dmrs
shuffDistToGenes = lapply(as.list(shuffGRlist[[1]]), function(f) .distance.named(geneGR, f));


#
shuffDistToGenes_summary = as.data.frame(cbind(do.call('rbind', 
                                                      lapply(shuffDistToGenes, summary)), 
                                              num=sapply(shuffDistToGenes, length)));

#
shuffDistToGenes1mb = lapply(shuffDistToGenes, function(f) f[f <= 1e6]);

#
shuffDistToGenes1mb_summary = as.data.frame(cbind(do.call('rbind', 
                                                       lapply(shuffDistToGenes1mb, summary)), 
                                               num=sapply(shuffDistToGenes1mb, length)));

#
shuffDistToGenes50kb = lapply(shuffDistToGenes, function(f) f[f <= 50000]);

#
shuffDistToGenes50kb_summary = as.data.frame(cbind(do.call('rbind', 
                                                          lapply(shuffDistToGenes50kb, summary)), 
                                                  num=sapply(shuffDistToGenes50kb, length)));

#########################
shuffDistToGenesAcrossRuns = lapply(shuffGRlist, 
                                    function(ff) lapply(as.list(ff), 
                                                        function(f) .distance.named(geneGR, f)));

shuffDistToGenesAcrossRunsNumOverlaps = lapply(shuffDistToGenesAcrossRuns, 
                                               function(ff) table(sapply(ff, 
                                                                         function(f) sum(f==0))));
tmp = sapply(shuffDistToGenesAcrossRunsNumOverlaps, length);
shuffDistToGenesAcrossRunsNumOverlaps_df = as.data.frame(do.call('rbind',shuffDistToGenesAcrossRunsNumOverlaps));
shuffDistToGenesAcrossRunsNumOverlaps_df$'3'[tmp!=4] = 0; rm(tmp)

shuffDistToGenesAcrossRuns50kb = lapply(shuffDistToGenesAcrossRuns, 
                                        function(f) lapply(f, 
                                                           function(ff) ff[ff <= 50000]));

shuffDistToGenesAcrossRuns50kb_summary = lapply(shuffDistToGenesAcrossRuns50kb, 
                                                function(f) as.data.frame(cbind(do.call('rbind', 
                                                                                        lapply(f, summary)), 
                                                                                num=sapply(f, length))));

shuffDistToGenesAcrossRuns50kb_summary_df = do.call('rbind', 
                                                    lapply(shuffDistToGenesAcrossRuns50kb_summary, 
                                                           function(f) data.frame(dist.mean=mean(f$Mean,na.rm=T), 
                                                                                  dist.median=median(f$Median,na.rm=T), 
                                                                                  num.mean=mean(f$num), 
                                                                                  num.median=median(f$num))));

shuffDistToGenesAcrossRuns100kb = lapply(shuffDistToGenesAcrossRuns, 
                                        function(f) lapply(f, 
                                                           function(ff) ff[ff <= 100000]));

shuffDistToGenesAcrossRuns100kb_summary = lapply(shuffDistToGenesAcrossRuns100kb, 
                                                function(f) as.data.frame(cbind(do.call('rbind', 
                                                                                        lapply(f, summary)), 
                                                                                num=sapply(f, length))));

shuffDistToGenesAcrossRuns100kb_summary_df = do.call('rbind', 
                                                    lapply(shuffDistToGenesAcrossRuns100kb_summary, 
                                                           function(f) data.frame(dist.mean=mean(f$Mean,na.rm=T), 
                                                                                  dist.median=median(f$Median,na.rm=T), 
                                                                                  num.mean=mean(f$num), 
                                                                                  num.median=median(f$num))));

shuffDistToGenesAcrossRuns250kb = lapply(shuffDistToGenesAcrossRuns, 
                                         function(f) lapply(f, 
                                                            function(ff) ff[ff <= 250000]));

shuffDistToGenesAcrossRuns250kb_summary = lapply(shuffDistToGenesAcrossRuns250kb, 
                                                 function(f) as.data.frame(cbind(do.call('rbind', 
                                                                                         lapply(f, summary)), 
                                                                                 num=sapply(f, length))));

shuffDistToGenesAcrossRuns250kb_summary_df = do.call('rbind', 
                                                     lapply(shuffDistToGenesAcrossRuns250kb_summary, 
                                                            function(f) data.frame(dist.mean=mean(f$Mean,na.rm=T), 
                                                                                   dist.median=median(f$Median,na.rm=T), 
                                                                                   num.mean=mean(f$num), 
                                                                                   num.median=median(f$num))));

shuffDistToGenesAcrossRuns500kb = lapply(shuffDistToGenesAcrossRuns, 
                                         function(f) lapply(f, 
                                                            function(ff) ff[ff <= 500000]));

shuffDistToGenesAcrossRuns500kb_summary = lapply(shuffDistToGenesAcrossRuns500kb, 
                                                 function(f) as.data.frame(cbind(do.call('rbind', 
                                                                                         lapply(f, summary)), 
                                                                                 num=sapply(f, length))));

shuffDistToGenesAcrossRuns500kb_summary_df = do.call('rbind', 
                                                     lapply(shuffDistToGenesAcrossRuns500kb_summary, 
                                                            function(f) data.frame(dist.mean=mean(f$Mean,na.rm=T), 
                                                                                   dist.median=median(f$Median,na.rm=T), 
                                                                                   num.mean=mean(f$num), 
                                                                                   num.median=median(f$num))));


shuffDistToGenesAcrossRuns1mb = lapply(shuffDistToGenesAcrossRuns, 
                                       function(f) lapply(f, 
                                                          function(ff) ff[ff <= 1e6]));

shuffDistToGenesAcrossRuns1mb_summary = lapply(shuffDistToGenesAcrossRuns1mb, 
                                                function(f) as.data.frame(cbind(do.call('rbind', 
                                                                                        lapply(f, summary)), 
                                                                                num=sapply(f, length))));

shuffDistToGenesAcrossRuns1mb_summary_df = do.call('rbind', 
                                                    lapply(shuffDistToGenesAcrossRuns1mb_summary, 
                                                           function(f) data.frame(dist.mean=mean(f$Mean,na.rm=T), 
                                                                                  dist.median=median(f$Median,na.rm=T), 
                                                                                  num.mean=mean(f$num), 
                                                                                  num.median=median(f$num))));


##############
load('null_GR_OV_DMRstats_Oct2016.RData')
nullDistToGenesAcrossRuns = lapply(nullGRlist, 
                                   function(ff) lapply(as.list(ff), 
                                                       function(f) .distance.named(geneGR, f)));

nullDistToGenesAcrossRunsNumOverlaps = lapply(nullDistToGenesAcrossRuns, 
                                              function(ff) table(sapply(ff, 
                                                                        function(f) sum(f==0))));

tmp = sapply(nullDistToGenesAcrossRunsNumOverlaps, length);
nullDistToGenesAcrossRunsNumOverlaps_df = as.data.frame(do.call('rbind',nullDistToGenesAcrossRunsNumOverlaps));
for (i in (1:nrow(nullDistToGenesAcrossRunsNumOverlaps_df))[-which(tmp==ncol(nullDistToGenesAcrossRunsNumOverlaps_df))]) {
  if (i %% 1000 == 0) { cat(i,'...') }
  numz = ncol(nullDistToGenesAcrossRunsNumOverlaps_df) - tmp[i];
  nullDistToGenesAcrossRunsNumOverlaps_df[i, (tmp[i]+1):ncol(nullDistToGenesAcrossRunsNumOverlaps_df)] = rep(0, numz);
}; rm(i);



