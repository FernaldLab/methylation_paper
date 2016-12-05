dirPrefix = 'CpGcombined.CG_4xCov'

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# repeat sections 5-8) - t-stat prep and calculation, get dmrs, plot dmrs
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# set variables that will stay constant throughout runs -
# for looping and filtering loci
nsVec=c(25,50,70); hVec=c(500,750,1000); 
# for t-stat calculation
reqCov=4;
groupND=c('3157_TENNISON','3677_MONK'); groupD=c('3165_BRISCOE','3581_LYNLEY');

groupNDnull1=c('3157_TENNISON','3165_BRISCOE'); groupDnull1=c('3677_MONK','3581_LYNLEY');
groupNDnull2=c('3157_TENNISON','3581_LYNLEY'); groupDnull2=c('3677_MONK','3165_BRISCOE');

# for dmr detection
centile=0.95; dmrMG=300; STAT='tstat.corrected';
# for dmr filtering
fCols=c('n','meanDiff'); fVals=c(3, 0.1); ord='areaStat';

for (ns in nsVec) {
  for (h in hVec) {
    # ----------------------------------------------------------------------------------------------------------
    # 5) load prepared smoothed data for t-statistic calculation
    # ----------------------------------------------------------------------------------------------------------
    smoothedDir = paste0(dirPrefix, '_ns', ns, 'h', h, 'mg', maxgap); #print(smoothedDir)

    FILEsm_filt = paste0(smoothedDir, '_fits_', reqCov, 'xCov.RData'); #print(FILEsm_filt)
#     save(cpgfits, cpgfits_f, cpg_forTtest, file=FILEsm_filt);
    cat(paste0('loading ',FILEsm_filt,'...\n'));
    load(FILEsm_filt);
#     
#     # ----------------------------------------------------------------------------------------------------------
#     # 6) compute t-statistics for each scaffold
#     # ----------------------------------------------------------------------------------------------------------
#     # run t-stat computation
#     # will use estimate.var="same" for all runs so no need to check methylation variation
#     # other parameters to BSmooth.tstat() that are hard-coded:
#     #  local.correct=T, mc.cores=4, verbose=F
    cat(paste0('computing t-stats...\n'));
    cat(paste0('  null1...\n'));
    cpgTstat_null1 = .BSmooth.tstatList(cpg_forTtest[-which(names(cpg_forTtest)=='scaffold_1713')], group1=groupNDnull1, group2=groupDnull1, estimate.var='same');
    cat(paste0('  null2...\n'));
    cpgTstat_null2 = .BSmooth.tstatList(cpg_forTtest[-which(names(cpg_forTtest)=='scaffold_1713')], group1=groupNDnull2, group2=groupDnull2, estimate.var='same');
#     
    # save out as RData object
    FILEt = gsub('.RData', '', FILEsm_filt);
    FILEt = paste0(FILEt, '_Ttests_nulls.RData');
    cat(paste0('saving ',FILEt,'...\n'));
    save(cpgTstat_null1, cpgTstat_null2, file=FILEt);
#     
#     # ----------------------------------------------------------------------------------------------------------
#     # 7) find differentially methylated regions (DMRs) on each scaffold
#     # ----------------------------------------------------------------------------------------------------------
    # get DMRs using default settings 
    cat(paste0('getting dmrs...\n'));
    cat(paste0('  null1...\n'));
    CUTnull1 = .getTstatQuantilesAcrossScaffolds(cpgTstat_null1$tstats, centile=centile);
    dmrsList_null1 = .dmrFinderList(cpgTstat_null1$tstats, cutoff=CUTnull1, maxGap=dmrMG, stat=STAT);
#     
    cat(paste0('  null2...\n'));
    CUTnull2 = .getTstatQuantilesAcrossScaffolds(cpgTstat_null2$tstats, centile=centile);
    dmrsList_null2 = .dmrFinderList(cpgTstat_null2$tstats, cutoff=CUTnull2, maxGap=dmrMG, stat=STAT);
#     
    # filter DMR list by n (number of loci) and meanDiff as suggested in BSmooth user guide
    cat(paste0('filtering dmrs...\n'));
    cat(paste0('  null1...\n'));
    dmrsList_f_null1 = .filterAndOrderDMRsList(dmrsList_null1$dmrList, filterOn=fCols, filterValues=fVals, orderOn=ord);
    dmrs_f_null1 = dmrsList_f_null1$filteredListCollapsed;
#     
    cat(paste0('  null2...\n'));
    dmrsList_f_null2 = .filterAndOrderDMRsList(dmrsList_null2$dmrList, filterOn=fCols, filterValues=fVals, orderOn=ord);
    dmrs_f_null2 = dmrsList_f_null2$filteredListCollapsed;
#      
    # save out as RData object
    filtParams = paste0(c(fCols, fVals), collapse='_');
    FILEdmrs = gsub('.RData', '', FILEt);
    FILEdmrs = paste0(FILEdmrs, '_DMRs_cut', centile, '_mg', dmrMG, '_', STAT, '_', filtParams,'_nulls.RData');
    cat(paste0('saving ',FILEdmrs,'...\n--------------------------------------------------------------------------\n\n'));
    save(dmrsList_f_null1, dmrs_f_null1, dmrsList_f_null2, dmrs_f_null2, file=FILEdmrs);
#     
#     # ----------------------------------------------------------------------------------------------------------
#     # 8) plot DMRs with t-stats
#     # ----------------------------------------------------------------------------------------------------------
#     plotDIR = paste(gsub('.RData', '', FILEdmrs), '_plots', sep='');
#     .plotDMRsAcrossList(DMRsList=dmrsList_f$filteredList, BSDList=cpg_forTtest, 
#                         BSD.TstatList=cpgTstat$tstats,
#                         colors=c('blue','red','red','blue'), extend=5000, addPoints=T, 
#                         DIR=plotDIR);
  }
} 



####################

rm(list=ls()); setwd('~/Documents/_BS-seq_analysis/');
nullDmrFiles = grep('DMRs.*nulls', list.files(), value=T);

# nullints = c();
# for (f in nullDmrFiles) {
#   load(f); rm(dmrsList_f_null1, dmrsList_f_null2);
#   nullints = unique(c(nullints, 
#                       c(paste0(dmrs_f_null1$chr, ':', dmrs_f_null1$start, '-', dmrs_f_null1$end),
#                         paste0(dmrs_f_null2$chr, ':', dmrs_f_null2$start, '-', dmrs_f_null2$end))));
# }; rm(f,dmrs_f_null1, dmrs_f_null2);
# 
# nullintssplit = strsplit(nullints, ':');
# nullintschr = sapply(nullintssplit, function(f) f[1]);
# nullintssplit = strsplit(sapply(nullintssplit, function(f) f[2]), '-');
# nullintsstart = as.numeric(sapply(nullintssplit, function(f) f[1]));
# nullintsend = as.numeric(sapply(nullintssplit, function(f) f[2]));
# nullGR = GRanges(seqnames=nullintschr, ranges=IRanges(start=nullintsstart, end=nullintsend));


pattern = 'mg1e+08_fits_4xCov_Ttests_DMRs_cut0.95_mg300_tstat.corrected_n3_meanDiff0.1.RData';
dmrFiles = grep(pattern, list.files(), value=T, fixed=T);
dmrFiles = grep('All|WORKSPACE', dmrFiles, value=T, invert=T);

dmrFilesRuns = sapply(strsplit(dmrFiles, '_'), function(f) f[3]);
nullDmrFilesRuns = sapply(strsplit(nullDmrFiles, '_'), function(f) f[3]);

if (all(dmrFilesRuns  ==  nullDmrFilesRuns)) {
  dmrf = list();
  nullints = list();
  for (f in 1:length(dmrFilesRuns)) {
    cat(dmrFilesRuns[f],'...')
    load(dmrFiles[f]);
    load(nullDmrFiles[f]);
    
    dmrf[[dmrFilesRuns[f]]] = dmrs_f;
    
    nullints[[dmrFilesRuns[f]]] = unique(c(paste0(dmrs_f_null1$chr, ':', dmrs_f_null1$start, '-', dmrs_f_null1$end), 
                                          paste0(dmrs_f_null2$chr, ':', dmrs_f_null2$start, '-', dmrs_f_null2$end)));
    
  }
}
rm(f, dmrs_f, dmrsList_f, dmrsList_f_null1, dmrsList_f_null2, dmrs_f_null1, dmrs_f_null2);


dmrf2 = dmrf;
for (i in 1:length(dmrf)) {
  dd = dmrf[[i]];
  dmrf2[[i]] = dd[!(paste0(dd$chr, ':', dd$start, '-', dd$end) %in% nullints[[i]]), ];
}
rm(i,dd);

dmrf0 = dmrf;
dmrf = dmrf2;

##################
allnullints = unique(unlist(nullints));


.filter_nullints_fromDMRs = function (dmrs, allnullints) {
  nov = list();
  dmrsints = paste0(dmrs$chr, ':', dmrs$start, '-', dmrs$end);
  allnullintsGR=.buildGR_for_intervals(allnullints);
  mcolCols = which(names(dmrs)=='n'):ncol(dmrs);
  dmrsGR = GRanges(seqnames=dmrs$chr, 
                   ranges=IRanges(start=dmrs$start, end=dmrs$end), 
                   mcols=dmrs[,mcolCols]);
  for (d in 1:nrow(dmrs)) {
    thisw = dmrs$width[d];
    nov[[dmrsints[d]]] = as.data.frame(findOverlaps(dmrsGR[d], allnullintsGR, minoverlap=(thisw/2)))
  }
  
  return(dmrs[!(dmrsints %in% names(nov)[sapply(nov, nrow)!=0]), ])
}


nov = list()
for (d in 1:nrow(dmrs)) {
  thisw = dmrs$width[d];
  nov[[dmrsints[d]]] = as.data.frame(findOverlaps(dmrsGR[d], allnullintsGR, minoverlap=(thisw/2)))
}; rm(d,thisw);

dmrs = dmrs[!(dmrsints %in% names(nov)[sapply(nov, nrow)!=0]), ];
dmrsints = paste0(dmrs$chr, ':', dmrs$start, '-', dmrs$end);
mcolCols = which(names(dmrs)=='n'):ncol(dmrs);
dmrsGR = GRanges(seqnames=dmrs$chr, ranges=IRanges(start=dmrs$start, end=dmrs$end), mcols=dmrs[,mcolCols]);