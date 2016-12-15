# ==========================================================================================================
# ==========================================================================================================
# run BSmooth with default settings on symmetric CpG data
# ==========================================================================================================
# ==========================================================================================================
# files named "aligned_trimmed4-98.adapters.q30.m0_bsmap2.9.bam_methratio_samtools0.1.19-m5-CpGcombined.CG" 
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# 1) setup
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# clear workspace, set options and working directory
rm(list=ls()); options(stringsAsFactors=F); setwd('~/Documents/_BS-seq_data/BROAD_genome/');

# load required functions
library('bsseq'); source('/Volumes/fishstudies-1/_code/run_bsmooth_functions.R');

# # create directory to hold processed data
# dir.create('../_BS-seq_analysis');

# set names of subjects and whether they were dominant or non-dominant
subjects = list.files()[grep('^3', list.files())]; 
groups = c('ND','D','D','ND'); 
names(groups) = subjects;

# set filename to load output of methratio.py for each subject
file_cg = 'aligned_trimmed4-98.adapters.q30.m0_bsmap2.9.bam_methratio_samtools0.1.19-m4-CpGcombined.CG';

# define classes and names of columns to load 
colClasses = c('character', 'numeric', rep('NULL',3), rep('numeric',2), rep('NULL',5));
colnames = c('chr','pos','Cov','M');

# set coverage required for a locus to be included in analysis
reqCov = 4;

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# 2) load data for each subject and convert to BSseq objects
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# make a list where each element is a dataframe containing methratio.py output for a subject
dfList = .loadMethRatioData(filename=file_cg, subjects=subjects, 
                            colClasses=colClasses, colNames=colnames, 
                            checkTime=T, header=F);

# convert dataframes into BSseq objects
bsdList = .makeBSseqListFromDFList(dfList=dfList, reqCov=reqCov, groups=groups, checkTime=T);

# combine all subject BSseq objects into a single object
bsdAll = .combineBSDList(bsdList=bsdList, checkTime=T);

# save out as R objects
save(dfList, file=paste0('../../_BS-seq_analysis/', file_cg, '_dfList.RData'));
save(bsdList, file=paste0('../../_BS-seq_analysis/', file_cg, '_bsdList.RData'));
save(bsdAll, file=paste0('../../_BS-seq_analysis/', file_cg, '_bsdAll.RData'));

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# 3) prepare for smoothing
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# reset workspace and switch working directory
rm(list=ls()); setwd('~/Documents/_BS-seq_analysis/');

# load raw symmetric CpG data from saved BSseq object containing all 4 animals
file_cg = 'aligned_trimmed4-98.adapters.q30.m0_bsmap2.9.bam_methratio_samtools0.1.19-m4-CpGcombined.CG';
load(paste0(file_cg, '_bsdAll.RData'));

# filter to loci that had non-zero coverage in all 4 subjects
# in effect, removes loci with <4x coverage in all subjects since .makeBSseqListFromDFList ran with reqCov=4 
#  and 0's were only introduced into the combined BSseq object made with .combineBSDList
cpgcommon = .filterLowCoverage2(bsdAll, reqCov=1, fix.seqlevels=T);

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# 4) smooth each scaffold separately using default BSmooth parameters
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# save smoothed scaffolds into a directory
# skip scaffolds with <70 CpG loci by setting skipLowNS=T
#  if skipLowNS=F, would reduce the ns parameter to smooth scaffolds that would otherwise be skipped
#                  these scaffolds would be saved into subdirectory '_modified_ns' 
ns=70; h=1000; maxgap=10^8; dirPrefix='CpGcombined.CG_4xCov';
.BSmoothAcrossScaffolds(BSD=cpgcommon, NS=ns, H=h, MAXGAP=maxgap, skipLowNS=T, dirPrefix=dirPrefix);

# ----------------------------------------------------------------------------------------------------------
# may need to do the following if .BSmoothAcrossScaffolds hangs...
# ----------------------------------------------------------------------------------------------------------
# sometimes in RStudio .BSmoothAcrossScaffolds will hang for unknown reasons
# I'm sure the function could be optimized and written to be more robust but here is a workaround for now...
#
# # .filterBeyondScaffold creates a subset of a BSseq object by removing all scaffolds below a given number
# # e.g. if .BSmoothAcrossScaffolds hangs on scaffold_1847 would run:
# cpgcommonPostRemaining = .filterBeyondScaffold(cpgcommon, 1847);
# .BSmoothAcrossScaffolds(BSD=cpgcommonPostRemaining, NS=ns, H=h, MAXGAP=maxgap);
#
# cpgcommonPostRemaining needs to be created so the original output directory isn't overwritten
# could also use dirPrefix to avoid overwritting original directory
# either way will need to copy/paste smoothed scaffolds into original directory
#
# # however, as of June2016 .filterBeyondScaffold has weird behavior if scaffold number is <1000
# # if .BSmoothAcrossScaffolds hangs on a scaffold under scaffold_1000 should do the following:
# # set name of directory with smoothed scaffolds
# smoothedDir = paste0(dirPrefix, '_ns', ns, 'h', h, 'mg', maxgap);
# # get names of scaffolds that have already been smoothed
# smoothedScaffolds0 = gsub('.RData', '', grep('^scaffold', list.files(smoothedDir), value=T));
# # in case skipLowNS was FALSE
# smoothedScaffoldsModified = grep('^scaffold', 
#                                  unlist(strsplit(list.files(paste0(smoothedDir,'/_modified_ns')), '_ns')), 
#                                  value=T);
# smoothedScaffolds = c(smoothedScaffolds0, smoothedScaffoldsModified);
# # get names of scaffolds in the BSseq object that have yet to be smoothed
# if (!all(seqlevels(cpgcommon) %in% smoothedScaffolds)) {
#   scaffoldsToSmooth = seqlevels(cpgcommon)[which(!(seqlevels(cpgcommon) %in% smoothedScaffolds))];
# }
# # create subset of the BSseq object with remaining scaffolds to smooth
# cpgcommonStillToSmooth = .fixSeqLevels(cpgcommon[seqnames(cpgcommon) %in% scaffoldsToSmooth]);
# .BSmoothAcrossScaffolds(BSD=cpgcommonStillToSmooth, NS=NS, H=H, MAXGAP=MAXGAP, skipLowNS=T);
# ----------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# 5) prepare smoothed data for t-statistic calculation
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# load files output by .BSmoothAcrossScaffolds(), print progress every 200 files
# will generate a list where each element is a smoothed scaffold stored in a BSSeq object
smoothedDir = paste0(dirPrefix, '_ns', ns, 'h', h, 'mg', maxgap);
cpgfits = .loadSmoothedScaffolds(smoothedDir, printEvery=200);

# filter to loci with at least 4x coverage in every subject
# is actually redundant since <4x loci were filtered above, 
#  but will run here for completeness and to double-check information on scaffolds with <ns loci
reqCov=4; 
cpgfits_f = .filterLowCoverageList(cpgfits, reqCov=reqCov, ns.BSmooth=ns);

# if all smoothed scaffolds had at least ns loci,
#  get the filtered list from the .filterLowCoverageList output
# otherwise remove scaffolds with fewer loci than BSmooth ns parameter 
if (is.null(cpgfits_f$numLociBelowNS)) {
  cpg_forTtest = cpgfits_f$filtered;
} else {
  toRemove = match(cpgfits_f$numLociBelowNS, names(cpgfits_f$filtered));
  cpg_forTtest = cpgfits_f$filtered[-toRemove];
}

# save out as RData object
FILEsm_filt = paste0(smoothedDir, '_fits_', reqCov, 'xCov.RData');
save(cpgfits, cpgfits_f, cpg_forTtest, file=FILEsm_filt);
#load(FILEsm_filt)

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# 6) compute t-statistics for each scaffold
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# reset workspace and re-load the smoothed data
rm(list=ls());
ns=70; h=1000; maxgap=10^8; dirPrefix='CpGcombined.CG_4xCov'; reqCov=4; 
smoothedDir = paste0(dirPrefix, '_ns', ns, 'h', h, 'mg', maxgap);
FILEsm_filt = paste0(smoothedDir, '_fits_', reqCov, 'xCov.RData');
load(FILEsm_filt);

# decide which group to use for variance estimation
methCVraw = .getMethylationCVList(cpg_forTtest, type='raw');
methCVsmooth = .getMethylationCVList(cpg_forTtest, type='smooth');
apply(methCVraw, 2, mean, na.rm=T);
apply(methCVsmooth, 2, mean, na.rm=T);
# will use estimate.var="same" since all roughly equal

# define groups based on sampleNames
groupND = c('3157_TENNISON','3677_MONK');
groupD = c('3165_BRISCOE','3581_LYNLEY');
 
# run t-stat computation
# other parameters to BSmooth.tstat() that are hard-coded:
#  local.correct=T, mc.cores=4, verbose=F
cpgTstat = .BSmooth.tstatList(cpg_forTtest, group1=groupND, group2=groupD, estimate.var='same');

# save out as RData object
FILEt = gsub('.RData', '', FILEsm_filt);
FILEt = paste0(FILEt, '_Ttests.RData');
save(cpgTstat, file=FILEt);

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# 7) find differentially methylated regions (DMRs) on each scaffold
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# reset workspace and re-load the smoothed data and t-statistics
rm(list=ls());
ns=70; h=1000; maxgap=10^8; dirPrefix='CpGcombined.CG_4xCov'; reqCov=4; 
smoothedDir = paste0(dirPrefix, '_ns', ns, 'h', h, 'mg', maxgap);
FILEsm_filt = paste0(smoothedDir, '_fits_', reqCov, 'xCov.RData');
FILEt = gsub('.RData', '', FILEsm_filt);
FILEt = paste0(FILEt, '_Ttests.RData');
load(FILEsm_filt);
load(FILEt);

# get DMRs using default settings 
centile=0.95; dmrMG=300; STAT='tstat.corrected';
CUT = .getTstatQuantilesAcrossScaffolds(cpgTstat$tstats, centile=centile);
dmrsList = .dmrFinderList(cpgTstat$tstats, cutoff=CUT, maxGap=dmrMG, stat=STAT);

# filter DMR list by n (number of loci) and meanDiff as suggested in BSmooth user guide
fCols=c('n','meanDiff'); fVals=c(3, 0.1); ord='areaStat';
dmrsList_f = .filterAndOrderDMRsList(dmrsList$dmrList, filterOn=fCols, filterValues=fVals, orderOn=ord);
dmrs_f = dmrsList_f$filteredListCollapsed;

# save out as RData object
filtParams = paste0(c(fCols, fVals), collapse='_');
FILEdmrs = gsub('.RData', '', FILEt);
FILEdmrs = paste0(FILEdmrs, '_DMRs_cut', centile, '_mg', dmrMG, '_', STAT, '_', filtParams,'.RData');
save(dmrsList_f, dmrs_f, file=FILEdmrs);

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# 8) plot DMRs with t-stats
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
#load(FILEdmrs);
plotDIR = paste(gsub('.RData', '', FILEdmrs), '_plots', sep='');
.plotDMRsAcrossList(DMRsList=dmrsList_f$filteredList, BSDList=cpg_forTtest, 
                    BSD.TstatList=cpgTstat$tstats,
                    colors=c('blue','red','red','blue'), extend=5000, addPoints=T, 
                    DIR=plotDIR);

# ==========================================================================================================
# ==========================================================================================================
# re-run with different smoothing parameters but same filters for t-stats and dmr filtering
# ==========================================================================================================
# ==========================================================================================================
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# repeat section 3) - prepare for smoothing
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
rm(list=ls()); setwd('~/Documents/_BS-seq_analysis/');
file_cg = 'aligned_trimmed4-98.adapters.q30.m0_bsmap2.9.bam_methratio_samtools0.1.19-m4-CpGcombined.CG';
load(paste0(file_cg, '_bsdAll.RData'));
cpgcommon = .filterLowCoverage2(bsdAll, reqCov=1, fix.seqlevels=T);

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# repeat section 4 - smoothing
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# see section 4) above for how to handle potential .BSmoothAcrossScaffolds hanging
#  may need to do for every run here, otherwise would just loop through ns and h values
# no need to vary maxgap (see 'run_bsmooth_compare_maxGaps.R')
maxgap=10^8; dirPrefix='CpGcombined.CG_4xCov'; skip=TRUE;
.BSmoothAcrossScaffolds(BSD=cpgcommon, NS=25, H=500, MAXGAP=maxgap, skipLowNS=skip, dirPrefix=dirPrefix);
.BSmoothAcrossScaffolds(BSD=cpgcommon, NS=25, H=750, MAXGAP=maxgap, skipLowNS=skip, dirPrefix=dirPrefix);
.BSmoothAcrossScaffolds(BSD=cpgcommon, NS=25, H=1000, MAXGAP=maxgap, skipLowNS=skip, dirPrefix=dirPrefix);
.BSmoothAcrossScaffolds(BSD=cpgcommon, NS=50, H=500, MAXGAP=maxgap, skipLowNS=skip, dirPrefix=dirPrefix);
.BSmoothAcrossScaffolds(BSD=cpgcommon, NS=50, H=750, MAXGAP=maxgap, skipLowNS=skip, dirPrefix=dirPrefix);
.BSmoothAcrossScaffolds(BSD=cpgcommon, NS=50, H=1000, MAXGAP=maxgap, skipLowNS=skip, dirPrefix=dirPrefix);
.BSmoothAcrossScaffolds(BSD=cpgcommon, NS=70, H=500, MAXGAP=maxgap, skipLowNS=skip, dirPrefix=dirPrefix);
.BSmoothAcrossScaffolds(BSD=cpgcommon, NS=70, H=750, MAXGAP=maxgap, skipLowNS=skip, dirPrefix=dirPrefix);

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
# for dmr detection
centile=0.95; dmrMG=300; STAT='tstat.corrected';
# for dmr filtering
fCols=c('n','meanDiff'); fVals=c(3, 0.1); ord='areaStat';

for (ns in nsVec) {
  for (h in hVec) {
# ----------------------------------------------------------------------------------------------------------
    # 5) prepare smoothed data for t-statistic calculation
# ----------------------------------------------------------------------------------------------------------
    # load files output by .BSmoothAcrossScaffolds(), print progress every 200 files
    smoothedDir = paste0(dirPrefix, '_ns', ns, 'h', h, 'mg', maxgap);
    cpgfits = .loadSmoothedScaffolds(smoothedDir, printEvery=200);
    
    # filter to loci with at least 4x coverage in every subject
    cpgfits_f = .filterLowCoverageList(cpgfits, reqCov=reqCov, ns.BSmooth=ns);
    
    # if all smoothed scaffolds had at least ns loci,
    #  get the filtered list from the .filterLowCoverageList output
    # otherwise remove scaffolds with fewer loci than BSmooth ns parameter 
    if (is.null(cpgfits_f$numLociBelowNS)) {
      cpg_forTtest = cpgfits_f$filtered;
    } else {
      toRemove = match(cpgfits_f$numLociBelowNS, names(cpgfits_f$filtered));
      cpg_forTtest = cpgfits_f$filtered[-toRemove];
    }
    
    # save out as RData object
    FILEsm_filt = paste0(smoothedDir, '_fits_', reqCov, 'xCov.RData');
    save(cpgfits, cpgfits_f, cpg_forTtest, file=FILEsm_filt);
    
# ----------------------------------------------------------------------------------------------------------
    # 6) compute t-statistics for each scaffold
# ----------------------------------------------------------------------------------------------------------
    # run t-stat computation
    # will use estimate.var="same" for all runs so no need to check methylation variation
    # other parameters to BSmooth.tstat() that are hard-coded:
    #  local.correct=T, mc.cores=4, verbose=F
    cpgTstat = .BSmooth.tstatList(cpg_forTtest, group1=groupND, group2=groupD, estimate.var='same');
    
    # save out as RData object
    FILEt = gsub('.RData', '', FILEsm_filt);
    FILEt = paste0(FILEt, '_Ttests.RData');
    save(cpgTstat, file=FILEt);
    
# ----------------------------------------------------------------------------------------------------------
    # 7) find differentially methylated regions (DMRs) on each scaffold
# ----------------------------------------------------------------------------------------------------------
    # get DMRs using default settings 
    CUT = .getTstatQuantilesAcrossScaffolds(cpgTstat$tstats, centile=centile);
    dmrsList = .dmrFinderList(cpgTstat$tstats, cutoff=CUT, maxGap=dmrMG, stat=STAT);
    
    # filter DMR list by n (number of loci) and meanDiff as suggested in BSmooth user guide
    dmrsList_f = .filterAndOrderDMRsList(dmrsList$dmrList, filterOn=fCols, filterValues=fVals, orderOn=ord);
    dmrs_f = dmrsList_f$filteredListCollapsed;
    
    # save out as RData object
    filtParams = paste0(paste(fCols, fVals), collapse='_');
    FILEdmrs = gsub('.RData', '', FILEt);
    FILEdmrs = paste0(FILEdmrs, '_DMRs_cut', centile, '_mg', dmrMG, '_', STAT, '_', filtParams,'.RData');
    save(dmrsList_f, dmrs_f, file=FILEdmrs);
    
# ----------------------------------------------------------------------------------------------------------
    # 8) plot DMRs with t-stats
# ----------------------------------------------------------------------------------------------------------
    plotDIR = paste(gsub('.RData', '', FILEdmrs), '_plots', sep='');
    .plotDMRsAcrossList(DMRsList=dmrsList_f$filteredList, BSDList=cpg_forTtest, 
                        BSD.TstatList=cpgTstat$tstats,
                        colors=c('blue','red','red','blue'), extend=5000, addPoints=T, 
                        DIR=plotDIR);
  }
} 

# ==========================================================================================================
# ==========================================================================================================
# identify dmrs that were found in all 9 runs
# ==========================================================================================================
# ==========================================================================================================
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# load dmrs from each run
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# clear workspace, set options and working directory, load functions
rm(list=ls()); options(stringsAsFactors=F); setwd('~/Documents/_BS-seq_analysis/');
source('/Volumes/fishstudies-1/_code/run_bsmooth_functions.R');

# set filename pattern to get all maxGap=1e+08 runs with n3_meanDiff0.1 filtering
pattern = 'mg1e+08_fits_4xCov_Ttests_DMRs_cut0.95_mg300_tstat.corrected_n3_meanDiff0.1.RData';
files = grep(pattern, list.files(), value=T, fixed=T);
files = grep('All|WORKSPACE', files, value=T, invert=T);

# load dmrs from each run as dataframes in a list and save out
dmrf = list();
for (f in files) { print(f); load(f); dmrf[[f]] = dmrs_f }
rm(f, dmrs_f, dmrsList_f);
names(dmrf) = sapply(strsplit(names(dmrf), '_'), function(f) f[3]);  
newprefix = 'CpGcombined.CG_4xCov_ns70_50_25hAll';
#save(dmrf, file=paste0(newprefix, pattern));

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# compute summary statistics for each run
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# 
statcols=7:16;
dmrfstats = list();
for (d in 1:length(dmrf)) { 
  dmrfstats[[names(dmrf)[d]]] = .computeStatsForDMRs(dmrf[[d]], statcols);
}; rm(d);

bigstatmatf = do.call('rbind', lapply(dmrfstats, function(f) f$statmat));
hypratiosf = sapply(dmrfstats, function(f) f$hypRatio);
numdmrsf = sapply(dmrfstats, function(f) f$numDmrs);

stattype = 'mean';
statrows = grep(paste0(stattype, '$'), rownames(bigstatmatf));
tmpmatf = bigstatmatf[statrows, ];
tmpmatf = as.data.frame(cbind(tmpmatf, 
                              numDmrs=sapply(dmrfstats, function(f) f$numDmrs),
                              hypRatio=sapply(dmrfstats, function(f) f$hypRatio)));

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# make plots of summary statistics
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
#
labels = sapply(strsplit(rownames(tmpmatf), 'mg'), function(f) f[1]);
par(mfrow=c(2,7))
for (i in 1:ncol(tmpmatf)) {
  xlim = range(tmpmatf[, i]);
  xrange = xlim[2] - xlim[1];
  xlim[1] = xlim[1] - (.1 * xrange);
  xlim[2] = xlim[2] + (.2 * xrange);
  plot(tmpmatf[, i], 1:nrow(tmpmatf), 
       type='n', main=names(tmpmatf)[i], 
       xlab='', ylab='', xlim=xlim, yaxt='n'); 
  text(tmpmatf[, i], 1:nrow(tmpmatf), labels=labels);
}; rm(i, xlim, xrange);

#
nslab = substr(sapply(strsplit(rownames(tmpmatf), 'h'), function(f) f[1]), 3, 4);
tmphlab = sapply(strsplit(rownames(tmpmatf), 'h'), function(f) f[2]);
tmphlab = gsub(paste0('.',stattype), '', tmphlab, fixed=T);
hlab = sapply(strsplit(tmphlab, 'mg'), function(f) f[1]); rm(tmphlab);

#
par(mfrow=c(4,7));
for (s in 1:ncol(tmpmatf)) {
  WGCNA::verboseBoxplot(tmpmatf[,s], nslab, frame.plot=F, notch=F, 
                        main=names(tmpmatf)[s], xlab='', ylab='', col='grey', border='darkgrey');
  WGCNA::verboseBoxplot(tmpmatf[,s], hlab, frame.plot=F, notch=F, 
                        main=names(tmpmatf)[s], xlab='', ylab='', col='grey', border='darkgrey');
}; rm(s);

#
res = 300;
for (d in 1:length(dmrf)) {
  filebase = paste0('June-CpGcombined.CG_4xCov_', 
                    names(dmrf)[d], 
                    '_fits_4xCov_Ttests_DMRs_cut0.95_mg300_tstat.corrected_n3_meanDiff0.1_postNullFilter');
  .plotDmrStatsHistogramsJpg(dmrf[[d]], res=res, file=paste0(filebase, '_histograms.jpg'), breaks=100);
  .plotDmrStatsHeatmap(dmrf[[d]], 7:15, file=paste0(filebase, '_heatmap.jpg'), res=res);
  jpeg(file=paste0(filebase, '_boxplots_hypoVShyper.jpg'), 
       width=14, height=8.5, units='in', quality=100, type='quartz', res=res);
  .plotDmrStatsAcrossFactor(dmrf[[d]], dmrf[[d]]$direction=='hypo')
  dev.off()
}; rm(d, filebase);

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# get dmrs that were common to all runs
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# get dmrs that were exactly the same across all runs
exdmrs0 = .getCommonDmrsAcrossDmrListExact(dmrf);
exdmrs = exdmrs0$merge;
exstats = .computeStatsForDMRs(exdmrs, statcols);
exints = paste0(exdmrs$chr, ':', exdmrs$start, '-', exdmrs$end);

# get dmrs that had at least 100bp overlap across all runs
ovdmrs = .getCommonDmrsAcrossDmrListOverlap(dmrf, minoverlap=100L);
ovstats = .computeStatsForDMRs(ovdmrs, statcols);
ovints = paste0(ovdmrs$chr, ':', ovdmrs$start, '-', ovdmrs$end);

# check for exact dmr matches that were missed due to having width < minoverlap
mex = setdiff(exints, ovints);
mexints = grep('^scaffold', unlist(strsplit(mex, ':')), value=T, invert=T);
mexwidths = sapply(strsplit(mexints, '-'), function(f) abs(as.numeric(f[1])-as.numeric(f[2])));
summary(mexwidths);
mdmrs = exdmrs[which(!(exints %in% ovints)), ];

#
filebase = gsub('.RData', '', paste0(newprefix, pattern), fixed=T);
save(exdmrs, ovdmrs, mdmrs, file=paste0(filebase, 'exactAnd100bpOverlaps.RData'));

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# make plots of statistics for common dmrs
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# 
dev.off(); res=300; wt=14; ht=8.5; un='in'; qual=100; gdev='quartz'; br=100;

#
fhist = paste0(filebase, '_histograms.jpg');
fheat = paste0(filebase, '_heatmap.jpg');
fbox1 = paste0(filebase, '_boxplots_hypoVShyper.jpg');
fbox2 = paste0(filebase, '_boxplots_widthLessThan120.jpg');

#
.plotDmrStatsHistogramsJpg(exdmrs, res=res, file=fhist, breaks=br);
.plotDmrStatsHeatmap(exdmrs, statcols[-length(statcols)], file=fheat, res=res);
jpeg(file=fbox1, width=wt, height=ht, units=un, quality=qual, type=gdev, res=res);
.plotDmrStatsAcrossFactor(exdmrs, exdmrs$direction=='hypo');
jpeg(file=fbox2, width=wt, height=ht, units=un, quality=qual, type=gdev, res=res);
.plotDmrStatsAcrossFactor(exdmrs, exdmrs$width < 120);
dev.off();

#
filebase = paste0(newprefix, strsplit(pattern,'_')[[1]][1], 
                  '_overlap_100bp', strsplit(pattern,'mg1e+08',fixed=T)[[1]][2]);
fhist = paste0(filebase, '_histograms.jpg');
fheat = paste0(filebase, '_heatmap.jpg');
fbox1 = paste0(filebase, '_boxplots_hypoVShyper.jpg');

#
.plotDmrStatsHistogramsJpg(ovdmrs, res=res, file=fhist, breaks=br);
.plotDmrStatsHeatmap(ovdmrs, statcols[-length(statcols)], file=fheat, res=res);
jpeg(file=fbox1, width=wt, height=ht, units=un, quality=qual, type=gdev, res=res);
.plotDmrStatsAcrossFactor(ovdmrs, ovdmrs$direction=='hypo');
dev.off();

#
.plotDmrStatsHistogramsJpg(dmrs, res=res, file=fhist, breaks=br);
.plotDmrStatsHeatmap(dmrs, statcols[-length(statcols)], file=fheat, res=res);
jpeg(file=fbox1, width=wt, height=ht, units=un, quality=qual, type=gdev, res=res);
.plotDmrStatsAcrossFactor(dmrs, dmrs$direction=='hypo');
dev.off();

# ==========================================================================================================
# ==========================================================================================================
# assign genomic features to common dmrs
# ==========================================================================================================
# ==========================================================================================================
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# load common dmrs and keep only those that passed hand curation
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# clear workspace, set options and working directory, load functions
rm(list=ls()); options(stringsAsFactors=F); setwd('~/Documents/_BS-seq_analysis/');
source('/Volumes/fishstudies-1/_code/run_bsmooth_functions.R');
library('GenomicRanges');

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

# get ids of dmrs that passed curation
goodstatus = c('ok','ok_s','ok_d');
goodhc2ints = hc2ints[hc2$status %in% goodstatus];

dmrs = dmrs0[dmrs0ints %in% goodhc2ints, ];
dmrsints = paste0(dmrs$chr, ':', dmrs$start, '-', dmrs$end);

# ----------------------------------------------------------------------------------------------------------
# if (!all(dmrsints %in% hc2ints)) {
#   missingdmrsints = dmrsints[!(dmrsints %in% hc2ints)];
#   missingdmrs = dmrs[dmrsints %in% missingdmrsints, ];
#   missingdmrsGR = GRanges(seqnames=missingdmrs$chr, 
#                           ranges=IRanges(start=missingdmrs$start, end=missingdmrs$end));
#   missingdmrsOv = do.call('rbind', .getOverlapsWithAmountsUnstranded(missingdmrsGR, hc2GR));
#   # scaffold_302:422660-423113, scaffold_60:1836444-1836766, 
#   #  scaffold_715:117539-117942, scaffold_567:200482-200697
#   #missingdmrsOv$mcols[is.na(missingdmrsOv$mcols)] = c('','','','')
#   goodhc2ints = c(goodhc2ints, rownames(missingdmrsOv)[missingdmrsOv$mcols %in% goodstatus]);
# }

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# load burtoni genome information and features
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# load a list with the following elements:
#  $gffGenes - gene gff dataframe
#  $gffGenesGR - gene gff GRanges
#  $gff - full gff in dataframe
#  $recipBL - dataframe mapping burtoni genes to other fish and human ids
#  $gois - vector of ids for typical Fernald lab genes of interest
#  $lookup - dataframe mapping burtoni transcript and protein ids to genes
an = .loadBurtoniGenesApr2016();

# load a list containing other types of genomic features
an2 = .loadBurtoniOtherFeatures('/Volumes/fishstudies-1/_Burtoni_annotations/');

# load a dataframe with mappings of burtoni genes to human ids that aren't in an$recipBL
handannos = read.table('~/Documents/_BS-seq_analysis/dmr_handAnnotations_hsaEntrez.txt', sep='\t', header=F, quote='');
handannos = handannos[!duplicated(handannos[,2]), ];
handannos = handannos[grep('^[0-9]', handannos[,1]), ];
names(handannos) = c('hsaHomologEntrez','gene','description');

# combine reciprocal BLAST results and hand annotations for mapping burtoni genes to human entrez ids
annoCombo = .combineHandAndMachineAnnos(an$recipBL, handannos);
abHsMap = strsplit(annoCombo$new$hsaHomologEntrez, ',');
names(abHsMap) = annoCombo$new$gene;
abHsMap = abHsMap[!sapply(abHsMap, function(f) all(is.na(f)))];

# load scaffold lengths into a named vector
sl = read.table('/Volumes/fishstudies-1/_Burtoni_genome_files/scaffold_lengths'); 
sl0 = sl; sl = sl0[, 2]; names(sl) = sl0[, 1];

# build table of statistics about genes on each scaffold
slGeneStats = as.data.frame(matrix(nrow=length(sl), ncol=7));
statnames = c('length','numGenes','avgGeneWidth','avgGeneGC','avgNumIsoforms','pctBpGene','genesPerMb');
dimnames(slGeneStats) = list(names(sl), statnames);
slGeneStats$length = sl;

# count genes on each scaffold
tmpgenenum = table(an$gffGenesDF$seqnames);
tmprows = match(rownames(slGeneStats), names(tmpgenenum));
slGeneStats$numGenes[which(!is.na(tmprows))] = tmpgenenum[na.omit(tmprows)];
rm(tmpgenenum, tmprows);
slGeneStats$numGenes[which(is.na(slGeneStats$numGenes))] = 0;

# get other statistics about genes by scaffold
for (i in which(slGeneStats$numGenes != 0)) {
  sc = rownames(slGeneStats)[i];
  if (i %% 1000 == 0) {cat(sc,'...')}
  scrows = an$gffGenesDF$seqnames == sc;
  slGeneStats$avgGeneWidth[i] = mean(an$gffGenesDF$width[scrows]);
  slGeneStats$avgGeneGC[i] = signif(mean(an$gffGenesDF$GC[scrows]), 4);
  slGeneStats$avgNumIsoforms[i] = signif(mean(an$gffGenesDF$isoforms[scrows]), 4);
  slGeneStats$pctBpGene[i] = signif(slGeneStats$avgGeneWidth[i] / slGeneStats$length[i], 4);
  slGeneStats$genesPerMb[i] = slGeneStats$numGenes[i] / slGeneStats$length[i] * 1000000
}; rm(i,sc,scrows);

#
.verboseScatterplotAllColumnPairs(slGeneStats, mfrow=c(3,7), abline.col='red',
                                  bg='grey', col='black', pch=21, cex=1.8);

#
.verboseScatterplotAllColumnPairs(an$gffGenesDF[,c(4,9,10)], mfrow=c(1,3), abline.col='red',
                                  bg=WGCNA::labels2colors(an$gffGenesDF$strand), col='black', 
                                  pch=21, cex=1.8);
.verboseScatterplotAllColumnPairs(an$gffGenesDF[,c(4,9,10)], mfrow=c(1,3), abline.col='red',
                                  bg=WGCNA::labels2colors(an$gffGenesDF$strand), col='black', 
                                  pch=21, cex=1.8,
                                  corOptions='method="s"');


# load dataframe of non-coding regions conserved from drerio
zcne = read.table('~/Documents/_annotationsDec2015/hBurtoni_zCNE_alignments.out.bed', sep='\t');

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# get/make GRanges objects for dmrs and genome features of interest
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# dmrs
mcolCols = which(names(dmrs)=='n'):ncol(dmrs);
dmrsGR = GRanges(seqnames=dmrs$chr, ranges=IRanges(start=dmrs$start, end=dmrs$end), mcols=dmrs[,mcolCols]);
dmrsGR20kb = .getSymmetricWindowsAroundFeatures(gr=dmrsGR, bpToAdd=20000, chrLengthVec=sl, sortOutput=F);
mcols(dmrsGR20kb) <- as.data.frame(cbind(as.data.frame(mcols(dmrsGR20kb)), dmrsints));

# gene bodies, transcription start sites, basal regulatory regions
geneGR = an$gffGenesGR;
tssGR = .getFeatureStartsFromGR(geneGR, stranded=T)$gr;
basalReg = .getWindowsAroundGeneStarts(gr=geneGR, up=5000, down=1000, sortOutput=F);
basalReg = .fixInvalidPositions(basalReg$gr, chrLengthVec=sl);
brGR = basalReg$gr;

# 3' basal regulatory region
.getWindowsAroundGeneEnds = function (gr, down=5000, up=1000, sortOutput=FALSE) {
  ofwd = which(strand(gr)=='+');
  orev = which(strand(gr)=='-');
  strand(gr)[ofwd] <- '-';
  strand(gr)[orev] <- '+';
  newgr = .getWindowsAroundGeneStarts(gr=gr, up=down, down=up, sortOutput=sortOutput)$gr;
  strand(newgr)[ofwd] <- '+';
  strand(newgr)[orev] <- '-';
  return(list(gr=newgr, df=as.data.frame(newgr)));
}

br3primeGR = .getWindowsAroundGeneEnds(gr=geneGR, down=5000, up=1000, sortOutput=F);
br3primeGR = .fixInvalidPositions(br3primeGR$gr, chrLengthVec=sl)$gr;

# .buildFeatureGR doesn't work as expected!
# # exons, coding sequence, mRNA (### utrs??)
# exonGR = .buildFeatureGR(an$gff, feature='exon');
# cdsGR = .buildFeatureGR(an$gff, feature='CDS');
# mrnaGR = .buildFeatureGR(an$gff, feature='mRNA');
# 
# # mRNA start sites to help find dmrs over specific transcripts
# mrnaSsGR = .getFeatureStartsFromGR(mrnaGR, stranded=T)$gr;

# make exon GR
gffExon0 = subset(an$gff, V3=='exon');
#gffExon = gffExon0[!duplicated(paste0(gffExon0$V1, gffExon0$V4, gffExon0$V5)), ];
gffExon = gffExon0
exonGR = GRanges(seqnames=gffExon$V1, 
                 ranges=IRanges(start=gffExon$V4, end=gffExon$V5), 
                 strand=gffExon$V7, 
                 mcols=as.data.frame(cbind(.parseGffMetaCol(gffExon, pattern='gene='),
                                           gffExon$V11,
                                           gffExon$V18)));

# make cds GR
gffCDS0 = subset(an$gff, V3=='CDS');
#gffCDS = gffCDS0[!duplicated(paste0(gffCDS0$V1, gffCDS0$V4, gffCDS0$V5)), ];
gffCDS = gffCDS0;
cdsGR = GRanges(seqnames=gffCDS$V1, 
                 ranges=IRanges(start=gffCDS$V4, end=gffCDS$V5), 
                 strand=gffCDS$V7, 
                 mcols=as.data.frame(cbind(.parseGffMetaCol(gffCDS, pattern='gene='),
                                           gffCDS$V11,
                                           gffCDS$V18)));

# make mRNA GR 
gffmrna0 = subset(an$gff, V3=='mRNA');
#gffmrna = gffmrna0[!duplicated(paste0(gffmrna0$V1, gffmrna0$V4, gffmrna0$V5)), ];
gffmrna = gffmrna0;
mrnaGR = GRanges(seqnames=gffmrna$V1, 
                ranges=IRanges(start=gffmrna$V4, end=gffmrna$V5), 
                strand=gffmrna$V7, 
                mcols=as.data.frame(cbind(.parseGffMetaCol(gffmrna, pattern='gene='),
                                          gffmrna$V11,
                                          gffmrna$V18)));

# make intron GR
load('~/Documents/_annotationsDec2015/intron.gff.RData');
x.int_rnanum = .parseGffMetaCol(x.int, pattern='Parent=');
gffIntron = as.data.frame(cbind(x.int, gene=an$lookup$gene[match(x.int_rnanum, an$lookup$rna_num)]));
intronGR = GRanges(seqnames=gffIntron$V1, 
                 ranges=IRanges(start=gffIntron$V4, end=gffIntron$V5), 
                 strand=gffIntron$V7, 
                 mcols=as.data.frame(cbind(geneSym=gffIntron$gene,
                                           gffIntron$V9)));

# transposons, microRNAs, conserved non-coding regions
teGR = an2$gr$te;
miGR = an2$gr$mi;
zranges = IRanges(start=zcne[,2], end=zcne[,3]);
zGR = GRanges(seqnames=zcne[,1], ranges=zranges, strand=zcne[,6], mcols=zcne[,4:5]);


gene_z_OV = .getOverlapsWithAmountsUnstranded(zGR, geneGR)
gene_z_OVhits = gene_z_OV[!is.na(gene_z_OV)]
gene_z_OVhitsGenes = lapply(gene_z_OVhits, function(f) f$mcols.geneSym)
gene_z_OVhitsGenes_numz = table(unlist(gene_z_OVhitsGenes));


br_te_OV = as.data.frame(findOverlaps(brGR,teGR))
br_with_te = mcols(brGR)[unique(br_te_OV$queryHits), 1];

gene_te_OV = as.data.frame(findOverlaps(geneGR,teGR))
gene_with_te = mcols(geneGR)[unique(gene_te_OV$queryHits), 1];

# gene_gene_OVnested = as.data.frame(findOverlaps(geneGR,ignoreSelf=T));
gene_gene_OV = as.data.frame(findOverlaps(geneGR,ignoreSelf=T,ignore.strand=T));
gene_gene_OVlist = apply(gene_gene_OV, 1, function(f) geneGR[unlist(f)]);
gene_gene_OVlist_biotypes = lapply(gene_gene_OVlist, function(f) mcols(f)[,3]);


# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# get dmr overlaps and nearest neighbors for all features
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
####
.getFeatureOverlapsAcrossGroupsNew = function (qGR, featureGRlist, groupVec, outNames=NULL, ...) {
  qname = deparse(substitute(qGR));
  groupVec = as.character(groupVec);
  groupnames = names(table(groupVec));
  for (feat in 1:length(featureGRlist)) {
    thisfeat = names(featureGRlist)[feat];
    fGR = featureGRlist[[feat]];
    cat(paste0('\n------------------------------------ ',thisfeat,'...\n'));
    thisname = paste0(qname,'_',thisfeat,'OV');
#     cat(paste0('\n----------------------------------------- ',thisname,'...\n'));
    tmp = .getOverlaps(qGR, fGR, ...);
    if (is.character(outNames) & length(outNames)==length(tmp)) {
      names(tmp) = outNames;
    }
    assign(thisname, tmp, pos=1);
    for (gr in groupnames) {
      assign(paste0(thisname,'.',gr), tmp[groupVec==gr], pos=1);
    }
  }
}

featureGRs = list(gene=geneGR, te=teGR, z=zGR,
                  br=brGR, br3prime=br3primeGR, 
                  exon=exonGR, cds=cdsGR,intron=intronGR);

.getFeatureOverlapsAcrossGroupsNew(dmrsGR, featureGRs, dmrs$direction);
.getFeatureOverlapsAcrossGroupsNew(dmrsGR20kb, featureGRs, dmrs$direction, 
                                   outNames=as.data.frame(mcols(dmrsGR20kb))$dmrsints);



dmrs_snps_OV = as.data.frame(findOverlaps(dmrsGR, 
                                          .buildGR_for_intervals(paste0(genomesnps$V1,':',genomesnps$V3,'-',genomesnps$V3))))
####





# loop through different feature types to create lists of nearest neighbors and overlaps
features = c('gene', 'tss', 'br', 'br3prime', 'mrna', 'exon', 'cds', 'intron', 'te', 'mi', 'z');
groups = list(hypo=dmrsints[dmrs$direction == 'hypo'], 
              hyper=dmrsints[dmrs$direction == 'hyper']);




#
.getDmrFeatureOverlapsAcrossGroups(dmrsGR, features, groups);
.getDmrFeatureOverlapsAcrossGroupsNoNN(dmrsGR20kb, features, groups);

.combineOvAndNnByDmr = function (features) {
  ov = paste0(features,'OV');
  nn = paste0(features,'NN');
  features = intersect(c(ov,nn), ls(pos=1));
  dmrints = names(get(features[1]));
  cat('Checking that all objects contain exact same dmrs...\n    ');
  for (feat in features) { 
    cat(feat,'...'); 
    if (any(names(get(feat)) != dmrints)) { stop() } 
  }
  dmrsGRhits = list();
  for (dmrint in dmrsints) {
    dmrsGRhits[[dmrint]] = list(nn=list(),ov=list());
    for (feat in grep('OV$',features,value=T)) {
      f = gsub('OV','',feat);
      dmrsGRhits[[dmrint]]$ov[[f]] = get(feat)[[dmrint]];
    }
    for (feat in grep('NN$',features,value=T)) {
      f = gsub('NN','',feat);
      dmrsGRhits[[dmrint]]$nn[[f]] = get(feat)[[dmrint]];
    }
  }
  return(dmrsGRhits);
}

.combineOvByDmr = function (features) {
  ov = paste0(features,'OV');
  features = intersect(c(ov), ls(pos=1));
  dmrints = names(get(features[1]));#print(head(dmrints))
  cat('Checking that all objects contain exact same dmrs...\n    ');
  for (feat in features) { 
    cat(feat,'...'); 
    if (any(names(get(feat)) != dmrints)) { stop() } 
  }
  dmrsGRhits = list();
  for (i in 1:length(dmrsints)) {#print('------------------------------------------------------')
    dmrsGRhits[[dmrints[i]]] = list();#print(dmrsGRhits[dmrint])
    for (feat in grep('OV$',features,value=T)) {
      f = gsub('OV','',feat);#print(f)
     
      tmp = get(feat);
      dmrsGRhits[[dmrints[i]]][[f]] = tmp[[i]];
      
    }
  }
  return(dmrsGRhits);
}


# 
dmrsGRhits = .combineOvAndNnByDmr(paste0('dmrsGR_',features));
dmrsGRhitsOV = lapply(dmrsGRhits, function(f) f$ov);
dmrsGRhitsNN = lapply(dmrsGRhits, function(f) f$nn);

dmrsGRhitsOVgenes = list()
for (i in 1:length(dmrsGRhitsOV)) {
  this = dmrsGRhitsOV[[i]];
  this = this[names(this) %in% c('dmrsGR_gene','dmrsGR_br','dmrsGR_br3prime')]; 
  this = this[!is.na(this)];
  dmrsGRhitsOVgenes[[names(dmrsGRhitsOV)[i]]] = unique(c(this$dmrsGR_gene$mcols.geneSym,
                                                         this$dmrsGR_br$mcols.geneSym,
                                                         this$dmrsGR_br3prime$mcols.geneSym));
}
rm(i,this)

.make_numOV = function (dmrsGRhitsOV, features, prefix='dmrsGR_') {
  numOV = as.data.frame(matrix(nrow=nrow(dmrs), ncol=length(features)));
  dimnames(numOV) = list(dmrsints, features);
  for (feat in features) {
    featOV = lapply(dmrsGRhitsOV, function(f) f[[which(names(f) == paste0(prefix,feat))]] );
    featOVnums = sapply(featOV, nrow);
    featOVnums[which(sapply(featOVnums, is.null))] = 0;
    featOVnums = unlist(featOVnums);
    numOV[, match(feat, names(numOV))] = featOVnums;
  }
  return(numOV)
}

# get num overlap hits for each feature type
numOV = .make_numOV(dmrsGRhitsOV, features);
numOV$br.only = rep(0,nrow(numOV));
numOV$br.same = rep(0,nrow(numOV));
numOV$br.diff = rep(0,nrow(numOV));
numOV$br3prime.only = rep(0,nrow(numOV));
numOV$br3prime.same = rep(0,nrow(numOV));
numOV$br3prime.diff = rep(0,nrow(numOV));
numOV$fc.anti = rep(0,nrow(numOV));  

brDMRsCheckGene = .checkOVmatchGenes(dmrsGRhitsOV, numOV, 'br');
br3primeDMRsCheckGene = .checkOVmatchGenes(dmrsGRhitsOV, numOV, 'br3prime');

numOV$br.only[rownames(numOV) %in% brDMRsCheckGene$only] = 1;
numOV$br.same[rownames(numOV) %in% brDMRsCheckGene$same] = 1;
numOV$br.diff[rownames(numOV) %in% brDMRsCheckGene$diff] = 1;
numOV$br3prime.only[rownames(numOV) %in% br3primeDMRsCheckGene$only] = 1;
numOV$br3prime.same[rownames(numOV) %in% br3primeDMRsCheckGene$same] = 1;
numOV$br3prime.diff[rownames(numOV) %in% br3primeDMRsCheckGene$diff] = 1;

for (i in 1:nrow(numOV)) {
  if (numOV$gene[i] > 0) {
    thisdmr = rownames(numOV)[i];
    thisdmrgene = dmrsGRhitsOVgenes[[thisdmr]];
    samecheck = any(thisdmrgene %in% rownames(allDmrGenesFcSame));
    anticheck = any(thisdmrgene %in% rownames(allDmrGenesFcAnti));
    if (samecheck & anticheck) {
      numOV$fc.anti[i] = 0;
    } else if (samecheck) {
      numOV$fc.anti[i] = -1;
    } else if (anticheck) {
      numOV$fc.anti[i] = 1;
    } else {
      numOV$fc.anti[i] = 0;
    }
  }
}
rm(i,thisdmr,thisdmrgene,samecheck,anticheck)

numOVyn = as.data.frame(apply(apply(numOV, 2, as.logical), 2, as.numeric));
rownames(numOVyn) = rownames(numOV);

numOVyntmp = numOVyn[, -match(c('tss','exon','mi','gene','mrna','br','br3prime'), names(numOVyn))];

numOVfeat=numOVyntmp
for(i in 1:ncol(numOVfeat)) { numOVfeat[(numOVfeat[,i] > 0) & !is.na(numOVfeat[,i]), i] = names(numOVfeat)[i] }; rm(i);
#numOVfeat$hypo[numOVfeat$hypo=='0'] = 'hyper';
fg = apply(numOVfeat, 1, function(f) paste0(f[f!='0' & !is.na(f)], collapse=':'))
fg[fg==''] = 'none';


# numOV2 = numOV;
# mrnarows = numOV2$mrna > 0;
# exoncdsintroncols = which(names(numOV2) %in% c('exon','cds','intron'))
# numOV2[mrnarows, exoncdsintroncols] = numOV2[mrnarows, exoncdsintroncols] / numOV2$mrna[mrnarows];
# fthresh = nrow(numOV2) * .05;
# numOV2 = numOV2[, apply(numOV2,2,sum) > fthresh];
# numOV2 = numOV2[, -which(names(numOV2) %in% c('mrna','exon'))];
# numOV2 = as.data.frame(cbind(numOV2 
#                              ,hypo=hyphyp
#                            #  ,w125=as.numeric(dmrs$width>125)
#                              ));
# 
# numOV2yn = numOV2;
# for(i in 1:ncol(numOV2yn)) { numOV2yn[numOV2yn[,i] > 0, i] = 1 }; rm(i);
# 
# 
# #thisnumOV = numOV2[,-1]
# thisnumOV = as.data.frame(cbind(numOV2yn[,-1], 
#                                 cds.intron.diff=numOV2yn$intron+numOV2yn$cds,
#                                 gene.br.diff=numOV2yn$br+numOV2yn$gene,
#                                 gene.br3prime.diff=numOV2yn$br3prime+numOV2yn$gene,
#                                 gene.te.diff=numOV2yn$te+numOV2yn$gene,
#                                 gene.z.diff=numOV2yn$z+numOV2yn$gene,
#                                 hypo=as.numeric(gsub('hyper','0',gsub('hypo','1',dmrs$direction)))))

fgroupings = names(table(fg));
numOVpca = as.data.frame(matrix(0,nrow=nrow(numOV),ncol=length(fgroupings)))
dimnames(numOVpca) = list(rownames(numOV), fgroupings);
for (thisfg in fgroupings) {
  numOVpca[[thisfg]][fg==thisfg] = 1;
}; rm(thisfg)



fnames = c('cds.only','intron.only','cds.intron',
           'br.same','br.diff','br.only',
           'br3prime.same','br3prime.diff','br3prime.only',
           'te.only','te.br','te.br3prime','te.cds','te.intron',
           'z.only','z.br','z.br3prime','z.te','z.intron');



numOVpca = as.data.frame(matrix(rep(0, nrow(numOVyntmp)*length(fnames)),
                                nrow=nrow(numOVyntmp),ncol=length(fnames)));
dimnames(numOVpca) = list(rownames(numOVyntmp), fnames);

#numhits = apply(numOVyntmp[,-which(names(numOVyntmp)=='fc.anti')], 1, sum);
numhits = apply(numOVyntmp, 1, sum);

numOVpca$cds.only[numhits==1 & numOVyntmp$cds==1] = 1;
numOVpca$intron.only[numhits==1 & numOVyntmp$intron==1] = 1;
numOVpca$cds.intron[numhits>1 & numOVyntmp$intron==1 & numOVyntmp$cds==1] = 1;

numOVpca$br.same = numOVyntmp$br.same;
numOVpca$br.diff = numOVyntmp$br.diff;
numOVpca$br.only = numOVyntmp$br.only;

numOVpca$br3prime.same = numOVyntmp$br3prime.same;
numOVpca$br3prime.diff = numOVyntmp$br3prime.diff;
numOVpca$br3prime.only = numOVyntmp$br3prime.only;

brcheck = numOVyntmp$br.only==1 | numOVyntmp$br.same==1 #| numOVyntmp$br.diff==1;
br3primecheck = numOVyntmp$br3prime.only==1 | numOVyntmp$br3prime.same==1 #| numOVyntmp$br3prime.diff==1;

numOVpca$te.only[numhits==1 & numOVyntmp$te==1] = 1;
numOVpca$te.br[numhits>1 & numOVyntmp$te==1 & brcheck] = 1;
numOVpca$te.br3prime[numhits>1 & numOVyntmp$te==1 & br3primecheck] = 1;
numOVpca$te.cds[numhits>1 & numOVyntmp$te==1 & numOVyntmp$cds==1] = 1;
numOVpca$te.intron[numhits>1 & numOVyntmp$te==1 & numOVyntmp$intron==1] = 1;

numOVpca$z.only[numhits==1 & numOVyntmp$z==1] = 1;
numOVpca$z.br[numhits>1 & numOVyntmp$z==1 & brcheck] = 1;
numOVpca$z.br3prime[numhits>1 & numOVyntmp$z==1 & br3primecheck] = 1;
numOVpca$z.te[numhits>1 & numOVyntmp$z==1 & numOVyntmp$te==1] = 1;
numOVpca$z.intron[numhits>1 & numOVyntmp$z==1 & numOVyntmp$intron==1] = 1;





thisnumOV = numOVpca
dmrpca = prcomp(thisnumOV, scale.=T, center=T); summary(dmrpca); 
dmrpcacor=cor(thisnumOV, dmrpca$x); 
dmrpcacorpval=WGCNA::corPvalueFisher(dmrpcacor,ncol(thisnumOV));
plot(dmrpca$x[,1], dmrpca$x[,2], bg=WGCNA::numbers2colors(thisnumOV$intron), pch=21);

intron_genes = unique(unlist(dmrsGR_geneOVhitsGenesAll[names(dmrsGR_geneOVhitsGenesAll) %in% names(fg)[fg=='intron']]))



# how many dmrs overlapped a gene?
dmrsGR_geneOVhits = dmrsGR_geneOV[numOV$gene > 0];
length(dmrsGR_geneOVhits) / length(dmrsGR_geneOV);

dmrsGR_brOVhits = dmrsGR_brOV[numOV$br > 0];
length(dmrsGR_brOVhits) / length(dmrsGR_brOV);

dmrsGR_br3primeOVhits = dmrsGR_br3primeOV[numOV$br3prime > 0];
length(dmrsGR_br3primeOVhits) / length(dmrsGR_br3primeOV);

# what were the genes?
dmrsGR_geneOVhitsGenesAll = lapply(dmrsGR_geneOVhits, function(f) f$mcols.geneSym);
dmrsGR_brOVhitsGenesAll = lapply(dmrsGR_brOVhits, function(f) f$mcols.geneSym);
dmrsGR_br3primeOVhitsGenesAll = lapply(dmrsGR_br3primeOVhits, function(f) f$mcols.geneSym);
dmrsGR_OVhitsGenesAll = c(dmrsGR_geneOVhitsGenesAll, dmrsGR_brOVhitsGenesAll, dmrsGR_br3primeOVhitsGenesAll)

allDmrGenes = unique(c(unlist(dmrsGR_geneOVhitsGenesAll), 
                       unlist(dmrsGR_brOVhitsGenesAll), 
                       unlist(dmrsGR_br3primeOVhitsGenesAll)));

allDmrGenes.hypo = unique(c(unlist(dmrsGR_geneOVhitsGenesAll[names(dmrsGR_geneOVhitsGenesAll) %in% groups$hypo]), 
                       unlist(dmrsGR_brOVhitsGenesAll[names(dmrsGR_brOVhitsGenesAll) %in% groups$hypo]), 
                       unlist(dmrsGR_br3primeOVhitsGenesAll[names(dmrsGR_br3primeOVhitsGenesAll) %in% groups$hypo])));

allDmrGenes.hyper = unique(c(unlist(dmrsGR_geneOVhitsGenesAll[names(dmrsGR_geneOVhitsGenesAll) %in% groups$hyper]), 
                            unlist(dmrsGR_brOVhitsGenesAll[names(dmrsGR_brOVhitsGenesAll) %in% groups$hyper]), 
                            unlist(dmrsGR_br3primeOVhitsGenesAll[names(dmrsGR_br3primeOVhitsGenesAll) %in% groups$hyper])));

dmrsGR_geneOVhitsGenes = unique(unlist(dmrsGR_geneOVhitsGenesAll));
gffGenesDF_dmrs = subset(an$gffGenesDF, geneSym %in% dmrsGR_geneOVhitsGenes);
gffGenesDF_dmrs = gffGenesDF_dmrs[match(dmrsGR_geneOVhitsGenes, gffGenesDF_dmrs$geneSym), ];

# what parts of the genes did they hit?
tmp = subset(numOV, gene>0);
tmpdmrs = apply(tmp, 2, function(f) rownames(tmp)[f>0]);
tmpdmrs$intron = setdiff(tmpdmrs$gene, tmpdmrs$exon);

dev.off()
xxx=.draw_quad_venn_from_groups(tmpdmrs$intron, tmpdmrs$brsame, setdiff(tmpdmrs$exon, tmpdmrs$cds), tmpdmrs$cds, 
                                fill=c('grey60','green','orangered','purple'), col='grey',
                                groupNames=c('intron','br.same','utr.only','cds'));

dev.off()
xxx=.draw_quad_venn_from_groups(tmpdmrs$intron, tmpdmrs$brsame, rownames(subset(numOV, exon!=cds)), tmpdmrs$cds, 
                                fill=c('grey60','green','orangered','purple'), col='grey',
                                groupNames=c('intron','br.same','utr','cds'));


# assumes gene name is always first column in mcols
.checkOVmatchGenes = function (dmrsGRhitsOV, numOV, feature='br') {
 # listfeatures = sapply(strsplit(names(dmrsGRhitsOV[[1]]),'_'), function(f) f[2]);
 # if (all(listfeatures))
  numOVgeneCol = which(names(numOV)=='gene');
  numOVfeatureCol = which(names(numOV)==feature);
  multirows = which(numOV[,numOVfeatureCol] > 1);
  numOVmulti = numOV[multirows, ];
  numOV = numOV[-multirows, ];
  fdmrs = rownames(numOV)[numOV[,numOVgeneCol]>0 & numOV[,numOVfeatureCol]>0];#print(fdmrs)
  fsamegene = c();
  for (i in 1:length(fdmrs)) {
    #print('-------------')
    this = dmrsGRhitsOV[[fdmrs[i]]];
    this = this[which(names(this) %in% c('dmrsGR_gene',paste0('dmrsGR_',feature)))];#print(this)
    gnames = sapply(this, function(f) f[,6]);       # assumes gene name is always first column in mcols
    if (length(unique(gnames))==1) {
      fsamegene = c(fsamegene, i);
    }
  }
  res = list(same=fdmrs[fsamegene], 
             diff=fdmrs[-fsamegene],
             only=rownames(numOV)[numOV[,numOVgeneCol]==0 & numOV[,numOVfeatureCol]>0]);
  for (i in 1:nrow(numOVmulti)) {
    if (numOVmulti[i, numOVgeneCol] == 0) {
      res$only = c(res$only, rownames(numOVmulti)[i]);
    } else if (numOVmulti[i, numOVgeneCol] %in% c(1,2)) {
      res$same = c(res$same, rownames(numOVmulti)[i]);
      res$diff = c(res$diff, rownames(numOVmulti)[i]);
    } else {
      res$same = c(res$same, rownames(numOVmulti)[i]);
      res$diff = c(res$diff, rownames(numOVmulti)[i]);
      warning(paste0('check ',rownames(numOVmulti)[i], ', hit >2 genes'));
    }
  }
  return(res);
}


# what other features did they hit?


# brDMRs = rownames(subset(numOV, br>0 & gene>0));
# brDMRsSameGene = c();
# for (i in 1:length(brDMRs)) {
#   print('-------------')
#   this = dmrsGRhitsOV[[brDMRs[i]]];
#   this = this[which(names(this) %in% c('dmrsGR_gene','dmrsGR_br'))];
#   print(this)
#   gnames = sapply(this, function(f) f$mcols.geneSym);
#   if (length(unique(gnames)) == 1) {
#     brDMRsSameGene = c(brDMRsSameGene, i);
#   }
# }
# rm(i,this,gnames);
# 
# tmpdmrs$brsame = brDMRs[brDMRsSameGene];
# tmpdmrs$brdiff = unique(c(brDMRs[-brDMRsSameGene], rownames(subset(numOV, br>0 & gene==0))))

dev.off()
xxx=.draw_quad_venn_from_groups(rownames(subset(numOV,gene>0)), rownames(subset(numOV, te>0)),
                                tmpdmrs$brsame, tmpdmrs$brdiff,
                                fill=c('grey','blue','green','yellow'), col='grey', 
                                groupNames=c('gene','te','br.same','br.diff'));

# how many dmrs hit each gene?
# how many genes were overlapped by >1 dmr?
# dmrsGR_geneOVhitsGenes_numDMRs = table(unlist(lapply(dmrsGR_geneOVhits, function(f) f$mcols.geneSym)));
# dmrsGR_geneOVhitsGenes_numDMRs = dmrsGR_geneOVhitsGenes_numDMRs[match(dmrsGR_geneOVhitsGenes, names(dmrsGR_geneOVhitsGenes_numDMRs))]
# dmrsGR_geneOVhitsGenes_numDMRsMulti = dmrsGR_geneOVhitsGenes_numDMRs[dmrsGR_geneOVhitsGenes_numDMRs > 1]
tmp = do.call('rbind',dmrsGR_geneOVhits)
dmrsByGene = split(tmp, tmp$mcols.geneSym); rm(tmp);
dmrsByGeneMulti = dmrsByGene[sapply(dmrsByGene, nrow) > 1];

# build table of multi-dmr genes
gffgenesyms = .parseGffMetaCol(subset(an$gff, V3=='mRNA'), pattern='gene=');
gffgenedesc = .parseGffMetaCol(subset(an$gff, V3=='mRNA'), pattern='product=');

tmplist = list()
for (g in 1:length(dmrsByGeneMulti)) {
  this = dmrsByGeneMulti[[g]];
  gdmrs = sapply(strsplit(rownames(this),'.',fixed=T), function(f) f[1]);
  tmp = as.data.frame(cbind(DMR=gdmrs, 
                            gene=this$mcols.geneSym, 
                            biotype=this$mcols.biotype,
                            desc=gffgenedesc[match(unique(this$mcols.geneSym), gffgenesyms)],
                            hypermethylated=dmrs[match(gdmrs, dmrsints), ]$direction,
                            DMR.fc=signif(abs(dmrFc[match(gdmrs, names(dmrFc))]), 3)));
  tmplist[[g]] = tmp;
}; rm(g,this,tmp);
names(tmplist) = names(dmrsByGeneMulti)
tmplistfcmeans = sapply(tmplist, function(f) mean(as.numeric(f$DMR.fc)))
tmplist2 = tmplist[order(-abs(tmplistfcmeans))];
tmplistfcmeans = sapply(tmplist2, function(f) mean(as.numeric(f$DMR.fc)))
tmp = as.data.frame(cbind(gene=sapply(tmplist2, function(f) unique(f$gene)),
                          biotype=sapply(tmplist2, function(f) unique(f$biotype)), 
                          desc=sapply(tmplist2, function(f) unique(f$desc)),
                          hypermethylated.DMRs=sapply(sapply(tmplist2, function(f) paste0(names(table(f$hypermethylated)), ':', table(f$hypermethylated))), function(f) paste0(f,collapse=','))
                          ));
tmp$hypermethylated.DMRs = gsub('hypo','D',tmp$hypermethylated.DMRs);
tmp$hypermethylated.DMRs = gsub('hyper','ND',tmp$hypermethylated.DMRs);
tmp$desc = gsub('%2C','',tmp$desc);
tmp$desc = sapply(strsplit(tmp$desc, ' variant'), function(f) f[1]);
tmp$desc = sapply(strsplit(tmp$desc, ' transcript'), function(f) f[1]);

dmrsByGeneMultiSummary = tmp; rm(tmp)
dmrsByGeneMultiSummaryFull = do.call('rbind', tmplist2); rm(tmplist2)
dmrsByGeneMultiSummaryFull$hypermethylated = gsub('hypo','D',dmrsByGeneMultiSummaryFull$hypermethylated)
dmrsByGeneMultiSummaryFull$hypermethylated = gsub('hyper','ND',dmrsByGeneMultiSummaryFull$hypermethylated)


WGCNA::verboseBoxplot(gffGenesDF_dmrs$isoforms, gffGenesDF_dmrs$geneSym %in% names(dmrsByGeneMulti))




# how many dmrs overlapped more than one gene?
sum(sapply(dmrsGR_geneOVhits, function(f) nrow(f)>1));
numOVgene = subset(numOV, gene>0);
numOVgeneSums = as.data.frame(cbind(multi=apply(subset(numOVgene, gene>1), 2, sum), 
                                    single=apply(subset(numOVgene, gene==1), 2, sum)));
for (i in c(2,3,8,10)) {
  print(rownames(numOVgeneSums)[i]);
  tmpmat = matrix(c(numOVgeneSums$multi[i],
                    numOVgeneSums$single[i],
                    sum(numOVgene$gene>1) - numOVgeneSums$multi[i],
                    sum(numOVgene$gene==1) - numOVgeneSums$single[i]),ncol=2);
  print(tmpmat)
  print(fisher.test(tmpmat))
}
rm(i,tmpmat);

# build info table about multi-gene DMRs


dmrsGRhitsGeneMulti0 = do.call('rbind',lapply(dmrsGRhits[names(dmrsGRhits) %in% dmrsints[numOV$gene==2]], function(f) f$ov$dmrsGR_gene))
dmrsGRhitsGeneMulti = as.data.frame(cbind(DMR=sapply(strsplit(rownames(dmrsGRhitsGeneMulti0), '.', fixed=T), function(f) f[1]),
                                          dmrsGRhitsGeneMulti0[,c(6,8)]));
dmrsGRhitsGeneMulti = as.data.frame(cbind(dmrsGRhitsGeneMulti, 
                                          desc=gffgenedesc[match(dmrsGRhitsGeneMulti$mcols.geneSym, gffgenesyms)]))
dmrsGRhitsGeneMulti$desc = gsub('%2C','',dmrsGRhitsGeneMulti$desc);
dmrsGRhitsGeneMulti$desc = sapply(strsplit(dmrsGRhitsGeneMulti$desc, ' variant'), function(f) f[1]);
dmrsGRhitsGeneMulti$desc = sapply(strsplit(dmrsGRhitsGeneMulti$desc, ' transcript'), function(f) f[1]);
dmrsGRhitsGeneMulti=as.data.frame(cbind(dmrsGRhitsGeneMulti, hypermethylated=dmrs$direction[match(dmrsGRhitsGeneMulti$DMR, dmrsints)]))
dmrsGRhitsGeneMulti$hypermethylated = gsub('hypo','D',dmrsGRhitsGeneMulti$hypermethylated)
dmrsGRhitsGeneMulti$hypermethylated = gsub('hyper','ND',dmrsGRhitsGeneMulti$hypermethylated)

# 7/15 gene pairs that share a dmr are a protein_coding gene paired with a lnc
#  is this significantly more than expected?
gene_gene_OV = as.data.frame(findOverlaps(geneGR,ignoreSelf=T,ignore.strand=T));
gene_gene_OVlist = apply(gene_gene_OV, 1, function(f) geneGR[unlist(f)]);
gene_gene_OVlist = gene_gene_OVlist[sapply(gene_gene_OVlist, 
                                           function(f) length(unique(as.character(strand(f)))))==2]
gene_gene_OVlist_biotypes = lapply(gene_gene_OVlist, function(f) mcols(f)[,3]);

ratios = c();
for (tr in 1:10000) {
  if (tr %% 1000 == 0) {cat(tr,'...')}
  rtmp = sapply(gene_gene_OVlist_biotypes[sample(1:length(gene_gene_OVlist_biotypes), 
                                                 nrow(dmrsGRhitsGeneMulti)/2)], 
                function(f) length(unique(f)));
  ratios = c(ratios, sum(rtmp==2) / sum(rtmp==1));
}; rm(tr,rtmp);





dmrsGRhitsGeneSingle0 = do.call('rbind',lapply(dmrsGRhits[names(dmrsGRhits) %in% dmrsints[numOV$gene==1]], function(f) f$ov$dmrsGR_gene))



# how many dmrs were in a gene basal regulatory region?
dmrsGR_brOVhits = dmrsGR_brOV[sapply(dmrsGR_brOV, function(f) !all(is.na(f)))];
# how many of these hits were to different genes than the overlap hit?
for (i in 1:length(dmrsGR_brOVhits)) {
  govind = which(names(dmrsGR_geneOVhits) == names(dmrsGR_brOVhits)[i]);
  intersect(dmrsGR_brOVhits[[i]]$mcols.geneSym, dmrsGR_geneOVhits[[govind]]$mcols.geneSym)
}; rm(i,govind);

# how many dmrs are within 20kb of a gene?
dmrsGR_geneNNhits = dmrsGR_geneNN[sapply(dmrsGR_geneNN, function(f) !all(is.na(f)))];

dmrsGR_geneNNhitsUp = lapply(dmrsGR_geneNNhits, function(f) subset(f, qstream=='up'));
dmrsGR_geneNNhitsUp = dmrsGR_geneNNhitsUp[sapply(dmrsGR_geneNNhitsUp, nrow) > 0];
dmrsGR_geneNNhitsUp = lapply(dmrsGR_geneNNhitsUp, function(f) f[f$dist==min(f$dist), ]);

dmrsGR_geneNNhits20kb = lapply(dmrsGR_geneNNhits, function(f) subset(f, dist<20000));
dmrsGR_geneNNhits20kb = dmrsGR_geneNNhits20kb[sapply(dmrsGR_geneNNhits20kb, nrow) > 0];

# 20kb upstream?
dmrsGR_geneNNhits20kbUp = lapply(dmrsGR_geneNNhits20kb, function(f) subset(f, any(qstream=='up')));
dmrsGR_geneNNhits20kbUp = dmrsGR_geneNNhits20kbUp[sapply(dmrsGR_geneNNhits20kbUp, nrow) > 0];

# how many dmrs are within 5kb of a gene?
dmrsGR_geneNNhits5kb = lapply(dmrsGR_geneNNhits, function(f) subset(f, dist<5000));
dmrsGR_geneNNhits5kb = dmrsGR_geneNNhits5kb[sapply(dmrsGR_geneNNhits5kb, nrow) > 0];

# 5kb upstream?
dmrsGR_geneNNhits5kbUp = lapply(dmrsGR_geneNNhits5kb, function(f) subset(f, any(qstream=='up')));
dmrsGR_geneNNhits5kbUp = dmrsGR_geneNNhits5kbUp[sapply(dmrsGR_geneNNhits5kbUp, nrow) > 0];

# draw venn diagram of DMRs overlapping a gene, DMRs within 20kb of a gene, DMRs within 20kb upstream
xxx=.draw_triple_venn_from_groups(names(dmrsGR_geneOVhits),
                                  names(dmrsGR_geneNNhits20kb),
                                  names(dmrsGR_geneNNhits20kbUp),
                                  groupNames=c('overlap','within 20kb', 'within 20kb up'))

# for dmrs that didn't overlap a gene, how close were their nearest neighbors?
dmrsGR_geneOVnohits = names(dmrsGR_geneOV[sapply(dmrsGR_geneOV, function(f) all(is.na(f)))])
dmrsGR_geneOVnohitsNN = dmrsGR_geneNN[names(dmrsGR_geneNN) %in% dmrsGR_geneOVnohits];
dmrsGR_geneOVnohitsNNhits = dmrsGR_geneOVnohitsNN[sapply(dmrsGR_geneOVnohitsNN, nrow) > 0];
summary(sapply(dmrsGR_geneOVnohitsNNhits, function(f) min(f$dist,na.rm=TRUE)));

dmrsGR_geneOVnohitsNNhits20kb = lapply(dmrsGR_geneOVnohitsNNhits, function(f) subset(f, dist<20000));
dmrsGR_geneOVnohitsNNhitsOutside20kb = names(dmrsGR_geneOVnohitsNNhits20kb)[which(sapply(dmrsGR_geneOVnohitsNNhits20kb, nrow) == 0)]
#dmrsGR_geneOVnohitsNNhits20kbUp = dmrsGR_geneOVnohitsNNhits[sapply(dmrsGR_geneOVnohitsNNhits, function(f) any(f$dist<20000 & f$qstream=='up', na.rm=T))];

dmrsGR_geneOVnohitsNNnohits = names(dmrsGR_geneOVnohitsNN[sapply(dmrsGR_geneOVnohitsNN, function(f) all(is.na(f)))])







# what types of TEs did DMRs overlap?
tenames0 = unlist(sapply(dmrsGR_teOVhits, function(f) f$mcols))
tenames = sapply(strsplit(tenames0, '_'), function(f) f[1])
tecounts=table(tenames);
pie(sort(tecounts), labels=paste0(names(sort(tecounts)),':',sort(tecounts)));

# how close were the TEs to genes?  not working...
teGR2 = teGR[mcols(teGR)[,1] %in% tenames0];
teGR_geneGR_OV = .getOverlapsStranded(teGR2, geneGR)
tenamesGene = names(teGR_geneGR_OV[!is.na(teGR_geneGR_OV)]);

teGR_geneGR_NN = .getNearestNeighborsStrandedNoOverlap(teGR2[!(mcols(teGR2)[,1] %in% tenamesGene)], geneGR)

brGR_teGR_OV = .getOverlapsStranded(brGR, teGR)



################

dmrsGR20kbhitsOV = .combineOvByDmr(paste0('dmrsGR20kb_',features));



numOV20kb = .make_numOV(dmrsGR20kbhitsOV, features, 'dmrsGR20kb_');


nogene20kb = rownames(subset(numOV20kb, gene==0));
nogene20kbscaffolds = unique(sapply(strsplit(nogene20kb, ':'), function(f) f[1]));










# check nn distances to genes for dmrs that didn't directly overlap a gene
nogeneOV_NN = dmrsGR_geneNN[numOV$gene==0];
nogeneOV_NN20kb = .subsetNearestNeighborHitListOnDistance(nogeneOV_NN, up=20000, down=20000);
nogeneOV_NN20kb_numhits = sapply(nogeneOV_NN20kb, nrow);
nogeneOV_NN20kb_numhits[which(sapply(nogeneOV_NN20kb_numhits, is.null))] = 0;
nogeneOV_NN20kb_numhits = unlist(nogeneOV_NN20kb_numhits);
nogeneOV_NN20kb0 = names(nogeneOV_NN20kb[nogeneOV_NN20kb_numhits == 0]);
nogeneOV_NN20kb = nogeneOV_NN20kb[nogeneOV_NN20kb_numhits > 0];

yesgeneOV_NN = dmrsGR_geneNN[numOV$gene > 0];
yesgeneOV_NN20kb = .subsetNearestNeighborHitListOnDistance(yesgeneOV_NN, up=20000, down=20000);
yesgeneOV_NN20kb_numhits = sapply(yesgeneOV_NN20kb, nrow);
yesgeneOV_NN20kb_numhits[which(sapply(yesgeneOV_NN20kb_numhits, is.null))] = 0;
yesgeneOV_NN20kb_numhits = unlist(yesgeneOV_NN20kb_numhits);
yesgeneOV_NN20kb0 = names(yesgeneOV_NN20kb[yesgeneOV_NN20kb_numhits == 0]);
yesgeneOV_NN20kb = yesgeneOV_NN20kb[yesgeneOV_NN20kb_numhits > 0];

numOVlog = apply(numOV, 2, as.logical);
TOPLOT = as.data.frame(numOVlog[,-c(4,5,9)]); rownames(TOPLOT) = rownames(numOV)
TOPLOT = as.data.frame(cbind(hypo=rownames(TOPLOT) %in% groups$hypo, TOPLOT))
ORDER = order(TOPLOT$hypo, TOPLOT$gene, TOPLOT$br, TOPLOT$exon, TOPLOT$cds, TOPLOT$te, TOPLOT$z, TOPLOT$tss, decreasing=T);
TOPLOT = TOPLOT[ORDER, match(c('hypo','gene','br','exon','cds','te','z','tss'), names(TOPLOT))];
TOPLOT = apply(TOPLOT, 2, function(f) as.logical(gsub(' ','',f)));
image(t(as.matrix(TOPLOT[nrow(TOPLOT):1, ])), col=c('white','greenyellow'), axes=F); 
axis(side=1, labels=colnames(TOPLOT), at=seq(0,1,(1/(ncol(TOPLOT)-1))), las=3);



# dmrsGRhitsOV_gene = do.call('rbind',lapply(dmrsGRhitsOV, function(f) f$dmrsGR_gene));
# dmrsGRhitsOV_gene = dmrsGRhitsOV_gene[!apply(dmrsGRhitsOV_gene, 1, function(f) all(is.na(f))), ];
# 
# dmrsGRhitsOV_tss = do.call('rbind',lapply(dmrsGRhitsOV, function(f) f$dmrsGR_tss));
# dmrsGRhitsOV_tss = dmrsGRhitsOV_tss[!apply(dmrsGRhitsOV_tss, 1, function(f) all(is.na(f))), ];

dmrsGRhitsOV_gene0 = lapply(dmrsGRhitsOV, function(f) f$dmrsGR_gene);
dmrsGRhitsOV_gene = dmrsGRhitsOV_gene0[sapply(dmrsGRhitsOV_gene0, function(f) !all(is.na(f)))];

dmrsGRhitsOV_tss0 = lapply(dmrsGRhitsOV, function(f) f$dmrsGR_tss);
dmrsGRhitsOV_tss = dmrsGRhitsOV_tss0[sapply(dmrsGRhitsOV_tss0, function(f) !all(is.na(f)))];


# # 
# # #
# # 
# #
# ovDF = do.call('rbind', lapply(dmrsGRhits, 
#                                function(f) sapply(f$ov, is.data.frame)));
# ovDF = as.data.frame(cbind(hypo=dmrs$direction=='hypo', ovDF));
# 
# toplot = ovDF;
# colOrder = c(1, (order(-apply(toplot[,-1],2,sum))+1));
# toplot = toplot[, colOrder];
# plotOrder = eval(parse(text=paste0('order(', 
#                                    paste0(paste0('!toplot$',names(toplot)), collapse=','), 
#                                    ')')));
# toplot = toplot[plotOrder, ];
# image(t(as.matrix(toplot[nrow(toplot):1, ])), col=c('white','greenyellow'), axes=F);
# axis(side=1, labels=colnames(toplot), at=seq(0,1,(1/(ncol(toplot)-1))), las=3);
# for (i in seq(0, 1, 1/nrow(toplot))){abline(h=i, lwd=.25, col='grey1')}



.buildDmrFeatureHitCountNullDist = function (dmrs, scaffold, sl, geneGR, dmrsGR_geneOV, trials, plot=T) {
  dmrswidths = round(subset(dmrs, chr==scaffold)$width);
  scaffoldlength = sl[names(sl)==scaffold];
  realnum = sum(!sapply(dmrsGR_geneOV[grepl(paste0(scaffold,':'),
                                            names(dmrsGR_geneOV),fixed=T)], 
                        function(f) all(is.na(f))));
  nhits = c();
  for (tr in 1:trials) {
    if (tr %% 100 == 0) { cat(tr,'...') }
    starts = sample(1:(scaffoldlength-max(dmrswidths)), length(dmrswidths));
    thisgr = GRanges(seqnames=scaffold, ranges=IRanges(start=starts, end=(starts+dmrswidths-1)));
    ov = as.data.frame(findOverlaps(thisgr, geneGR[seqnames(geneGR) == scaffold]));
    nhits = c(nhits, length(unique(ov$queryHits)));
  }
  pval = sum(nhits >= realnum) / trials;
  if (plot) {
    hist(nhits, col='grey', border='darkgrey', main=paste0(pval), xlab=''); 
    abline(v=realnum, col='red');
  }
  return(list(realhits=realnum, pval=pval, nullhits=nhits));
}

dmrscaffolds = names(table(dmrs$chr));
nlist = list();
for (scaf in dmrscaffolds) {
  cat(scaf,'...');
  nlist[[scaf]] = .buildDmrFeatureHitCountNullDist(dmrs, scaf, sl, geneGR, dmrsGR_geneOV, 1000, FALSE);
}; rm(scaf);


.buildDmrFeatureHitCountNullDistWholeGenome = function (dmrs, sl, geneGR, dmrsGR_geneOV, trials, plot) {
  dmrsWidthsByScaffold = sapply(split(dmrs, dmrs$chr), function(f) round(f$width))
  
  for (i in 1:trials) {
    # generate intervals for each scaffold
    for (ii in 1:length(dmrsWidthsByScaffold)) 
    
  }
  
}



# ==========================================================================================================
# ==========================================================================================================
# get GO results for gene-based features
# ==========================================================================================================
# ==========================================================================================================
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# get different gene lists
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
#
bg = unique(unlist(strsplit(annoCombo$new$hsaHomologEntrez, ',')));
org.Hs.egGO2ALLEGSlist = as.list(org.Hs.egGO2ALLEGS);
limit=2; limitMethod='first.sorted';
suffixes = c('NN','NN.hypo','NN.hyper','OV','OV.hypo','OV.hyper');

# genes
feat = 'dmrsGR_gene';

# overlap anywhere
for (nnov in grep('^OV',suffixes,value=T)) {
  inname = paste0(feat,nnov);
  cat(inname,'...');
  assign(paste0(inname,'.ids'), 
         .buildIdsFromOV_NN(fOV_NN=get(inname), featureName='mcols.geneSym', 
                            maplist=abHsMap, limit=limit, limitMethod=limitMethod, smartUnlist=TRUE));
}; rm(nnov,inname);

feat = 'dmrsGR20kb_gene';
for (nnov in grep('^OV',suffixes,value=T)) {
  inname = paste0(feat,nnov);print(inname)
  cat(inname,'...');
  assign(paste0(inname,'.ids'), 
         .buildIdsFromOV_NN(fOV_NN=get(inname), featureName='mcols.geneSym', 
                            maplist=abHsMap, limit=limit, limitMethod=limitMethod, smartUnlist=TRUE));
}; rm(nnov,inname);



# basal regulatory regions
feat='dmrsGR_br'
for (nnov in grep('^OV',suffixes,value=T)) {
  inname = paste0(feat,nnov);
  cat(inname,'...');
  assign(paste0(inname,'.ids'), 
         .buildIdsFromOV_NN(fOV_NN=get(inname), featureName='mcols.geneSym', 
                            maplist=abHsMap, limit=limit, limitMethod=limitMethod, smartUnlist=TRUE));
}; rm(nnov,inname);

# within 20kb 
feat = 'dmrsGR_gene';
for (nnov in grep('^NN',suffixes,value=T)) {
  inname = paste0(feat,nnov);
  cat(inname,'...');
  assign(paste0(inname,'.20kbUpDown'), 
         .subsetNearestNeighborHitListOnDistance(dmrsNN=get(inname), up=20000, down=20000));
  assign(paste0(inname,'.20kbUpDown.ids'), 
         .buildIdsFromOV_NN(fOV_NN=get(paste0(inname,'.20kbUpDown')), featureName='name', 
                            maplist=abHsMap, limit=limit, limitMethod=limitMethod, smartUnlist=TRUE));
}; rm(nnov,inname);

# within 50kb 
feat = 'dmrsGR_gene';
for (nnov in grep('^NN',suffixes,value=T)) {
  inname = paste0(feat,nnov);
  cat(inname,'...');
  assign(paste0(inname,'.50kbUpDown'), 
         .subsetNearestNeighborHitListOnDistance(dmrsNN=get(inname), up=50000, down=50000));
  assign(paste0(inname,'.50kbUpDown.ids'), 
         .buildIdsFromOV_NN(fOV_NN=get(paste0(inname,'.50kbUpDown')), featureName='name', 
                            maplist=abHsMap, limit=limit, limitMethod=limitMethod, smartUnlist=TRUE));
}; rm(nnov,inname);




# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# compute GO and KEGG enrichments
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
#
# IDS = list(gene.ids=c(dmrsGR_geneOV.ids$ids$newvec, dmrsGR_geneNN.50kbUpDown.ids$ids$newvec),
#            gene.hypo.ids=c(dmrsGR_geneOV.hypo.ids$ids$newvec, dmrsGR_geneNN.hypo.50kbUpDown.ids$ids$newvec),
#            gene.hyper.ids=c(dmrsGR_geneOV.hyper.ids$ids$newvec, dmrsGR_geneNN.hyper.50kbUpDown.ids$ids$newvec));
# IDS = list(geneOV.ids=dmrsGR_geneOV.ids$ids$newvec,
#            geneOV.hypo.ids=dmrsGR_geneOV.hypo.ids$ids$newvec,
#            geneOV.hyper.ids=dmrsGR_geneOV.hyper.ids$ids$newvec);
# IDS = list(geneNN.ids=dmrsGR_geneNN.20kbUpDown.ids$ids$newvec,
#            geneNN.hypo.ids=dmrsGR_geneNN.hypo.20kbUpDown.ids$ids$newvec,
#            geneNN.hyper.ids=dmrsGR_geneNN.hyper.20kbUpDown.ids$ids$newvec);
# IDS = list(brOV.ids=dmrsGR_brOV.ids$ids$newvec,
#            brOV.hypo.ids=dmrsGR_brOV.hypo.ids$ids$newvec,
#            brOV.hyper.ids=dmrsGR_brOV.hyper.ids$ids$newvec);
# IDS = list(brOV_geneNN.ids=c(dmrsGR_brOV.ids$ids$newvec, dmrsGR_geneNN.20kbUpDown.ids$ids$newvec),
#            brOV_geneNN.hypo.ids=c(dmrsGR_brOV.hypo.ids$ids$newvec, dmrsGR_geneNN.hypo.20kbUpDown.ids$ids$newvec),
#            brOV_geneNN.hyper.ids=c(dmrsGR_brOV.hyper.ids$ids$newvec, dmrsGR_geneNN.hyper.20kbUpDown.ids$ids$newvec));
IDS = list(gene_br.ids=unique(c(dmrsGR_geneOV.ids$ids$newvec, dmrsGR_brOV.ids$ids$newvec)),
           gene_br.hypo.ids=unique(c(dmrsGR_geneOV.hypo.ids$ids$newvec, dmrsGR_brOV.hypo.ids$ids$newvec)),
           gene_br.hyper.ids=unique(c(dmrsGR_geneOV.hyper.ids$ids$newvec, dmrsGR_brOV.hyper.ids$ids$newvec)));
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

for (i in grep('go.geneOV',ls(),fixed=T,value=T) ) { print(i);print(get(i)[-8]) }; rm(i)

# kegg
PVALCUT = 0.1;
ADJMETH = 'fdr';
UNI = bg;
MINGS = 10;
kegg.all0 = summary(enrichKEGG(c(dmrsGR_brOV.ids$ids$newvec, dmrsGR_geneNN.20kbUpDown.ids$ids$newvec), universe=UNI, readable=T, 
                               pvalueCutoff=PVALCUT, minGSSize=MINGS, pAdjustMethod=ADJMETH)); 
kegg.hypo0 = summary(enrichKEGG(ids0$hypo, universe=UNI, readable=T, 
                                pvalueCutoff=PVALCUT, minGSSize=MINGS, pAdjustMethod=ADJMETH)); 
kegg.hyper0 = summary(enrichKEGG(ids0$hyper, universe=UNI, readable=T, 
                                 pvalueCutoff=PVALCUT, minGSSize=MINGS, pAdjustMethod=ADJMETH)); 


# get dmrs for a term
thisgo = go.brOV_geneNN.ids_BH.1;
thisgoinput = c(dmrsGR_brOV.ids$ids$newvec, dmrsGR_geneNN.20kbUpDown.ids$ids$newvec);

termidsHs = .getGOFunctionTermGenes(thisgo, row=2, geneCol=8);
termidsAb = thisgoinput[thisgoinput %in% termidsHs];
tmp = c(dmrsGR_brOV.ids$processedInput$fvec[dmrsGR_brOV.ids$processedInput$fvec %in% names(termidsAb)],
        dmrsGR_geneNN.20kbUpDown.ids$processedInput$fvec[dmrsGR_geneNN.20kbUpDown.ids$processedInput$fvec %in% names(termidsAb)]);
goodnames = names(tmp) %in% dmrsints;
goodints0 = names(tmp)[goodnames];
badints = names(tmp)[!goodnames];
goodints = c(goodints0, substr(badints, 1, nchar(badints)-1));


tmpall = c(dmrsGR_brOV.ids$processedInput$fvec[dmrsGR_brOV.ids$processedInput$fvec %in% names(thisgoinput)],
           dmrsGR_geneNN.20kbUpDown.ids$processedInput$fvec[dmrsGR_geneNN.20kbUpDown.ids$processedInput$fvec %in% names(thisgoinput)]);

goodnamesall = names(tmpall) %in% dmrsints;
goodints0all = names(tmpall)[goodnamesall];
badintsall = names(tmpall)[!goodnamesall];
goodintsall = c(goodints0all, substr(badintsall, 1, nchar(badintsall)-1));



# check for GO enrichments within term genes
.GOFunctionAllOntologies(termidsHs, bg, fdrmethod='BH', fdrth=.1);

# check overlaps and neighbors for dmrs associated with term genes
termgeneHits = dmrsGRhits[names(dmrsGRhits) %in% goodints];
termgeneTfHits0 = tf$hitsWithCpGs[names(tf$hitsWithCpGs) %in% goodints];
termgeneTfHits = lapply(termgeneTfHits0, function(f) unique(f$TF));

allgeneHits = dmrsGRhits[names(dmrsGRhits) %in% goodintsall];
allgeneTfHits0 = tf$hitsWithCpGs[names(tf$hitsWithCpGs) %in% goodintsall];
allgeneTfHits = lapply(allgeneTfHits0, function(f) unique(f$TF));

# atf = sort(table(unlist(tf$hitsWithCpGs)[grep('.TF', names(unlist(tf$hitsWithCpGs)), fixed=T)]), decreasing=T);
atf = sort(table(unlist(allgeneTfHits)), decreasing=T);
ttf = sort(table(unlist(termgeneTfHits)), decreasing=T);



for (tt in 1:length(ttf)) {
  #print(names(ttf)[tt]);
  mat = matrix(c(ttf[tt], 
                 length(termgeneHits)-ttf[tt], 
                 atf[match(names(ttf)[tt], names(atf))], 
                 length(allgeneHits) - atf[match(names(ttf)[tt], names(atf))]), 
               ncol=2)
  ft = fisher.test(mat)
  if (ft$p.value < .05) { print(names(ttf)[tt]);print(ft); print(mat);print('------------------------------------------') }
}
rm(tt,ft,mat)












testgenes = geneOV$fvec;
hsa = useMart('ENSEMBL_MART_ENSEMBL','hsapiens_gene_ensembl',host='www.ensembl.org');
getBM(attributes=c('go_id'), filters='entrezgene', values=testgenesHs[1:10], mart=hsa);






# .convertIDsWithBiomaRt = function (biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', fromID='ensembl_gene_id', toID='entrezgene', ids) {
#   mart = useMart(biomart=biomart, dataset=dataset, host='www.ensembl.org');
#   return(getBM(attributes=c(fromID, toID), filters=fromID, values=ids, mart=mart));
# }

# hsa = useMart('ENSEMBL_MART_ENSEMBL','hsapiens_gene_ensembl',host='www.ensembl.org');
# library(GO.db);
# GOterms = as.list(GOTERM);

# getModGOfromEntrez = function(entrezIDsList) {
# GOlist = list();
# for (m in 1:length(entrezIDsList)) {
# print(names(entrezIDsList)[m])
# GOlist[[m]] = getBM(attributes=c('go_id'), filters='entrezgene', values=entrezIDsList[[m]], mart=hsa);
# names(GOlist)[m] = names(entrezIDsList)[m];
# }

# return(GOlist);
# }

# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
# parse GOFunction results
# ----------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------
#


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

#
tfidlookup = .buildTfIdNameLookup(tfhitsDmrsMulti);

#
tf = .analyzeTfHitsByDmr(tfhitsDmrsMulti, cpgcommon);

#

ar_hits = do.call('rbind',lapply(tf$hitsWithCpGs, function(f) f[f$TF=='AR',]))
ar_hits_numOV = numOV[match(ar_hits$seqnames, rownames(numOV)), ];

.verboseScatterplotAllColumnPairs(tf$stats[,4:7],mfrow=c(2,3),col=labels2colors(tf$stats$class),pch=20,cex=3);
.verboseBoxplotColumns(tf$stats[4:7],tf$stats$class,c(4,1),col=names(table(labels2colors(tf$stats$class))));
.verboseScatterplotVecAgainstColumns(abs(dmrs[,7:15]), tf$maxCpgByHit, col=labels2colors(dmrs$direction),cex=2,pch=20);
.verboseBoxplotColumns((dmrs[,7:15]), as.factor(tf$maxCpgByHit), col='grey',frame.plot=F);


ind = 1
ccci = ccc[[ind]]
maxhits = ccci$ov$hits[ccci$ov$numHits == max(ccci$ov$numHits)]
thismaxhit = maxhits[[1]];
thismaxhitint = unlist(strsplit(strsplit(names(maxhits), ':')[[1]], '-'));
tmp[[ind]][tmp[[ind]]$start==thismaxhitint[2] & tmp[[ind]]$end==thismaxhitint[3], ]

thisdmrpos = thismaxhit$start
thisabspos = ccci$dmrloci$abspos[ccci$dmrloci$dmrpos %in% thismaxhit$start]





# most numerous top hits?
numtop = sort(table(unlist(sapply(tmp, function(f) f$ID[1]))), decreasing=T);
tfid = names(numtop)[1];
attributes(tfhitsDmrsMulti[[tfid]]$tfm[[1]])$name
attributes(tfhitsDmrsMulti[[tfid]]$tfm[[1]])$tag$family




tmp1 = tmp[[1]]
tmp1loci = .countLociInTfHit(tmp1,'scaffold_1314',57386,58160,cpgcommon)
tmpov = .getTfHitOverlaps(tmp1);
tmpovmat = .buildTfHitOverlapTable(tmpov)
i=apply(tmpovmat, 1, function(f) sum(f)>0);   heatmap(cor(tmpovmat[i, i]), symm=T); rm(i)

####
patternTtest = 'mg1e+08_fits_4xCov_Ttests.RData';
filesTtest = grep(pattern, list.files(), value=T, fixed=T);
ttests = list();
for (i in filesTtest) {
  cat(i,'...')
  ttests[[i]] = load(i);
  ttests[[i]] = get(ttests[[i]][[1]]);
  ttests[[i]] = ttests[[i]]$tstats;
}
rm(i);
names(ttests) = gsub('CpGcombined.CG_4xCov_', '', names(ttests), fixed=T);
names(ttests) = gsub('mg1e+08_fits_4xCov_Ttests.RData', '', names(ttests), fixed=T);


# ==========================================================================================================
# ==========================================================================================================
# get gene expression information for dmr-related genes
# ==========================================================================================================
# ==========================================================================================================
exprdir = '/Volumes/fishstudies-1/_LYNLEY_RNAseq/data/_kallisto_092915/';
nd1exp = read.table(paste0(exprdir, 'ATCACG/abundance.tsv'), header=T, sep='\t');
rownames(nd1exp) = sapply(strsplit(nd1exp$target_id, '|', fixed=T), function(f) f[length(f)]);
nd2exp = read.table(paste0(exprdir, 'TGACCA/abundance.tsv'), header=T, sep='\t');
rownames(nd2exp) = sapply(strsplit(nd2exp$target_id, '|', fixed=T), function(f) f[length(f)]);
d1exp = read.table(paste0(exprdir, 'CGATGT/abundance.tsv'), header=T, sep='\t');
rownames(d1exp) = sapply(strsplit(d1exp$target_id, '|', fixed=T), function(f) f[length(f)]);
d2exp = read.table(paste0(exprdir, 'TTAGGC/abundance.tsv'), header=T, sep='\t');
rownames(d2exp) = sapply(strsplit(d2exp$target_id, '|', fixed=T), function(f) f[length(f)]);

tpm = as.data.frame(cbind('3157_TENNISON'=nd1exp$tpm, 
                          '3165_BRISCOE'=d1exp$tpm, 
                          '3581_LYNLEY'=d2exp$tpm, 
                          '3677_MONK'=nd2exp$tpm));
rownames(tpm) = rownames(nd1exp);
rownames(tpm) = paste(rownames(tpm), 
                      an$lookup$gene[match(rownames(nd1exp), an$lookup$transcript_id)], 
                      sep=':');
tpmGeneNames = sapply(strsplit(rownames(tpm), ':'), function(f) f[2]);
tpmGeneList = split(tpm, tpmGeneNames);
tpmGeneMat = do.call('rbind', lapply(tpmGeneList, function(f) apply(f,2,sum)));
#tpmGeneMat = tpmGeneMat[apply(tpmGeneMat, 1, sum) != 0, ];
tpmGeneDfQnorm = as.data.frame(normalize.quantiles(tpmGeneMat))
dimnames(tpmGeneDfQnorm) = dimnames(tpmGeneMat);
tpmD = apply(tpmGeneDfQnorm[, which(names(tpmGeneDfQnorm) %in% c('3165_BRISCOE','3581_LYNLEY'))], 
             1, mean);
tpmND = apply(tpmGeneDfQnorm[, which(names(tpmGeneDfQnorm) %in% c('3157_TENNISON','3677_MONK'))], 
             1, mean);
#tpmFc = as.numeric(log2( (tpmD+ (1e-200)) / (tpmND+ (1e-200)) ));
#tpmFc = as.numeric(log2( (tpmD) / (tpmND) ));
names(tpmFc) = rownames(tpmGeneDfQnorm);




tpmFc_dmrGenes = tpmFc[match(dmrsGR_geneOVhitsGenes, names(tpmFc))]


dmrFc = as.numeric(log2( dmrs$group2.mean / dmrs$group1.mean ));
names(dmrFc) = dmrsints

#ovhits = dmrsGR_geneOVhits
.buildFcVecs = function (ovhits, tpmFc, dmrFc) {
  dmrFcVec = c();
  tpmFcVec = c();
  for (i in 1:length(ovhits)) {
    cdscheck = ovhits[[i]]$mcols.biotype=='protein_coding';
    if (any(cdscheck)) {
      thisdmr = names(ovhits)[i];
      genes = ovhits[[i]]$mcols.geneSym[cdscheck];
      tpmmatch = match(genes, names(tpmFc))
      tpmFcVec = c(tpmFcVec, tpmFc[tpmmatch]);
      dmrFcVec = c(dmrFcVec, rep(dmrFc[names(dmrFc)==thisdmr], length(tpmmatch)));
    } else {
      tpmFcVec = c(tpmFcVec, NA);
      dmrFcVec = c(dmrFcVec, rep(dmrFc[names(dmrFc)==thisdmr], length(tpmmatch)));
    }
  }
 #return(as.data.frame(cbind(tpm=tpmFcVec, dmr=dmrFcVec)))
 return(list(tpm=tpmFcVec, dmr=dmrFcVec))
}

dmrsGR_brOVhits_fc = .buildFcVecs(dmrsGR_brOVhits,tpmFc,dmrFc);
dmrsGR_brOVhits_fc = as.data.frame(cbind(tpm=dmrsGR_brOVhits_fc$tpm, dmr=dmrsGR_brOVhits_fc$dmr))

dmrsGR_geneOVhits_fc = .buildFcVecs(dmrsGR_geneOVhits,tpmFc,dmrFc);
dmrsGR_geneOVhits_fc = as.data.frame(cbind(tpm=dmrsGR_geneOVhits_fc$tpm, dmr=dmrsGR_geneOVhits_fc$dmr))


signmatchcheck=(cbind(neg=apply(dmrsGR_geneOVhits_fc, 1, function(f) all(f<0)), 
                      pos=apply(dmrsGR_geneOVhits_fc, 1, function(f) all(f>0))))


dmrsGR_geneOVhits_fcAnti = dmrsGR_geneOVhits_fc[apply(signmatchcheck, 1, sum)==0, ]
dmrsGR_geneOVhits_fcSame = dmrsGR_geneOVhits_fc[apply(signmatchcheck, 1, sum)==1, ]




## make bed files for expression
tpmGeneDfQnormAnno = as.data.frame(cbind(an$gffGenesDF, tpmGeneDfQnorm[match(an$gffGenesDF$geneSym, rownames(tpmGeneDfQnorm)), ]))

for (i in 11:14) {
  write.table(as.data.frame(cbind(as.character(tpmGeneDfQnormAnno$seqnames), 
                                  tpmGeneDfQnormAnno$start, tpmGeneDfQnormAnno$end,
                                  tpmGeneDfQnormAnno[,i])), 
              file = paste0('tpm_',names(tpmGeneDfQnormAnno)[i],'.bedGraph'),
              quote=F, row.names=F, col.names=F, sep='\t');
}; rm(i)




write.table(tpmGeneDfQnormAnno$seqnames, 
            tpmGeneDfQnormAnno$start, tpmGeneDfQnormAnno$end,
            tpmGeneDfQnormAnno, quote=F, row.names=F, col.names=F, sep='\t');

#tpm = tpm[apply(tpm, 1, sum) > 0, ];


# xxx = subset(hc2, status=='ok')
# xxx = .addTPM(xxx, tpmqn, which(colnames(xxx)=='gene'));
# 
# dmrs = xxx;
# dmrs = cbind(dmrs, 
#              genefc=as.numeric(log2( (dmrs$expD+1)  / (dmrs$expND+1)  )),
#              dmrfc=as.numeric(log2( (dmrs$group2.mean)  / (dmrs$group1.mean)  ))
# );




# DESeq2
counts = as.data.frame(cbind('3157_TENNISON'=nd1exp$est_counts, 
                          '3165_BRISCOE'=d1exp$est_counts, 
                          '3581_LYNLEY'=d2exp$est_counts, 
                          '3677_MONK'=nd2exp$est_counts));

rownames(counts) = rownames(nd1exp);
counts = counts[match(an$lookup$transcript_id,rownames(counts)), ]
countsl = split(counts,an$lookup$gene)

counts = as.data.frame(t(sapply(countsl, function(f) apply(f, 2, sum))));
coldata = as.data.frame(matrix(c('ND','D','D','ND'), ncol=1));
rownames(coldata) = names(counts);
names(coldata) = 'group'


dds = DESeqDataSetFromMatrix(countData=round(counts), 
                             colData=coldata, 
                             design = ~ group
);
dds = DESeq(dds, parallel=TRUE);

res = results(dds, contrast=c('group','D','ND'));

deseq_dmrGenes = res[match(allDmrGenes, rownames(res)), ]
deseq_dmrGenesFc = deseq_dmrGenes$log2FoldChange
names(deseq_dmrGenesFc) = rownames(deseq_dmrGenes);

allDmrGenesFc = as.data.frame(cbind(tpm=deseq_dmrGenesFc, dmr=rep(NA,length(deseq_dmrGenesFc))))
tmp = c(dmrsGR_geneOVhitsGenesAll, dmrsGR_brOVhitsGenesAll, dmrsGR_br3primeOVhitsGenesAll);
for (g in 1:nrow(allDmrGenesFc)) {
  thisgene = rownames(allDmrGenesFc)[g];
  for (d in 1:length(tmp)) {
    if (thisgene %in% tmp[[d]]) {
      thisdmr = names(tmp)[d];
      allDmrGenesFc$dmr[g] = dmrFc[names(dmrFc)==thisdmr];
      next;
    }
  }
}
rm(g,thisgene,d,thisdmr,tmp);




signmatchcheckDESeq=(cbind(neg=apply(allDmrGenesFc, 1, function(f) all(f<0)), 
                           pos=apply(allDmrGenesFc, 1, function(f) all(f>0))))

allDmrGenesFcAnti = allDmrGenesFc[apply(signmatchcheckDESeq, 1, sum)==0, ];
allDmrGenesFcAnti = allDmrGenesFcAnti[!apply(allDmrGenesFcAnti, 1, function(f) all(is.na(f))), ]
allDmrGenesFcSame = allDmrGenesFc[apply(signmatchcheckDESeq, 1, sum)==1, ]
allDmrGenesFcSame = allDmrGenesFcSame[!apply(allDmrGenesFcSame, 1, function(f) all(is.na(f))), ]

par(mfrow=c(1,3));
xlim = c(-3.5,4.15);
ylim=c(-.725,.9)
WGCNA::verboseScatterplot(dmrsGR_geneOVhits_fcDESeq$dmr, 
                          dmrsGR_geneOVhits_fcDESeq$tpm, 
                          abline=T, abline.col='red', frame.plot=F, 
                        #  pch=21, 
                        #  bg=WGCNA::labels2colors(rownames(dmrsGR_geneOVhits_fcDESeq) %in% dmrsGR_geneOVhits.hypoGenes), 
                          xlab='fc.meth', ylab='fc.expr', main=paste0('n=',nrow(dmrsGR_geneOVhits_fcDESeq),'\n'),
                          corOptions='method="p", use="p"',
                          xlim=xlim, ylim=ylim); 
abline(h=0,v=0,col='grey')

WGCNA::verboseScatterplot(dmrsGR_geneOVhits_fcDESeqAnti$dmr, 
                          dmrsGR_geneOVhits_fcDESeqAnti$tpm, 
                          abline=T, abline.col='red', frame.plot=F, 
                        #  pch=21, 
                        #  bg=WGCNA::labels2colors(rownames(dmrsGR_geneOVhits_fcDESeqAnti) %in% dmrsGR_geneOVhits.hypoGenes), 
                          xlab='fc.meth', ylab='fc.expr', main=paste0('n=',nrow(dmrsGR_geneOVhits_fcDESeqAnti),'\n'),
                          corOptions='method="p", use="p"',
                          xlim=xlim, ylim=ylim); 
abline(h=0,v=0,col='grey')

WGCNA::verboseScatterplot(dmrsGR_geneOVhits_fcDESeqSame$dmr, 
                          dmrsGR_geneOVhits_fcDESeqSame$tpm, 
                          abline=T, abline.col='red', frame.plot=F, 
                        #  pch=21, 
                        #  bg=WGCNA::labels2colors(rownames(dmrsGR_geneOVhits_fcDESeqSame) %in% dmrsGR_geneOVhits.hypoGenes), 
                          xlab='fc.meth', ylab='fc.expr',main=paste0('n=',nrow(dmrsGR_geneOVhits_fcDESeqSame),'\n'),
                          corOptions='method="p", use="p"',
                          xlim=xlim, ylim=ylim); 
abline(h=0,v=0,col='grey')









.addGenesToGOFunctionResults(xxx,org.Hs.egGO2ALLEGSlist,.mapAbGenesHsEntrez(rownames(allDmrGenesFcAnti), abHsMap,'ab'))



# dmrOVgenes = unique(na.omit(unlist(lapply(dmrsGRhitsOV_gene, function(f) f[6]))));
# dmrOVgenes.hypo = unique(na.omit(unlist(lapply(dmrsGRhitsOV_gene[names(dmrsGRhitsOV_gene) %in% groups$hypo], 
#                                                function(f) f[6]))));
# dmrOVgenes.hyper = unique(na.omit(unlist(lapply(dmrsGRhitsOV_gene[names(dmrsGRhitsOV_gene) %in% groups$hyper], 
#                                                function(f) f[6]))));
